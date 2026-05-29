import pickle 
import matplotlib.pyplot as plt
import spine_calculations as calc
import utilities as utils
import pandas as pd
import numpy as np
import os 
import re

class Dendrite_Imaging:
    def __init__(self,analysis_dir):
        """
        Initialize Dendrite_Imaging analysis object.

        Parameters
        ----------
        analysis_dir : str
            Path to root analysis directory. Expected structure:
            analysis_dir/<mouse_id>/<fov_id>/<date>/<pickle_file>
        """
        self.analysis_dir = analysis_dir

    JL117_DAYS = ['250211', '250212', '250213', '250214', '250215', '250216', '250217', '250218', '250219', '250220', '250221', '250222', '250224', '250225']
    JL122_DAYS = ['250402', '250403', '250404', '250405', '250406', '250407', '250408', '250409', '250410', '250411', '250412', '250413', '250414', '250415']
    JL123_DAYS = ['250402', '250403', '250404', '250405', '250406', '250407', '250408', '250409', '250410', '250411', '250412', '250413', '250414', '250415']
    
    def load_grouped_data(self,reanalyze=False,save_grouped=True):
        """
        Load & aggregate spine data across mice, FOVs, & days. 

        Loads per-day pickle files, aligns dates to each mouse's full day list, 
        inserts "missing" as placeholder for missing days, interpolates missing days, 
        builds MultiIndex DataFrame reindexed to all (spine_ID, day) combinations. 
        Computes: denrite_ID, dendrite length, and a binary "exists" flag.

        Parameters
        ----------
        reanalyze : bool, optional
            If True, reprocess from raw files.
            Defaults to False.
        save_grouped : bool, optional
            If True, save resulting DataFrame to 
            analysis_dir/grouped/grouped_data.pickle. Defaults to True.

        Returns
        -------
        self.spine_df : pd.DataFrame
            Columns : spine_ID, day, dendrite_ID, layer, position, volume, flags, exists, dend_length
        """
        analysis_dir = self.analysis_dir
        
        # Path to the grouped data 
        savepath=os.path.join(analysis_dir,"grouped","grouped_data.pickle")
        file_exist = os.path.exists(savepath)
        if reanalyze or not file_exist:
            df_list = []

            # Loop across mice
            for mouse_id in os.listdir(analysis_dir):
                if mouse_id == "grouped":
                    continue
                mouse_path = os.path.join(analysis_dir, mouse_id)

                # Loop across FOVs
                for fov_id in os.listdir(mouse_path):
                    fov_path = os.path.join(mouse_path, fov_id)
                    date_list = os.listdir(fov_path)
                    full_date_list = getattr(self, f"{mouse_id}_DAYS")
                    temp = []
                    i = 0
                    j = 0
                    while i < len(full_date_list) and j < len(date_list):
                        if full_date_list[i] == date_list[j]:
                            temp.append(date_list[j])
                            i += 1
                            j += 1
                        else:
                            temp.append("missing")
                            i += 1
                    date_list = temp
                    # print(date_list)

                    # Load the pickle files in that folder
                    for i, date in enumerate(date_list):
                        if date == "missing":
                            df_list.append(None)
                            continue
                        date_path = os.path.join(fov_path, date)
                        pickle_files = [f for f in os.listdir(date_path) if f.endswith(".pickle")]
                        if len(pickle_files) == 0:
                            print(f"No pickle file in {date_path}")
                            continue
                        pickle_file = pickle_files[0]

                        match = re.search(r"(JL\d+).*_(L\d+)\.pickle", pickle_file) #extract mouse ID & layer
                        layer = match.group(2)

                        try:
                            df_out = utils.get_single_day_df(
                                fname=os.path.join(analysis_dir, mouse_id, fov_id, date, pickle_file),
                                fov_id=fov_id,
                                mouse_id=mouse_id
                            )
                        except ValueError as e:
                            print(f"ERROR in file: {pickle_file}")
                            print(f"  mouse={mouse_id}, fov={fov_id}, date={date}")
                            print(f"  {e}")
                            raise  # re-raise to still get the full traceback

                        # Assign day number and layer ID
                        df_out["day"] = i + 1
                        df_out["layer"] = layer

                        df_list.append(df_out)

            df_list = self.fix_missing_day(df_list)
            # Concatenate across FOVs + days 

            grouped_data = pd.concat(df_list,axis=0, ignore_index=True)

            all_spines = grouped_data["spine_ID"].unique() # get all unique spine IDs

            # Create a MultiIndex with every combination of spine ID x days (1-14), every spine has a row for every day, even if data does not exist
            new_index = pd.MultiIndex.from_product([all_spines,np.arange(1,15)],names=["spine_ID","day"]) 
            # Reindex the DataFrame to the full spine x day, insert NaNs
            grouped_data = grouped_data.set_index(["spine_ID","day"]).reindex(new_index) 
            # Restore spine_ID & day from the index back into regular columns
            grouped_data = grouped_data.reset_index()

            # Reconstruct dendrite_ID by dropping last part of spine_ID
            grouped_data["dendrite_ID"] = grouped_data["spine_ID"].str.split("_").str[:-1].str.join("_")
            # Fill layer label for each dendrite (All spines on same dendrite have the same layer)
            grouped_data["layer"] = grouped_data.groupby(["dendrite_ID"])["layer"].transform("first")

            # Compute dendrite length as maximum pairwise distance between spine positions
            dist = grouped_data.groupby(["dendrite_ID","day"])["position"].apply(lambda x: utils.find_max_pos_diff(np.array(x))).reset_index()
            dist = dist.rename(columns={"position": "dend_length"})

            # Merge dendrite length back into main DataFrame
            grouped_data = grouped_data.merge(dist, on=["day", "dendrite_ID"], how="left")

            # Add spines exist (0/1) to main DataFrame
            grouped_data['exists'] = grouped_data['flags'].apply(
                lambda x: 1 if (x == [] or (isinstance(x, list) and 'New Spine' in x))
                else 0
            )

            if save_grouped:
                with open(savepath, 'wb') as file:
                    pickle.dump(grouped_data, file)

            self.spine_df = grouped_data
            
        else:
            with open(savepath, 'rb') as file:
                self.spine_df = pickle.load(file)

    
    def fix_missing_day(self, df_list):
        """
        Interpolate data for days that were not imaged. 

        For each None entry in df_list (placeholder for missing day), 
        create a DataFrame by averaging "position" and "volume" columns of the previous & following days. 
        Only spines present in both are retained. Day index is incremented accordingly. 

        Parameters
        ----------
        df_list : list of pd.DataFrame or None
            Ordered list of per-day DataFrames. None indicates missing day.

        Returns
        -------
        df_list : list of pd.DataFrame
            Same list returned but with None entries replaced by interpolated DataFrames.
        """
        interpolate_cols = ["position", "volume"]
        for i, df in enumerate(df_list):
            if df is None:
                df_loc = i
                prev_df = df_list[i-1]
                next_df = df_list[i+1]

                # Select only spines that exist day 1 & day 3
                shared_spines = prev_df["spine_ID"].isin(next_df["spine_ID"])

                missing_df = prev_df[shared_spines].copy() # missing_df starts with data from day 1
                next_df_shared = next_df[next_df["spine_ID"].isin(missing_df["spine_ID"])].copy()
                missing_df = missing_df.reset_index(drop=True)
                next_df_shared = next_df_shared.reset_index(drop=True)

                missing_df[interpolate_cols] = (
                    missing_df[interpolate_cols] + next_df_shared[interpolate_cols]
                ) / 2

                missing_df["day"] = missing_df["day"] + 1 # Update day (is day 1 until updated)

                # print(df_loc)
                df_list[df_loc] = missing_df
                    
        return df_list
    
    def spine_turnover(self):
        """
        Calculate new & eliminated spines per µm across consecutive days

        Returns
        -------
        self.turnover : pd.DataFrame
            Columns: dendrite_ID, day, elim_spines, elim_per, new_spines, new_per,
                    turnover_per, neg_elim_per, layer
        """
        self.turnover = calc.spine_turnover(self.spine_df)
        return self

    def spine_lifetime(self, exist_only=True, formed_days=None):
        """
        Calculate raw spine lifetimes & their proportional distributions.

        Parameters
        ----------
        exist_only : bool, optional
            If True, restrict analysis to spines marked as existing. Defaults to True. 
        formed_days : list of int or None, optional
            Subset of days on which to calculate lifetimes. 
            If None, all days are used.

        Returns
        -------
        self.lifetimes : pd.DataFrame
            Per-spine lifetimes. 
        self.lifetime_props : pd.DataFrame
            Proportional distribution of lifetimes. 
        """
        self.lifetimes = calc.spine_lifetime(
            self.spine_df,
            exist_only=exist_only,
            formed_days=formed_days
        )
        self.lifetime_props = calc.get_lifetime_proportions(self.lifetimes)

        return self

    def spine_survival(self, select_days=[6,8,10,12,14], first_days=[1,2,3,4], new_days=[2,3,4]):
        """
        Calculate survival fractions for pre-existing & new spines.

        Parameters
        ----------
        select_days : list, optional
            Days to calculate survival across. Defaults to [6, 8, 10, 12, 14].
        first_days : list, optional
            Days used to define pre-existing spines. Defaults to [1, 2, 3, 4].
        new_days : list, optional
            Days used to define newly formed spines. Defaults to [2, 3, 4].

        Returns
        -------
        self.survival : pd.DataFrame
            Columns: dendrite_ID, layer, day, pre_frac, new_frac.
        """
        self.survival = calc.spine_survival(self.spine_df, select_days, first_days, new_days)
        return self

    def spine_density(self):
        """
        Calculate spine density (spines per µm) for each dendrite over time.

        Returns
        -------
        self.density : pd.DataFrame
            Columns: dendrite_ID, day, density, layer.
        """
        self.density = calc.spine_density(self.spine_df)
        return self

    def spine_volume_change(self):
        """
        Calculate pairwise spine volume change across all days.

        Returns
        -------
        self.vol_change_df : pd.DataFrame
            Spine-level volume changes. Columns: dendrite_ID, day2, vol_change, plasticity, layer.
        self.dend_volchange_df : pd.DataFrame
            Dendrite-level median volume changes. Columns: dendrite_ID, day2, vol_change, layer.
        """
        self.vol_change_df, self.dend_volchange_df = calc.spine_volume_change(self.spine_df)
        return self

    def spine_plasticity_count(self):
        """
        Count the number of spines per plasticity type per dendrite per day.
        Run spine_volume_change() first.

        Returns
        -------
        self.plasticity_df : pd.DataFrame
            Columns: dendrite_ID, day2, plasticity, plasticity_count, n_day1_spines.
        """
        self.plasticity_df = calc.spine_plasticity_count(self.spine_df, self.vol_change_df)
        return self
    
    def volume_lifetime(self):
        """
        Correlate spine volume on day 1 with spine lifetime.
        Requires spine_lifetime() to be run first.

        Returns
        -------
        self.volume_lifetime_df : pd.DataFrame
            Columns: spine_ID, dendrite_ID, volume, layer, lifetime.
        
        Raises
        ------
        AttributeError
            If spine_lifetime() has not been called prior to this method. 
        """
        if not hasattr(self, 'lifetimes'):
            raise AttributeError("Run spine_lifetime() before calling volume_lifetime()")
        
        day1_vol = self.spine_df.query("day == 1 & exists == 1")[["spine_ID", "dendrite_ID", "volume", "layer"]]
        self.volume_lifetime_df = day1_vol.merge(self.lifetimes[["spine_ID", "lifetime"]], on="spine_ID")
        return self