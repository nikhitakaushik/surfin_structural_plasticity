import pandas as pd
import numpy as np

def spine_turnover(df):
        """
        Calculate new & elimiated spines per µm across consecutive days. 

        For each pair of adjacent days (1-2, 2-3,...13-14), counts the number of spines that were : 
        Present on day A & absent on day B (eliminations) & Absent of day A & present on day B (additions)
        Also sums addition/elimination to calculate spine turnover. 

        Parameters
        ----------
        df : pd.DataFrame
            Spine-level DataFrame with columns: dendrite_ID, day, exists, dend_length, layer.

        Returns
        -------
        turnover_df : pd.DataFrame
            Columns: dendrite_ID, day, elim_spines, elim_per, new_spines, new_per, turnover_per, neg_elim_per, layer.
        """
        rows = []
        for day_b in range(2, 15):
            day_a = day_b - 1

            day1 = df[df["day"] == day_a].groupby("dendrite_ID")["exists"].apply(np.array)
            day2 = df[df["day"] == day_b].groupby("dendrite_ID")["exists"].apply(np.array)
            median_length = df[df["day"] == 1].groupby("dendrite_ID")["dend_length"].apply(np.nanmedian)

            for dendrite in day1.index.intersection(day2.index):
                exist_1    = day1[dendrite] == 1
                missing_1  = day1[dendrite] == 0
                eliminated_2 = day2[dendrite] == 0
                new_2        = day2[dendrite] == 1

                dendrite_length = median_length.loc[dendrite]
                elim_spines = np.sum(exist_1 & eliminated_2)
                new_spines  = np.sum(missing_1 & new_2)

                rows.append({
                    "dendrite_ID": dendrite,
                    "day":         day_b,
                    "elim_spines": elim_spines,
                    "elim_per":    elim_spines / dendrite_length,
                    "new_spines":  new_spines,
                    "new_per":     new_spines / dendrite_length,
                })

        turnover_df = pd.DataFrame(rows)

        layer_map = df.groupby("dendrite_ID")["layer"].first()
        turnover_df["layer"] = turnover_df["dendrite_ID"].map(layer_map)

        turnover_df["turnover_per"] = turnover_df["new_per"] + turnover_df["elim_per"]

        turnover_df["neg_elim_per"] = -turnover_df["elim_per"]

        return turnover_df

def calc_lifetime(x):
    """
    Calculate lifetime of a single spine. 

    Lifetime defined as the number of days from first appearance to first disappearance 
    (or end of imaged days if never lost).
    
    Parameters
    ----------
    x : array-like of int
        Binary existence array for one spine across days (0=absent, 1=present).

    Returns
    -------
    lifetime : int
        Number of days the spine was alive, including first & last day present.  
    """
    x = np.asarray(x, dtype=int)

    birth = np.argmax(x == 1)
    death = birth + np.argmax(x[birth:] == 0) if np.any(x[birth:] == 0) else len(x) - 1
    lifetime = death - birth + 1

    return lifetime

def spine_lifetime(df, exist_only=False, formed_days=None):
    """
    Calculate lifetime of each spine.
    
    Parameters
    ----------
    df : pd.DataFrame
        Spine-level DataFrame with columns spine_ID, dendrite_ID, day, exists, layer
    exist_only : bool, optinal
        If True, only include spines present on day 1.
    formed_days :list of int or None, optional
        If provided, only include spines whose first day of existence
        is in this collection (e.g. [2, 3, 4]). If None, all spines included.

    Returns
    -------
    lifetimes : pd.DataFrame
        Columns: dendrite_ID, spine_ID, lifetime, layer.
    """
    filtered_df = df.copy()

    if exist_only:
        day1_spines = filtered_df[(filtered_df['day'] == 1) & (filtered_df['exists'] == 1)]
        day1_spineids = day1_spines['spine_ID'].unique()
        filtered_df = filtered_df[filtered_df['spine_ID'].isin(day1_spineids)]

    if formed_days is not None:
        first_appearance = (
            filtered_df[filtered_df['exists'] == 1]
            .groupby(['dendrite_ID', 'spine_ID'], sort=False)['day']
            .min()
            .reset_index(name='first_day')
        )
        
        valid_spines = first_appearance[first_appearance['first_day'].isin(formed_days)]

        filtered_df = filtered_df.merge(
            valid_spines[['dendrite_ID', 'spine_ID']],
            on=['dendrite_ID', 'spine_ID'],
            how='inner'
        )

    x = filtered_df.groupby(["dendrite_ID", "spine_ID"], sort=False)["exists"]
    lifetimes = x.apply(calc_lifetime).reset_index()
    lifetimes.columns = ["dendrite_ID", "spine_ID", "lifetime"]

    layer_map = df.groupby("dendrite_ID")["layer"].first()
    lifetimes["layer"] = lifetimes["dendrite_ID"].map(layer_map)

    return lifetimes

def get_lifetime_proportions(lifetimes_df):
    """
    Calculate the percentage of spines falling into each lifetime bin, by layer. 

    Parameters
    ----------
    lifetimes_df : pd.DataFrame
        Output of spine_lifetime(). Columns: dendrite_ID, spine_ID, lifetime, layer. 

    Returns
    -------
    df_prop : pd.DataFrame
        Columns: layer, lifetime, count, total, proportion. 
        Proportion expressed as percentage (1-100%)
    """
    # Count occurrences per layer and lifetime
    counts = lifetimes_df.groupby(['layer', 'lifetime']).size().reset_index(name='count')
    
    # Calculate totals per layer to divide by
    totals = lifetimes_df.groupby('layer').size().reset_index(name='total')
    
    # Merge and calculate proportion
    df_prop = counts.merge(totals, on='layer')
    df_prop['proportion'] = (df_prop['count'] / df_prop['total']) * 100
    
    return df_prop

def spine_survival(df, select_days, first_days, new_days):
    """
    Calculate survival fractions for pre-exisiting & newly formed spines. 

    Pre-exisiting spines defined as present on day 1. 
    New spines defined as those absent on day 1 but first appearing within new_days. 
    Survival fraction is the proportion of each group still present of each of select_days. 

    Parameters
    ----------
    df : pd.DataFrame
        Spine-level DataFrame with columns spine_ID, dendrite_ID, day, exists, layer.
    select_days : list of int   
        Days at which survival fraction is evaluated
    first_days : list of int
        Days used to define spine population (pre & new)
    new_days : list of int
        Days within first_days on which spine must first appear to be "new"

    Returns
    -------
    survival_df : pd.DataFrame
        Columns: dendrite_ID, layer, day, pre_frac, new_frac.
    """
    all_days = first_days + select_days
    wide = (df.pivot_table(
                index=["dendrite_ID", "spine_ID", "layer"],
                columns="day",
                values="exists",
                aggfunc="max",
                fill_value=0)
            .reindex(columns=all_days, fill_value=0)
            .astype(int))

    rows = []
    for (dendrite_id, layer), group in wide.groupby(level=["dendrite_ID", "layer"]):
        select_spines = group[first_days].any(axis=1)
        group_select  = group.loc[select_spines]

        prespines_mask = group_select[1] == 1
        first_appear   = group_select[first_days].idxmax(axis=1)
        newspines_mask = (group_select[1] == 0) & first_appear.isin(new_days)

        pre_frac = group_select.loc[prespines_mask, select_days].mean(axis=0)
        new_frac = group_select.loc[newspines_mask, select_days].mean(axis=0)

        for day in select_days:
            rows.append({
                "dendrite_ID": dendrite_id,
                "layer": layer,
                "day": day,
                "pre_frac": pre_frac[day],
                "new_frac": new_frac[day]
            }

            )

    return pd.DataFrame(rows)

def spine_density(df):
    """
    Calculate spine density (spines per µm) for each dendrite over days. 

    Density calculated as number of existing spines divided by mean dendrite length. 
    Also returns a normalized density column scaled to each dendrite's value on day 1. 

    Parameters
    ----------
    df : pd.DataFrame
        Spine-level DataFrame with columns: spine_ID, dendrite_ID, day, exists, dend_length, layer.

    Returns
    -------
    density : pd.DataFrame
        Columns: dendrite_ID, day, density, layer, density_norm
    """
    exist_df = df[df['exists'] == 1] 
    spine_count = exist_df.groupby(["dendrite_ID", "day"])["spine_ID"].count()
    lengths = exist_df.groupby("dendrite_ID")["dend_length"].mean()

    density = (spine_count/lengths).reset_index()
    density.columns = ["dendrite_ID", "day", "density"]

    layer_map = df.groupby("dendrite_ID")["layer"].first()
    density["layer"] = density["dendrite_ID"].map(layer_map)

    density["density_norm"] = density.groupby("dendrite_ID")["density"].transform(
        lambda x: x / x.iloc[0]
    )

    return density

def pairwise_plasticity(df, days=range(1, 15), vol_col="spine_volume"):
    out = []

    for d1, d2 in zip(days[:-1], days[1:]):

        # --- filter valid spines ---
        sub = df[
            df["day"].isin([d1, d2]) &
            (df["exists"]) 
        ]

        # --- count spines on day 1 (per dendrite) ---
        day1_counts = (
            sub[sub["day"] == d1]
            .groupby("dendrite_ID")["spine_ID"]
            .nunique()
            .rename("n_day1_spines")
        )

        # --- pivot to get volumes ---
        pivot = (
            sub.pivot_table(
                index=["dendrite_ID", "spine_ID"],
                columns="day",
                values=vol_col
            )
            .dropna(subset=[d1, d2])
            .reset_index()
        )

        # --- compute change ---
        pivot["vol_change"] = pivot[d2] / pivot[d1]

        # --- classify plasticity ---
        conditions = [
            pivot["vol_change"] < 0.75,
            pivot["vol_change"] > 1.5
        ]
        choices = ["sLTD", "sLTP"]

        pivot["plasticity"] = np.select(conditions, choices, default="stable")

        # --- attach day1 counts ---
        pivot = pivot.merge(day1_counts, on="dendrite_ID", how="left")

        # --- metadata ---
        pivot["day1"] = d1
        pivot["day2"] = d2

        layer_map = df.groupby("dendrite_ID")["layer"].first()
        pivot["layer"] = pivot["dendrite_ID"].map(layer_map)

        out.append(pivot[[
            "dendrite_ID",
            "spine_ID",
            "day1",
            "day2",
            "vol_change",
            "plasticity",
            "n_day1_spines",
            "layer"
        ]])

    return pd.concat(out, ignore_index=True)

def spine_volume_change(df, days=range(1, 15), vol_col="volume"):
    """
    Calculate pairwise spine volume change across all days. 

    Calls pairwise_plasticity() to compute spine-level volume changes. 
    Then aggregates to dendrite-level median values.

    Parameters
    ----------
    df : pd.DataFrame 
        Spine-level DataFrame with columns: spine_ID, dendrite_ID, day, volume, layer.
    days : range or list of int, optional
        Days over which to compute pairwise changes. Defaults to range (1,15).
    vol_col : str, optional
        Name of the volume column. Defaults to "volume". 

    Returns
    -------
    vol_change_df : pd.DataFrame
        Spine-level volume changes. 
        Columns: dendrite_ID, day2, vol_change, plasticity, layer.
    dend_vol_change_df : pd.DataFrame
        Dendrite-level median volume changes
        Columns: dendrite_ID, day2, vol_change, layer.
    """
    vol_change_df = pairwise_plasticity(df, days=days, vol_col=vol_col)
    
    # aggregate to dendrite level
    layer_map = df.groupby("dendrite_ID")["layer"].first()
    vol_change_df["layer"] = vol_change_df["dendrite_ID"].map(layer_map)

    agg_function = {'layer': 'first', 'vol_change': np.median}
    dend_volchange_df = vol_change_df.groupby(["dendrite_ID", "day2"]).agg(agg_function).reset_index()

    return vol_change_df, dend_volchange_df


def spine_plasticity_count(df, vol_change_df):
    """
    Count the number of spines per plasticity type per dendrite per day. 
    Run spine_volume_change() first. 

    Parameters
    ----------
    df : pd.DataFrame
        Spine-level DataFrame used to look up layer for each dendrite. 
    vol_change_df : pd.DataFrame
        Output of spine_volume_change(). Columns: dendrite_ID, day2, plasticity, n_day1_spines, layer
    
    Returns
    -------
    plasticity_df : pd.DataFrame
        Columns: dendrite_ID, day2, plasticity, plasticity_count, plasticity_proportion, n_day1_spines, layer
    """
    agg_function = {'n_day1_spines': 'first', 'plasticity': 'count'}
    plasticity_df = (vol_change_df.groupby(["dendrite_ID", "day2", "plasticity"])
                                  .agg(agg_function)
                                  .rename(columns={"plasticity": "plasticity_count"})
                                  .reset_index())
    plasticity_df["plasticity_proportion"] = (
        plasticity_df["plasticity_count"] / plasticity_df["n_day1_spines"]
    )
    layer_map = df.groupby("dendrite_ID")["layer"].first()
    plasticity_df["layer"] = plasticity_df["dendrite_ID"].map(layer_map)

    return plasticity_df