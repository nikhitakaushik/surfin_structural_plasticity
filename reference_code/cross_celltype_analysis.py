import pandas as pd
from Plotting.cross_celltypes_plotter import CrossCelltypePlots

class CrossCelltypeAnalysis:
    def __init__(self, cell_type_dict: dict):
        """
        Combine data across cell types - each stored in SpineData object
        analyses: dict like {"L5_apical": SpineData(...),
                             "L23_apical": SpineData(...),
                             "L23_basal":  SpineData(...)}
        """
        self.analyses = cell_type_dict
        self.plot = CrossCelltypePlots(self)
        self.dendrite_types = list(cell_type_dict.keys())

    # ------------------------------------------------------------------
    # Generic helper: run any per-cell-type method and collect results
    # ------------------------------------------------------------------
    def run_all(self, func_name: str, **kwargs):
        """Call the same method on each dataset and concatenate results."""
        dfs = []
        for fov_type, dataset in self.analyses.items():
            func = getattr(dataset.calc, func_name)
            df = func(**kwargs)
            df["dendrite_type"] = fov_type
            dfs.append(df)
        return pd.concat(dfs, ignore_index=True)
    
    def find_variable_names(self, ds_name, startswith=None, endswith=None, contains=None):
        """
        Return the intersection of variable names (data_vars) matching given patterns
        across all cell-type datasets. If some variables are missing in certain datasets,
        they are printed as a warning (no error raised).

        Parameters
        ----------
        ds_name : str
            Dataset base name (e.g., "activity", "coactivity").
        startswith, endswith, contains : str or list of str, optional
            Match criteria for variable names.

        Returns
        -------
        list of str
            Intersection of matching variable names shared across all datasets.
        """

        # Normalize input patterns
        if isinstance(startswith, str):
            startswith = [startswith]
        if isinstance(endswith, str):
            endswith = [endswith]
        if isinstance(contains, str):
            contains = [contains]

        # --- Collect matching variable sets per dataset ---
        var_sets = {}
        for fov_type, obj in self.analyses.items():
            ds = getattr(obj, f"{ds_name}_data", None)
            if ds is None:
                print(f"⚠️  {fov_type}: No dataset named '{ds_name}_data'. Skipping.")
                continue

            filtered = []
            for var in list(ds.data_vars):
                start_ok = any(var.startswith(p) for p in startswith) if startswith else True
                end_ok = any(var.endswith(s) for s in endswith) if endswith else True
                contains_ok = any(sub in var for sub in contains) if contains else True
                if start_ok and end_ok and contains_ok:
                    filtered.append(var)

            var_sets[fov_type] = set(filtered)

        # --- Compute intersection across all datasets ---
        if not var_sets:
            print(f"No vars found for '{ds_name}'.")
            return []

        shared = set.intersection(*var_sets.values())
        all_unique = set.union(*var_sets.values())

        # --- Print warnings if any variables are missing ---
        if shared != all_unique:
            print(f"⚠️  Variable mismatch across cell types for '{ds_name}':")
            for k, v in var_sets.items():
                missing = all_unique - v
                if missing:
                    print(f"    - {k} missing: {sorted(list(missing))}")

        return sorted(shared)



    # convenience wrappers
    def get_density(self, session="Early", **kwargs):
        return self.run_all("density", session=session, scale=10, **kwargs)

    def get_spine_volumes(self, session="Early", **kwargs):
        return self.run_all("get_spine_properties", session=session, exist_only=True, exclude_shaft=True, **kwargs)

    def get_plasticity(self, days=["Early","Early_followup"], **kwargs):
        return self.run_all("plasticity_by_day",
                            days=days,
                            **kwargs)

