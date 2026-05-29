import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2
import numpy as np
import statsmodels.formula.api as smf
import scikit_posthocs as sp


def run_mixedlm_test(
                    df,
                    value_col,
                    factor,
                    random_effect,
                    other_factors=None,
                    include_interaction=False,
                    alpha=0.05,
                    name=None
    ):

        df = df.copy().dropna(subset=[value_col, factor, random_effect])
        df = df.reset_index(drop=True)

        # --- categorical ---
        df[factor] = df[factor].astype("category")
        df[random_effect] = df[random_effect].astype("category")

        if other_factors:
            if isinstance(other_factors,str):
                other_factors = [other_factors]
            for f in other_factors:
                df[f] = df[f].astype("category")

        # =====================================================
        # BUILD FORMULAS
        # =====================================================
        if other_factors:
            fixed_terms = [factor] + other_factors

            if include_interaction:
                full_formula = f"{value_col} ~ " + " * ".join(fixed_terms)
                main_formula = f"{value_col} ~ " + " + ".join(fixed_terms)
            else:
                full_formula = f"{value_col} ~ " + " + ".join(fixed_terms)
                main_formula = full_formula

        else:
            fixed_terms = [factor]
            full_formula = f"{value_col} ~ {factor}"
            main_formula = full_formula

        # --- fit full model (ML required for LRT) ---
        full_model = smf.mixedlm(
            full_formula,
            data=df,
            groups=df[random_effect]
        ).fit(reml=False)

        rows = []

        # =====================================================
        # LRT HELPER
        # =====================================================
        def lrt(full, reduced):
            lr = 2 * (full.llf - reduced.llf)
            df_diff = full.df_modelwc - reduced.df_modelwc
            return chi2.sf(lr, df_diff)

        # =====================================================
        # 1. TEST MAIN EFFECT (factor)
        # =====================================================
        if other_factors:
            reduced_terms = [t for t in fixed_terms if t != factor]

            if reduced_terms:
                reduced_formula = f"{value_col} ~ " + " + ".join(reduced_terms)
            else:
                reduced_formula = f"{value_col} ~ 1"
        else:
            reduced_formula = f"{value_col} ~ 1"

        reduced_model = smf.mixedlm(
            reduced_formula,
            data=df,
            groups=df[random_effect]
        ).fit(reml=False)

        p_factor = lrt(full_model, reduced_model)

        rows.append({
            "comparison": name,
            "test": "lrt",
            "term": factor,
            "p": p_factor,
            "significant": p_factor < alpha
        })

        # =====================================================
        # 2. TEST INTERACTION (if requested)
        # =====================================================
        p_interaction = None

        if include_interaction and other_factors:
            no_inter_model = smf.mixedlm(
                main_formula,
                data=df,
                groups=df[random_effect]
            ).fit(reml=False)

            p_interaction = lrt(full_model, no_inter_model)

            rows.append({
                "comparison": name,
                "test": "lrt",
                "term": "interaction",
                "p": p_interaction,
                "significant": p_interaction < alpha
            })

        # =====================================================
        # 3. WALD TEST (optional, nice to have)
        # =====================================================
        wald = full_model.wald_test_terms()

        for term, stat in wald.table.iterrows():
            rows.append({
                "comparison": name,
                "test": "wald",
                "term": term,
                "p": stat["pvalue"],
                "significant": stat["pvalue"] < alpha
            })

        results = pd.DataFrame(rows)

        # =====================================================
        # 4. POSTHOC (DUNN)
        # =====================================================
        posthoc = None

        if p_factor < alpha:
            dunn = sp.posthoc_dunn(
                df,
                val_col=value_col,
                group_col=factor,
                p_adjust=None
            )

            posthoc = (
                dunn.reset_index()
                .melt(id_vars="index", var_name="group_2", value_name="p_raw")
                .rename(columns={"index": "group_1"})
            )

            posthoc = posthoc.query("group_1 != group_2").copy()
            # Remove symmetric duplicates
            pairs = np.sort(posthoc[["group_1", "group_2"]].values, axis=1)
            posthoc = posthoc.loc[
                ~pd.DataFrame(pairs).duplicated().values
            ].copy()
            
            posthoc["p_fdr"] = multipletests(
                posthoc["p_raw"],
                method="fdr_tsbh"
            )[1]

            posthoc["comparison"] = name
            posthoc["test"] = "dunn"
            posthoc["term"] = factor

        # =====================================================
        # FINAL TABLE
        # =====================================================
        if posthoc is not None:
            final = pd.concat([results, posthoc], ignore_index=True)
        else:
            final = results

        return final


def lrt(full_model, reduced_model):
    lr = 2 * (full_model.llf - reduced_model.llf)
    df_diff = full_model.df_modelwc - reduced_model.df_modelwc
    return chi2.sf(lr, df_diff)