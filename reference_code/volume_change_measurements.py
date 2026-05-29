import pandas as pd
import numpy as np

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