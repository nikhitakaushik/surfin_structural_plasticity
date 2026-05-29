import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from Plotting.plot_parameters import palettes,prettify_axis_labels
from Utilities.stats_collector import StatsCollector

class CrossCelltypePlots:
    """Wraps a CrossCelltypeAnalysis object with plotting functions."""
    def __init__(self, analysis_obj):
        self.group_data = analysis_obj
        self.palettes = palettes
        self.stats = StatsCollector()

    def baseline_summary(self, session="Early", figsize=(10, 3),
                         ylims_list=None, savefig=None):
        """Plot density, volume, and plasticity across cell types."""

        if ylims_list is None:
            ylims_list = [(0, 12), (0, 6), (-0.01, 0.8), (-0.01, 0.8)]

        # Define the metrics to plot
        metrics = [
            {"label": "density",
            "getter": self.group_data.get_density,
            "kwargs": {"session": session},
            "value": "density"},

            {"label": "spine volume",
            "getter": self.group_data.get_spine_volumes,
            "kwargs": {"session": session},
            "value": "spine_volumes"},

            {"label": "sLTP",
            "getter": self.group_data.get_plasticity,
            "kwargs": {"days": [session, session + "_followup"],"group_level":"dendrite"},
            "value": "sLTP"},

            {"label": "sLTD",
            "getter": self.group_data.get_plasticity,
            "kwargs": {"days": [session, session + "_followup"],"group_level":"dendrite"},
            "value": "sLTD"},
        ]

        fig, axes = plt.subplots(1, len(metrics), figsize=figsize, sharex=False)
        plt.subplots_adjust(wspace=0.4)

        for i, m in enumerate(metrics):
            df = m["getter"](**m["kwargs"])
            sns.boxplot(df, x="dendrite_type", y=m["value"], ax=axes[i],
                        hue="dendrite_type", fliersize=0, palette=palettes["dendrite_type"])
            prettify_axis_labels(ax=axes[i],ylim=ylims_list[i])
            axes[i].set_title(m["label"])

            # Collect stats
            self.stats.run_stats(
                name=f"{session}_{m['label']}",
                df=df,
                group_col="dendrite_type",
                value_col=m["value"],
            )

        if savefig:
            plt.savefig(savefig, bbox_inches="tight", dpi=300)
        plt.show()

        return {m["label"]: m["getter"](**m["kwargs"]) for m in metrics}

    def activity_rates(self,session="Early",figsize=(10,4),periods=["session","mvmt"],n_bins=20,plot_type="hist"):
        dtype_list = self.group_data.dendrite_types

        dend_activity_df = self.group_data.run_all("get_activity_data",
                                                var_name="dendrite_activity_rate",
                                                session=session,
                                                period=periods,
                                                level="dendrite",
                                                agg_func="median"
                                                )

        spine_activity_df = self.group_data.run_all("get_activity_data",
                                                var_name="spine_activity_rate",
                                                session=session,
                                                period=periods,
                                                level="spine",              
                                                drop_missing=True             
                                            )

        # Rename for consistency
        spine_activity_df = spine_activity_df.rename(columns={"spine_activity_rate": "activity_rate"})
        dend_activity_df = dend_activity_df.rename(columns={"dendrite_activity_rate": "activity_rate"})

        # Set plotting params based on the data
        xmin = min(dend_activity_df["activity_rate"].min(), spine_activity_df["activity_rate"].min())
        xmax = max(dend_activity_df["activity_rate"].max(), spine_activity_df["activity_rate"].max())

        bins = np.linspace(xmin, xmax, n_bins + 1)
        
        def plot_activity(df, title_suffix, axes,plot_type):

            if plot_type=='hist':
                for i, dt in enumerate(dtype_list):
                    sub = df.query("dendrite_type == @dt")
                    sns.histplot(
                        data=sub,
                        x="activity_rate",
                        hue="period",
                        bins=bins,
                        element="step",
                        fill=False,
                        stat="density",
                        common_norm=False,
                        palette=palettes["move_period"],
                        ax=axes[i]
                    )
                    axes[i].set_title(f"{dt} {title_suffix}", fontsize=12)
                    axes[i].set_xlabel("Activity rate (min$^{-1}$)")
                    axes[i].set_ylabel("Density" if i == 0 else "")
                    axes[i].set_xlim(xmin, xmax)

                    leg = axes[i].get_legend()

                    if leg is not None:
                        leg.set_title("")  # optional: remove "period" or "hue" title
                        leg.get_frame().set_linewidth(0.0)  # remove border line
                        leg.get_frame().set_facecolor('none')  # transparent background
                        for text in leg.get_texts():
                            text.set_fontsize(8)  # smaller font size
            elif plot_type=='line':

                for i, dt in enumerate(dtype_list):
                    sub = df.query("dendrite_type == @dt")
                    sns.lineplot(data=sub,x='period',y='activity_rate',ax=axes[i])
                    axes[i].set_title(f"{dt} {title_suffix}", fontsize=12)


        fig, axes = plt.subplots(2, len(dtype_list),
                                figsize=figsize,
                                sharex=True, sharey=True)

        # Ensure consistent shape when only one dtype
        if len(dtype_list) == 1:
            axes = np.array(axes).reshape(2, 1)

        plot_activity(spine_activity_df, "Spine", axes[0],plot_type)
        plot_activity(dend_activity_df, "Dendrite", axes[1],plot_type)

        plt.tight_layout()
        plt.show()

        return spine_activity_df,dend_activity_df

    def MSI_index(self,session="Early",periods=["mvmt","nonmvmt"],figsize=(2,3),ylims=(-1,1),xlims=(-1,1),measure="index",plot_type="dots"):

        MSI_df = self.group_data.run_all("movement_selectivity_index",
                                         session=session,
                                         periods=periods)
        
        fig,ax = plt.subplots(1,1,figsize=figsize)

        if measure == 'index':
            data_col = 'MSI'
        elif measure == 'ratio':
            data_col = 'x_over_y'

        if plot_type == 'dots':
            sns.stripplot(MSI_df,
                        x="dendrite_type",
                        y=data_col,
                        alpha=0.4,
                        palette=palettes["dendrite_type"],
                        hue="dendrite_type",
                        s=8,
                        legend=None
                        )

            sns.pointplot(MSI_df,
                        x="dendrite_type",
                        y=data_col,
                        errorbar="se",
                        ax=ax,
                        dodge=0.5,
                        linestyles="none",
                        hue="dendrite_type",
                        palette=palettes["dendrite_type"],
                        markersize=10,
                        legend=None
                        )

            ax.set_ylim(ylims)
        elif plot_type == 'hist':

            sns.histplot(data=MSI_df,
                        x=data_col,
                        hue="layer",
                        stat="proportion",
                        common_norm=False,
                        palette=palettes["layer"]
                        )
            
            # Calculate the median
            median_values = MSI_df.groupby("layer")[data_col].median()
            # median_values
            # Add a vertical line for the median
            for layer, value in median_values.items():
                plt.vlines(x=value,
                        ymin=0, ymax=0.5,linestyles='--',
                        color=palettes["layer"][layer], lw=2)

            if ylims[0] > 0:
                ax.set_ylim(ylims)

            ax.set_xlim(xlims)

        # Collect stats
        self.stats.run_stats(
            name=f"{session}_{periods[0]}_{periods[1]}_selectivity_index",
            df=MSI_df,
            group_col="layer",
            value_col=data_col,
        )

        self.stats.run_mixedlm_test(
            name=f"{session}_{periods[0]}_{periods[1]}_selectivity_index",
            df=MSI_df,
            value_col=data_col,
            fixed_effect="layer",
            random_effect="FOV",
            )

        return MSI_df