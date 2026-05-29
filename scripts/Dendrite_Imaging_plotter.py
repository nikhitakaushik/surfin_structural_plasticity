import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scipy.stats as stats
 
 
# Plot style
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["figure.titlesize"] = 10

plotcolors = {
    'plasticity' : {
        'sLTP': "#351acc",
        'sLTD': "#cc3333",
        'stable': "#666666"
    },

    'add' : {
        'L23': "#62AD62",
        'L5' : "#b56bbe"
    },

    'elim' : {
        'L23': "#1E381E",
        'L5' : "#310936"
    },

    'layer' : {
        'L23': "#1D551D",
        'L5' : "#6e1f79"
    }
}

ALL_DAYS = list(range(1, 15))
TURNOVER_DAYS = list(range(2, 15))
 
 
class DendriteImagingPlots:
    """
    Plotting methods for a Dendrite_Imaging object
 
    Usage
    -----
    dend = Dendrite_Imaging(analysis_dir)
    dend.load_grouped_data()
    dend.plot = DendriteImagingPlots(dend)
 
    dend.plot.density()
    dend.plot.add_elim_line()
    dend.plot.add_elim_bar()
    dend.plot.turnover()
    dend.plot.survival()
    dend.plot.volume_change()
    dend.plot.plasticity_events()
    dend.plot.lifetime()
    """
 
    def __init__(self, imaging_obj):
        self.data = imaging_obj
 
# Plots
 
    def density(self, figsize=(4, 3), ylim=(0.1, 0.5), savefig=None):
        """
        Spine density per day, grouped by layer.
 
        Parameters
        ----------
        figsize : tuple
        ylim : tuple
        savefig : str or None
        """
        fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(
            data=self.data.density,
            x="day",
            y="density",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=ax
        )
        ax.set_title("Spine Density per Day")
        ax.set_xticks(ALL_DAYS)
        ax.set_xticklabels(ALL_DAYS)
        ax.set_ylim(ylim)
        ax.set(xlabel="Day", ylabel="Spine density (spines / µm)")
        ax.legend(title="Layer")
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()

    def density_norm(self, figsize=(4, 3), ylim=(0, 2), savefig=None):
        """
        Spine density per day, grouped by layer.
 
        Parameters
        ----------
        figsize : tuple
        ylim : tuple
        savefig : str or None
        """
        fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(
            data=self.data.density,
            x="day",
            y="density_norm",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=ax
        )
        ax.set_title("Spine Density per Day")
        ax.set_xticks(ALL_DAYS)
        ax.set_xticklabels(ALL_DAYS)
        ax.set_ylim(ylim)
        ax.set(xlabel="Day", ylabel="Normalized spine density")
        ax.legend(title="Layer")
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def add_elim_line(self, figsize=(10, 3), ylim=(-0.01, 0.1), savefig=None):
        """
        Line plot of new & eliminated spines per µm over days, grouped by layer.
 
        Parameters
        ----------
        figsize : tuple
        ylim : tuple
        savefig : str or None
        """
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True, sharex=True)
 
        sns.lineplot(
            data=self.data.turnover,
            x="day",
            y="new_per",
            hue="layer",
            palette=plotcolors["add"],
            errorbar="se",
            ax=axes[0]
        )
        axes[0].set_title("Number of Spines Added per Day")
        axes[0].set_xticks(TURNOVER_DAYS)
        axes[0].set_xticklabels(TURNOVER_DAYS)
        axes[0].set_ylim(ylim)
        axes[0].set(xlabel="Day", ylabel="Spines Added / µm")
        axes[0].legend(title="Layer")
 
        sns.lineplot(
            data=self.data.turnover,
            x="day",
            y="elim_per",
            hue="layer",
            palette=plotcolors["elim"],
            errorbar="se",
            ax=axes[1]
        )
        axes[1].set_title("Number of Spines Eliminated per Day")
        axes[1].set_xticks(TURNOVER_DAYS)
        axes[1].set_xticklabels(TURNOVER_DAYS)
        axes[1].set_ylim(ylim)
        axes[1].set(xlabel="Day", ylabel="Spines Eliminated / µm")
        axes[1].legend(title="Layer")
 
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def add_elim_bar(self, figsize=(13, 4), ylim=(-0.1, 0.1), savefig=None):
        """
        Mirrored bar chart of added & eliminated spines per µm, one panel per layer
    
        Parameters
        ----------
        figsize : tuple
        ylim : tuple
        savefig : str or None
        """
        layers = self.data.turnover["layer"].unique()
        fig, axes = plt.subplots(1, len(layers), figsize=figsize, sharey=True, sharex=True)
    
        for ax, layer in zip(axes, layers):
            layer_df = self.data.turnover.query("layer == @layer")
            
            # Plot eliminated spines (negative) using 'elim' colors
            sns.barplot(data=layer_df, x="day", y="neg_elim_per", ax=ax, 
                    color=plotcolors['elim'][layer])
            
            # Plot added spines (positive) using 'add' colors
            sns.barplot(data=layer_df, x="day", y="new_per", ax=ax, 
                    color=plotcolors['add'][layer])
            
            ax.axhline(0, linewidth=1, color="black")
            ax.set_title(f"Added and Eliminated Spines - {layer}")
            ax.set_xticks([x - 2 for x in TURNOVER_DAYS])
            ax.set_xticklabels(TURNOVER_DAYS)
            ax.set_ylim(ylim)
            ax.set(xlabel="Day", ylabel="Spines Added & Eliminated / µm")
    
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def turnover(self, figsize=(7, 3), ylim=(-0.03, 0.2), savefig=None):
        """
        Spine turnover rate per day, grouped by layer
 
        Parameters
        ----------
        figsize : tuple
        ylim : tuple
        savefig : str or None
        """
        fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(
            data=self.data.turnover,
            x="day",
            y="turnover_per",
            hue="layer",
            palette=plotcolors["layer"],
            ax=ax
        )
        ax.set_title("Spine Turnover Rate per Day")
        ax.set_xticks(TURNOVER_DAYS)
        ax.set_xticklabels(TURNOVER_DAYS)
        ax.set_ylim(ylim)
        ax.set(xlabel="Day", ylabel="Turnover rate / µm")
        ax.legend(title="Layer")
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def survival(self, figsize=(10, 3), savefig=None):
        """
        Survival curves for pre-existing and newly formed spines, grouped by layer
    
        Parameters
        ----------
        figsize : tuple
        savefig : str or None
        """
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True, sharex=True)
    
        sns.lineplot(
            data=self.data.survival,
            x="day",
            y="pre_frac",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=axes[0],
            linestyle="-"  # solid line for pre-existing
        )
        axes[0].set_title("Pre-existing Fraction Survived")
        axes[0].set(xlabel="Day", ylabel="Fraction of Spines Present")
        axes[0].legend(title="Layer")
    
        sns.lineplot(
            data=self.data.survival,
            x="day",
            y="new_frac",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=axes[1],
            linestyle="--"  # dashed line for newly-formed
        )
        axes[1].set_title("Newly-formed Fraction Survived")
        axes[1].set(xlabel="Day", ylabel="Fraction of Spines Present")
        axes[1].legend(title="Layer")
    
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def volume_change(self, figsize=(10, 3), savefig=None):
        """
        Spine volume change over days at spine and dendrite level, grouped by layer
 
        Parameters
        ----------
        figsize : tuple
        savefig : str or None
        """
        all_days_list = sorted(self.data.vol_change_df["day2"].unique())
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True, sharex=True)
 
        sns.lineplot(
            data=self.data.vol_change_df,
            x="day2",
            y="vol_change",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=axes[0]
        )
        axes[0].set_title("Volume Change (Spine Level)")
        axes[0].set_xticks(all_days_list)
        axes[0].set_xticklabels(all_days_list)
        axes[0].set(xlabel="Day", ylabel="Volume Change")
        axes[0].legend(title="Layer")
 
        sns.lineplot(
            data=self.data.dend_volchange_df,
            x="day2",
            y="vol_change",
            hue="layer",
            palette=plotcolors["layer"],
            errorbar="se",
            ax=axes[1]
        )
        axes[1].set_title("Volume Change (Dendrite Level)")
        axes[1].set_xticks(all_days_list)
        axes[1].set_xticklabels(all_days_list)
        axes[1].set(xlabel="Day", ylabel="Volume Change")
        axes[1].legend(title="Layer")
 
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def plasticity_events(self, figsize=(10, 3), savefig=None):
        """
        Spine plasticity type (sLTD/sLTP/stable), grouped by layer
 
        Parameters
        ----------
        figsize : tuple
        savefig : str or None
        """
        all_days_list = sorted(self.data.plasticity_df["day2"].unique())
        layers = sorted(self.data.plasticity_df["layer"].unique())
        
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True, sharex=True)

        for i, layer in enumerate(layers):
            layer_data = self.data.plasticity_df[self.data.plasticity_df["layer"] == layer]
            
            sns.lineplot(
                data=layer_data,
                x="day2",
                y="plasticity_proportion", # Changed from count
                hue="plasticity",
                palette=plotcolors["plasticity"],
                errorbar="se",
                ax=axes[i]
            )
            
            axes[i].set_title(f"Plasticity Proportion: {layer}")
            axes[i].set_xticks(all_days_list)
            axes[i].set_xticklabels(all_days_list)
            
            # Adjusting labels for clarity
            axes[i].set(xlabel="Day", ylabel="Proportion of Spines")
            axes[i].legend(title="Plasticity Type", loc='upper right')

        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()
 
    def lifetime(self, figsize=(12, 3), savefig=None):
        """
        Distribution of spine lifetimes (# days existed on dendrite), grouped by layer
 
        Parameters
        ----------
        figsize : tuple
        savefig : str or None
        """
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
        all_lifetimes = range(1, 15)
        for i, layer_name in enumerate(['L23', 'L5']):
            subset = (self.data.lifetime_props
                    .query(f"layer == '{layer_name}'")
                    .set_index("lifetime")
                    .reindex(all_lifetimes, fill_value=0)
                    .reset_index())
            sns.barplot(
                data=subset,
                x="lifetime",
                y="proportion",
                ax=axes[i],
                color=plotcolors['layer'][layer_name]
            )
            axes[i].set_title(f"Spine Lifetime Distribution - {layer_name}")
            axes[i].set_xlabel("Spine Lifetime (days)")
            axes[i].set_ylabel("Percentage of Spines (%)" if i == 0 else "")
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()

    def volume_lifetime_correlation(self, figsize=(10, 4), savefig=None):
        """
        Scatter plot of spine volume on day 1 vs lifetime, one panel per layer.

        Parameters
        ----------
        figsize : tuple
        savefig : str or None
        """
        layers = self.data.volume_lifetime_df["layer"].unique()
        fig, axes = plt.subplots(1, len(layers), figsize=figsize, sharey=True, sharex=True)

        titles = {"L23": "L2/3", "L5": "L5"}

        for ax, layer in zip(axes, layers):
            layer_df = self.data.volume_lifetime_df.query("layer == @layer")
            r, p = stats.pearsonr(layer_df["volume"], layer_df["lifetime"])

            sns.regplot(data=layer_df, x="volume", y="lifetime", ax=ax,
                        scatter=False, color=plotcolors["layer"][layer])
            sns.scatterplot(data=layer_df, x="volume", y="lifetime", ax=ax,
                            color=plotcolors["layer"][layer])
            ax.set_title(f"{titles[layer]}\nr = {r:.2f}, p = {p:.3f}")
            ax.set_xlabel("Spine Volume (day 1)")
            ax.set_ylabel("Spine Lifetime (days)" if ax == axes[0] else "")

        fig.suptitle("Spine Volume vs Lifetime")
        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=600, bbox_inches="tight")
        plt.show()