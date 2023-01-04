# pylint: disable=line-too-long

"""Plot class for ReCo."""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .attribute_descriptors import LoggerDesc


class PlotPanel:
    """
    Class representing a matplotlib plot panel.

    Parameters
    ----------
    logger : logging.Logger
        A logging.Logger object.
    sample_name : str
        The name of the sample.
    sample_type : str
        The type of the sample.
    read_count_df : pd.DataFrame
        A pandas dataframe containing the read counts.
    expected_reads : int
        The number of expected reads.
    observed_reads : int
        The number of observed reads.
    trimmed_reads : int
        The number of trimmed reads.
    aligned_reads : int
        The number of aligned reads.
    template_marker : str
        A string that appears as a substring in the name of the template sequence.
    output_dir : str
        The path to a directory to which the plot is saved as a pdf ana png file.

    Methods
    -------
    generate
        Generate and save plots.
    """
    logger = LoggerDesc()

    def __init__(self, logger=None, sample_name=None, sample_type=None, read_count_df=None,
                 expected_reads=0, observed_reads=0, trimmed_reads=0, aligned_reads=0,
                 template_marker="template", output_dir=None):
        self.logger = logger
        self.sample_name = sample_name
        self.sample_type = sample_type
        self.read_count_df = read_count_df

        self.expected_reads = expected_reads
        self.observed_reads = observed_reads
        self.trimmed_reads = trimmed_reads
        self.aligned_reads = aligned_reads

        self.template_marker = template_marker

        self.output_dir = output_dir
        self.pdf_file = os.path.join(self.output_dir, f"{self.sample_name}_plot_panel.pdf")
        self.png_file = os.path.join(self.output_dir, f"{self.sample_name}_plot_panel.png")

    def generate(self):
        """
        Generate plot panel and save to files.
        Returns
        -------

        """
        fig = plt.figure(figsize=(20, 20))
        fig.suptitle(f"{self.sample_name}\n"
                     f"Expected reads: {self.expected_reads:n}\n"
                     f"Observed reads: {self.observed_reads:n} ({round(self.observed_reads / (self.expected_reads / 100), 2)}% of expected)",
                     y=1.0)
        grid_spec = fig.add_gridspec(8, 16)

        self._rates_bar_plot(fig, grid_spec)
        self._dist_box_plot(fig, grid_spec)
        self._dist_plot(fig, grid_spec)
        self._cumul_plot(fig, grid_spec)
        self._histogram(fig, grid_spec)

        plt.tight_layout()

        plt.savefig(self.pdf_file, bbox_inches="tight")
        plt.savefig(self.png_file, bbox_inches="tight")
        plt.close()

    def _rates_bar_plot(self, fig, grid_spec):
        """
        Generate a bar plot representing the trimming and alignment rates.

        Parameters
        ----------
        fig : plt.Figure
            A matplotlib figure object
        grid_spec : plt.gridspec.GridSpect
            A matplotlib GridSpec object.

        Returns
        -------
        ax1 : plt.axes.Axes
            A matplitlib Axes object.
        """
        ax1 = fig.add_subplot(grid_spec[0:4, 0:4])
        ax1.set_title("Trimming/alignment rates")

        trimmed_rate = 100 - (self.trimmed_reads / self.observed_reads * 100)
        aligned_rate = 100 - (self.aligned_reads / self.observed_reads * 100)

        aligned_bar = [100 - aligned_rate]
        invalid_bar = [aligned_rate - trimmed_rate]
        valid_bar = [trimmed_rate]

        x_pos = [0]
        bar_width = 0.75
        ax1.set_xlim([-1, 1])
        names = [self.sample_name]

        ax1.bar(
            x_pos,
            aligned_bar,
            color="green",
            edgecolor="white",
            width=bar_width,
            label=f"Aligned: {round(100 - aligned_rate, 2)}%",
        )
        ax1.bar(
            x_pos,
            invalid_bar,
            bottom=aligned_bar,
            color="salmon",
            edgecolor="white",
            width=bar_width,
            label=f"Not aligned: {round(aligned_rate - trimmed_rate, 2)}%",
        )
        ax1.bar(
            x_pos,
            valid_bar,
            bottom=[i + j for i, j in zip(aligned_bar, invalid_bar)],
            color="#a3acff",
            edgecolor="white",
            width=bar_width,
            label=f"Not trimmed: {round(trimmed_rate, 2)}%",
        )

        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(names[0:])
        ax1.set_ylabel("%")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1_legend = ax1.legend(loc="best", frameon=False)
        for legend_handle in ax1_legend.legendHandles:
            legend_handle.set_alpha(1)

        return ax1

    def _dist_box_plot(self, fig, grid_spec):
        """
        Generate a boxplot representing the distribution of read counts, shows expected and observed median count.

        Parameters
        ----------
        fig : plt.Figure
            A matplotlib figure object
        grid_spec : plt.gridspec.GridSpect
            A matplotlib GridSpec object.

        Returns
        -------
        ax2 : plt.axes.Axes
            A matplitlib Axes object.
        """
        ax2 = fig.add_subplot(grid_spec[0:4, 4:8])
        my_median = self.read_count_df[self.sample_name].median()
        expected_mean = self.expected_reads / self.read_count_df.shape[0]
        observed_mean = self.aligned_reads / self.read_count_df.shape[0]
        ax2.set_title(
            f"Sequencing depth\nmean: {round(observed_mean, 2)}, median: {round(my_median, 2)}"
        )

        ax2.boxplot(
            self.read_count_df[self.sample_name].apply(np.log10),
            medianprops={"linestyle": "-", "linewidth": 1, "color": "black"},
            showfliers=False,
        )
        ax2.scatter(
            np.random.normal(1, 0.02, self.read_count_df.shape[0]),
            self.read_count_df[self.sample_name].apply(np.log10),
            color="lightgrey",
            alpha=0.4,
        )

        ax2.axhline(
            np.log10(expected_mean),
            color="lightgrey",
            linestyle="--",
            zorder=0,
            label=f"Expected mean: {round(expected_mean, 2)}",
        )
        ax2.axhline(
            np.log10(observed_mean),
            color="red" if observed_mean < expected_mean else "green",
            linestyle="--",
            zorder=12,
            label=f"Observed mean: {round(observed_mean, 2)}",
        )

        ax2.set_ylabel("log10(read count + 1)")
        ax2.set_xticklabels([self.sample_name])
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2_legend = ax2.legend(loc="best", frameon=False)
        for legend_handle in ax2_legend.legendHandles:
            legend_handle.set_alpha(1)

        return ax2

    def _dist_plot(self, fig, grid_spec):
        """
        Generate a plot representing the distribution of read counts, highlights abundance of template sequences.

        Parameters
        ----------
        fig : plt.Figure
            A matplotlib figure object
        grid_spec : plt.gridspec.GridSpect
            A matplotlib GridSpec object.

        Returns
        -------
        ax3 : plt.axes.Axes
            A matplitlib Axes object.
        """
        ax3 = fig.add_subplot(grid_spec[0:4, 8:16])
        ax3.set_title("Read counts")
        ax3.set_xlabel(
            f"gRNA{'s' if self.sample_type == 'single' else ' combinations' if self.sample_type == 'paired' else ''}"
        )
        ax3.set_ylabel("log10(read count + 1)")
        ax3.set_yscale("log")
        ax3.spines["top"].set_visible(False)
        ax3.spines["right"].set_visible(False)

        self.read_count_df.sort_values(by=self.sample_name, ascending=True, inplace=True)
        self.read_count_df["i"] = range(0, self.read_count_df.shape[0])

        if self.sample_type == "single":
            df_excl_template = self.read_count_df[~self.read_count_df.index.str.contains(self.template_marker)]
            df_only_template = self.read_count_df[self.read_count_df.index.str.contains(self.template_marker)]
            template_label = "template"
        elif self.sample_type == "paired":
            df_excl_template = self.read_count_df[~(self.read_count_df.index.get_level_values(0).str.contains(self.template_marker)) & ~(self.read_count_df.index.get_level_values(1).str.contains(self.template_marker))]
            df_only_template = self.read_count_df[(self.read_count_df.index.get_level_values(0).str.contains(self.template_marker)) & (self.read_count_df.index.get_level_values(1).str.contains(self.template_marker))]
            template_label = "template-template"

        ax3.scatter(
            df_excl_template["i"],
            df_excl_template[self.sample_name],
            color="royalblue", alpha=0.3
        )
        ax3.scatter(
            df_only_template["i"],
            df_only_template[self.sample_name],
            color="magenta", alpha=0.7, label=template_label
        )
        if df_only_template.shape[0] > 0:
            ax3_legend = ax3.legend(loc="upper left", frameon=False)
            for legend_handle in ax3_legend.legendHandles:
                legend_handle.set_alpha(1)

        return ax3

    def _cumul_plot(self, fig, grid_spec):
        """
        Generate a cumulative plot representing the distribution of read counts, provides the area under the curve to assess uniformity.
        Provides the number of found gRNAs.

        Parameters
        ----------
        fig : plt.Figure
            A matplotlib figure object
        grid_spec : plt.gridspec.GridSpect
            A matplotlib GridSpec object.

        Returns
        -------
        ax4 : plt.axes.Axes
            A matplitlib Axes object.
        """
        ax4 = fig.add_subplot(grid_spec[4:8, 0:8])

        self.read_count_df["ideal"] = 1

        self.read_count_df.sort_values(by=self.sample_name, ascending=False, inplace=True)

        for col in [self.sample_name, "ideal"]:
            x_values = [x / self.read_count_df.shape[0] for x in range(self.read_count_df.shape[0])]
            y_cumul = np.cumsum(self.read_count_df[col] / self.read_count_df[col].sum())
            auc = np.trapz(
                y_cumul,
                [i * (1.0 / self.read_count_df.shape[0]) for i in range(self.read_count_df.shape[0])],
            )

            x_values.append(1.0)
            y_values = list(
                pd.concat(
                    [
                        pd.Series([0.0]),
                        np.cumsum(self.read_count_df[col] / self.read_count_df[col].sum()),
                    ]
                ).values
            )

            if col == "ideal":
                my_color = "dimgray"
            else:
                my_color = "royalblue"

            ax4.plot(x_values, y_values, label=f"{col}, AUC={round(auc, 2)}", color=my_color)

        ax4.set_xlabel("fraction of NGS reads, ranked by abundance")
        ax4.set_ylabel("cumulative fraction of NGS reads")

        ax4.legend(frameon=False)
        found = self.read_count_df[self.read_count_df[self.sample_name] > 0].shape[0]
        total = self.read_count_df.shape[0]
        ax4.set_title(
            f"Cumulative distribution of read counts (Lorenz plot)\nfound {found:n} of {total:n} ({round(found / (total / 100), 2):n}% complete)"
        )
        ax4.spines["top"].set_visible(False)
        ax4.spines["right"].set_visible(False)

        return ax4

    def _histogram(self, fig, grid_spec):
        """
        Generate a histogram and a density plot representing the distribution of read counts, shows the skew of the distribution.

        Parameters
        ----------
        fig : plt.Figure
            A matplotlib figure object
        grid_spec : plt.gridspect.GridSpect
            A matplotlib GridSpec object.

        Returns
        -------
        ax5 : plt.axes.Axes
            A matplitlib Axes object.
        """
        ax5 = fig.add_subplot(grid_spec[4:8, 8:16])

        self.read_count_df.sort_values(by=self.sample_name, ascending=True, inplace=True)

        top10 = np.percentile(self.read_count_df[self.sample_name], 90)
        bottom10 = np.percentile(self.read_count_df[self.sample_name], 10)
        if bottom10 != 0:
            skew_ratio = top10 / bottom10
        else:
            skew_ratio = np.nan
        ax5.set_title(
            f"Histogram of read counts,\nlibrary distribution skew = {round(skew_ratio, 2)}"
        )

        bins = int(
            round(1 + 3.322 * np.log10(self.read_count_df.shape[0]))
        )  # determine number of bins using Sturge's rule

        sns.histplot(data=self.read_count_df, x=self.sample_name, kde=True, bins=bins, ax=ax5)
        #self.read_count_df[self.sample_name].hist(bins=b, grid=False, ax=ax5)
        plt.draw()

        ax5.axvline(bottom10, linestyle="--")
        ax5.axvline(top10, linestyle="--", label="10% & 90%")

        ax5.spines["top"].set_visible(False)
        ax5.spines["right"].set_visible(False)

        ax5_legend = ax5.legend(loc="best", frameon=False)
        for legend_handle in ax5_legend.legendHandles:
            legend_handle.set_alpha(1)

        return ax5

    def __repr__(self):
        """Representation."""
        return f"{self.__class__.__name__}"

    def __str__(self):
        """String representation."""
        return f"<{self.__class__.__name__}>"
