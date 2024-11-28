import re, sys
import logging
from ..lib import base
import random
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, median_abs_deviation
from sklearn.mixture import GaussianMixture

logger = logging.getLogger('main.ks_fitting')
class Kf:
    def __init__(self, config_pra, parameter):
        self.ref_length = ""
        self.query_length = ""
        self.ks_file = ""
        self.collinearity_file = ""
        self.output_file = ""
        self.overwrite = False
        self.components = 4
        self.correct_file = ""
        self.focal_species = ""
        self.disable_arrow = True

        self.ks_range = "0,3"
        self.figsize = "12,7"
        for i in config_pra.sections():
            if i == 'ks_fitting':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])

        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.minus_axis_ratio = 0.1
        self.move_distance = 0
        self.bins_number = 100
        self.components = int(self.components)
        self.col = "ks_NG86"
        self.min_block_length = 10

    def fitting_init(self):
        if not self.query_length:
            logger.error("Please specify your query species chromosome length file")
            sys.exit(1)
        if not self.ref_length:
            logger.error("Please specify your reference species chromosome length file")
            sys.exit(1)
        if not self.ks_file:
            logger.error("Please specify your focal species ks file")
            sys.exit(1)
        if not self.collinearity_file:
            logger.error("Please specify your focal species collinearity file(self vs self)")
            sys.exit(1)
        if not self.output_file:
            logger.error("Please specify your output file name")
            sys.exit(1)

        base.file_empty(self.ref_length)
        base.file_empty(self.query_length)
        base.file_empty(self.ks_file)
        base.file_empty(self.collinearity_file)
        base.output_file_parentdir_exist(self.output_file, self.overwrite)

    @staticmethod
    def read_ks_get_map(ks_file, col, ks_min, ks_max):
        df = pd.read_csv(ks_file, sep='\t', header=0, index_col=None, comment="#")
        df = df[(df[col] >= max(ks_min, 0)) & (df[col] <= ks_max) & (df[col] != -1)]
        df['id1'] = df['id1'].astype(str)
        df['id2'] = df['id2'].astype(str)
        df['key'] = df['id1'] + '_' + df['id2']
        df['reverse_key'] = df['id2'] + '_' + df['id1']
        df_keys = pd.concat([df[['key', col]], df[['reverse_key', col]].rename(columns={'reverse_key': 'key'})])
        df_keys.drop_duplicates(subset=["key"] ,inplace=True)
        dict_pair_ks = df_keys.set_index('key')[col].to_dict()

        return dict_pair_ks
    
    @staticmethod
    def read_collinearity(collinearity, ref_chr_list, query_chr_list):
        collinearity_dict = {}
        flag = False
        with open(collinearity)as f:
            _ = next(f)
            _ = next(f)
            for line in f:
                if line.startswith("#"):
                    flag = True
                    chr_pair = line.split()[4].split("&")
                    ref_chr = chr_pair[0]
                    query_chr = chr_pair[1]
                    if ref_chr not in ref_chr_list or query_chr not in query_chr_list:
                        flag = False
                        continue
                    block_number = str(int(line.split()[1]) - 1)
                    collinearity_dict[block_number] = []
                    continue
                else:
                    if flag:
                        pair_name = line.split()[0] + "_" + line.split()[5]
                        collinearity_dict[block_number].append(pair_name)

        return collinearity_dict

    @staticmethod
    def get_block_all_ks_and_median_and_no_outlier(ks_dict, collinearity_dict, min_block_length=5):
        block_dict = {}
        block_dict_rm_outlier = {}
        median_one = {}
        median_one_rm_outlier = {}
        last_key = -1
        for key, values in collinearity_dict.items():
            block_dict[key] = []
            if last_key != -1:
                if len(block_dict[last_key]) >= min_block_length:
                    # init
                    block_dict_rm_outlier[last_key] = []
                    median_one[last_key] = np.median(block_dict[last_key])
                    # remove outliers
                    mad = median_abs_deviation(block_dict[last_key])
                    for i in block_dict[last_key]:
                        if median_one[last_key] - mad <= i <= median_one[last_key] + mad: # type: ignore
                            block_dict_rm_outlier[last_key].append(i)
                    median_one_rm_outlier[last_key] = np.median(block_dict_rm_outlier[last_key])
                else:
                    block_dict.pop(last_key)
            for pair_name in values:
                try:
                    ks = ks_dict[pair_name]
                except KeyError:
                    # ks==-1 had been deleted
                    continue
                block_dict[key].append(ks)
            last_key = key
        if len(block_dict[last_key]) >= min_block_length:
            # init
            block_dict_rm_outlier[last_key] = []
            median_one[last_key] = np.median(block_dict[last_key])
            # remove outliers
            mad = median_abs_deviation(block_dict[last_key])
            for i in block_dict[last_key]:
                if median_one[last_key] - mad <= i <= median_one[last_key] + mad:
                    block_dict_rm_outlier[last_key].append(i)
            median_one_rm_outlier[last_key] = np.median(block_dict_rm_outlier[last_key])
        else:
            block_dict.pop(last_key)
        return block_dict_rm_outlier, median_one_rm_outlier

    def get_all_median_dict(self, ks_min, ks_max, ref_chr_list, query_chr_list):
        ks_dict = self.read_ks_get_map(self.ks_file, self.col, ks_min, ks_max)
        collinearity_dict = self.read_collinearity(self.collinearity_file, ref_chr_list, query_chr_list)
        block_dict_ks, block_dict_median = self.get_block_all_ks_and_median_and_no_outlier(ks_dict, collinearity_dict, self.min_block_length)
        
        return block_dict_ks, block_dict_median

    @staticmethod
    def gmm(all_medians_list, n_wgds, max_iter=618, n_init=1):
        all_medians_list_logtransformed = np.log(all_medians_list)
        # all_medians_list = np.array(all_medians_list)
        gmm = GaussianMixture(n_components=n_wgds, covariance_type="spherical", max_iter=max_iter, n_init=n_init)
        gmm.fit(all_medians_list_logtransformed.reshape(len(all_medians_list_logtransformed), 1))
        labels = gmm.predict(all_medians_list_logtransformed.reshape(len(all_medians_list_logtransformed), 1))
        return labels

    @staticmethod
    def get_cluster_color(cluster_to_all_ks):
        cluster_color = {}
        cluster_median = []
        for cluster in cluster_to_all_ks:
            cluster_median.append([cluster, np.median(cluster_to_all_ks[cluster])])
        sorted_cluster_median = sorted(cluster_median, key=lambda x: x[1])
        color_list = [("royalblue", "a"), ("red", "b"), ("green", "c"), ("black", "d"), ("yellow", "e")]
        for i in range(len(sorted_cluster_median)):
            cluster_color[sorted_cluster_median[i][0]] = color_list[i]
        return cluster_color

    def prepare_cluster_ks_info(self, block_label, block_dict_ks, block_dict_median, all_median_list):
        cluster_to_median = {}
        cluster_to_alignment_median_pair = {}
        cluster_to_all_ks = {}
        # init dict list
        cluster_set = set(block_label)
        for cluster in cluster_set:
            cluster_to_median[cluster] = []
            cluster_to_alignment_median_pair[cluster] = []
            cluster_to_all_ks[cluster] = []

        block_name_list = list(block_dict_median.keys())
        for i in range(len(block_label)):
            cluster_bk = block_label[i]
            block_name = block_name_list[i]
            block_list = block_dict_ks[block_name]
            median_bk = all_median_list[i]
            cluster_to_median[cluster_bk].append(median_bk)        
            cluster_to_alignment_median_pair[cluster_bk].append([block_name, median_bk])
            cluster_to_all_ks[cluster_bk].extend(block_list)
        # get cluster_color
        cluster_color = self.get_cluster_color(cluster_to_all_ks)
        return cluster_to_median, cluster_to_alignment_median_pair, cluster_to_all_ks, cluster_color

    @staticmethod
    def gaussian_fuc(x, center, sd):
        y = np.zeros_like(x) + 1 /(np.sqrt(2 * np.pi) * sd) * np.exp(-(x - center)**2/(2*sd**2))
        return y

    @staticmethod
    def gaussian_approximate_fuc(x, head, center, fake_sd):
        y = np.zeros_like(x) + head * np.exp(-((x - center)/fake_sd)**2)
        return y

    def kde_fit(self, data, x):
        kde = gaussian_kde(data, bw_method='scott')
        kde.set_bandwidth(bw_method=kde.factor/3.)
        p = kde(x)

        result = curve_fit(self.gaussian_approximate_fuc, x, p, maxfev=102400)
        parameter = result[0]
        y = self.gaussian_approximate_fuc(x, *parameter)

        return y, parameter

    @staticmethod
    def determine_marker_position(cluster_id, paras, cluster_info):
        head = paras[0]
        top_x = paras[1]

        top_y = head

        cluster_info[cluster_id] = [top_x, top_y]
        return cluster_info

    def determine_corrdinate_gap1(self, fig: plt.figure, y_max: float, div_peaks_number: int):
        # height unit inch
        height = fig.get_size_inches()[1]
        # default bottom 0.11  top 0.88
        bottom_ratio = fig.subplotpars.bottom
        top_ratio = fig.subplotpars.top
        # retain points
        retain_point = (top_ratio - bottom_ratio) * height * 72
        # actual_tick/point
        unit = y_max / retain_point
        # move actual tick
        move_tick = unit * self.move_distance
        move_gap_tick = round((move_tick / (div_peaks_number + 1)), 2)
        
        return move_tick, move_gap_tick
    
    def determine_coordinate_gap(self, y_max: float, div_peaks_number: int):
        
        move_tick = self.minus_axis_ratio * y_max
        move_gap_tick = round((move_tick / (div_peaks_number + 1)), 2)
        return move_tick, move_gap_tick

    def plot_wgd_marker(self, ax, cluster_id, cluster_color, cluster_info, y_max_limit, y_min_ratio):

        top_x = cluster_info[cluster_id][0]
        top_y = cluster_info[cluster_id][1]
        marker_color = cluster_color[cluster_id][0]
        marker_name = cluster_color[cluster_id][1]
        
        zorder = cluster_info[cluster_id][2]

        if self.correct_file:
            line_max = y_max_limit * ( 1 + self.minus_axis_ratio)
            line_ratio = (top_y + self.minus_axis_ratio * y_max_limit) / line_max
        else:
            line_max = y_max_limit
            line_ratio = top_y / line_max
        ax.axvline(x=top_x, ymin=y_min_ratio, ymax=line_ratio + 0.04, color=marker_color,
                 linestyle=(0, (3, 3)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                 zorder=zorder +1)

        # s= 256 = 16 point  = 16 linewidth = 16 marker size
        # (16 + 1 * 0.5) ** 2

        ax.scatter(x=top_x, y=top_y + 0.04 * line_max, s=256,
                    facecolor="white", edgecolor="white",
                    linewidth=1, zorder=zorder +2)

        ax.scatter(x=top_x, y=top_y + 0.04 * line_max, s=256,
                    facecolor=(0.64, 0.64, 0.64, 0.75), edgecolor="w",
                    linewidth=1, zorder=zorder +3)

        ax.scatter(x=top_x, y=top_y + 0.04 * line_max, s=256,
                    facecolor=to_rgba(marker_color, 0.4), edgecolor="w",
                    linewidth=1, zorder=zorder +4)
        ax.text(x=top_x, y=top_y + 0.04 * line_max, s=marker_name, fontsize=12,
                horizontalalignment='center', verticalalignment='center', clip_on=True, zorder=zorder +5, color='w')
    
    def determine_y_min(self):
        if self.correct_file and not self.disable_arrow:
            y_min_ratio = self.minus_axis_ratio / (1 + self.minus_axis_ratio)
        else:
            y_min_ratio = 0
        return y_min_ratio

    @staticmethod
    def sort_key_mode(item):
        match = re.search(r'(\d+\.\d+)', item[1])
        if match:
            return float(match.group(1))
        return 0.0

    def focal_species_first_run(self):
        self.fitting_init()
        ref_chr_list = pd.read_csv(self.ref_length, sep="\t", header=0, index_col=None)['chr'].astype(str).tolist()
        query_chr_list = pd.read_csv(self.query_length, sep="\t", header=0, index_col=None)['chr'].astype(str).tolist()

        ks_min, ks_max = float(self.ks_range.split(",")[0]), float(self.ks_range.split(",")[1])
        # block_dict_ks: key is alignment number-1(string 0 based) and value is ks value list
        # block_dict_median: key is alignment number-1(string 0 based) and value is median ks value
        block_dict_ks, block_dict_median = self.get_all_median_dict(ks_min, ks_max, ref_chr_list, query_chr_list)

        plt.rcParams['ytick.major.pad'] = 0
        plt.rcParams['mathtext.default'] = 'regular'
        plt.rcParams['font.size'] = 14

        all_median_list = list(block_dict_median.values())
        block_label = self.gmm(all_median_list, self.components)
        new_block_label = block_label.astype(str).tolist()

        cluster_to_median, cluster_to_alignment_median_pair, cluster_to_all_ks, cluster_color = self.prepare_cluster_ks_info(new_block_label, block_dict_ks, block_dict_median, all_median_list)
        figsize = [float(k) for k in self.figsize.split(',')]
        fig, ax = plt.subplots(1, 1, figsize=(figsize[0], figsize[1]))
        ax.set_xlim(ks_min, ks_max)

        ax.set_xlabel(r'${K_s}$')
        ax.set_ylabel('Density')
        plt.setp(ax.yaxis.get_majorticklabels(), rotation=90, va='center')        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # determine y axis limit
        cluster_info = {}
        zorder = 100

        map_handle_label = []
        for cluster_id in cluster_to_all_ks:
            x = np.linspace(ks_min, ks_max, 1024)
            y, paras = self.kde_fit(cluster_to_all_ks[cluster_id], x)
            # plot line
            ax.plot(x, y, color=cluster_color[cluster_id][0], linewidth=2, zorder=zorder-9)
            ctr = float(paras[1])
            # wgd marker zorder > histogram zorder > curve fitting line zorder
            _, _, patches = ax.hist(np.array(cluster_to_all_ks[cluster_id]), linewidth=2, bins=self.bins_number, alpha=0.3,
                                    color=cluster_color[cluster_id][0], density=True, label=f"Cluster {cluster_color[cluster_id][1]} (mode {round(ctr, 3)})", zorder=zorder-8)
            map_handle_label.append((patches[0], f"Cluster {cluster_color[cluster_id][1]} (mode {round(ctr, 3)})"))
            cluster_info = self.determine_marker_position(cluster_id, paras, cluster_info)
            
            cluster_info[cluster_id].append(zorder-7)
            zorder -= 10

        y_min_ratio = self.determine_y_min()
        top_y_max= 0
        for cluster_id in cluster_to_all_ks:
            top_y = cluster_info[cluster_id][1]
            top_y_max = max([top_y, top_y_max])
        ax.set_ylim(0, top_y_max * 1.2)
        y_max_limit = top_y_max * 1.2

        for cluster_id in cluster_to_all_ks:
            self.plot_wgd_marker(ax, cluster_id, cluster_color, cluster_info, y_max_limit, y_min_ratio)

        map_handle_label.sort(key=self.sort_key_mode)
        handles, labels = zip(*map_handle_label)
        handles, labels = list(handles), list(labels)

        if not self.correct_file:
            ax.spines['bottom'].set_position(('outward', self.move_distance))
            ax.legend(handles=handles, labels=labels, frameon=False)
            fig.suptitle(f'Clustering of syntenic pairs from {self.focal_species}')
            fig.savefig(self.output_file, transparent=True, dpi=300)
            plt.show()
        return fig, ax, y_min_ratio, handles, labels

    @staticmethod
    def plot_divergence_line(axis, peak, sd, line_color, label, z_order, y_min_ratio, height):

        # transparent rectangle with width: peak +/- sd
        axis.axvspan(xmin=peak - sd, xmax=peak + sd, ymin=y_min_ratio, ymax=height, facecolor=line_color,
                    alpha=0.3, capstyle='butt', joinstyle='miter', zorder=z_order - 2)
        # dashed line in correspondence of the distribution peak
        axis.axvline(x=peak, ymin=y_min_ratio, ymax=height, color=line_color, alpha=0.618,
                    linestyle=(0, (6, 6)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                    zorder=z_order - 1, label=label)
        # additional transparent white dashed line on top of previous to lighten line on plot but not the lines in legend
        axis.axvline(x=peak, ymin=y_min_ratio, ymax=height, color='w', alpha=0.3,
                    linestyle=(0, (6, 6)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                    zorder=z_order)

    @staticmethod
    def plot_arrow(axis, uncorrected_peak, corrected_peak, arrow_y, line_color):
        # if abs(corrected_peak - uncorrected_peak) > 0.2:
        axis.annotate("", xy=(corrected_peak, arrow_y), xytext=(uncorrected_peak, arrow_y),
                    arrowprops=dict(arrowstyle="->,head_length=0.08,head_width=0.05", shrinkA=0, shrinkB=0,
                                    edgecolor=line_color, alpha=1,
                                    connectionstyle="arc3", linewidth=0.9))
        # else:
        #     axis.annotate("", xy=(corrected_peak, arrow_y), xytext=(uncorrected_peak, arrow_y),
        #                 arrowprops=dict(arrowstyle="->,head_length=0.08,head_width=0.05", shrinkA=0, shrinkB=0,
        #                                 edgecolor=line_color, linewidth=0.9))

    @staticmethod
    def plot_divergence_id(axis, divergence_id, x, y_max, n_ids, circle_color, z_order):
        # so to reduce the chance that circles overlap.
        y_jitter_max = 0.008
        if n_ids > 1:
            y_jitter = random.uniform(-y_jitter_max * n_ids, y_jitter_max * n_ids)
        else:
            y_jitter = 0

        # white circle with colored border and number inside

        y = y_max * (0.922 + y_jitter)
        axis.scatter(x=x, y=y, s=220, facecolor='w', edgecolor=circle_color,
                    linewidth=1, zorder=z_order - 1)
        axis.text(x=x, y=y, s=divergence_id, fontsize=12, horizontalalignment='center', verticalalignment='center',
                clip_on=True, zorder=z_order, color='k')
    
    @staticmethod
    def define_legend_size(axis):
        __, labels = axis.get_legend_handles_labels()
        current_legend_width = 0
        prev = 0
        for label in labels:
            # Get a legend width that can host the longest latin name of the list 
            if len(label) <= 53:
                current_legend_width = max((len(label) + 17) / 100, prev)
                prev = current_legend_width
            else:
                current_legend_width = 0.7
                break
        final_legend_size = (1.01, 0, current_legend_width, 1)
        return final_legend_size

    @staticmethod
    def create_legend(axis, legend_size, insert_pos, wgd_handles, wgd_labels):
        handles, labels = axis.get_legend_handles_labels()
        sorted_handles, sorted_labels = handles.copy(), labels.copy()
        sorted_handles, sorted_labels = sorted_handles[insert_pos:], sorted_labels[insert_pos:]

        empty_rect = mpatches.Rectangle((0, 0), 0, 0, fill=False, edgecolor='none', visible=False)

        wgd_handles.extend([empty_rect, empty_rect])
        wgd_labels.extend(["", "Divergence with:"])

        sorted_handles, sorted_labels = wgd_handles + sorted_handles, wgd_labels + sorted_labels
        lgd = axis.legend(sorted_handles, sorted_labels, loc="upper left", handlelength=1.618, facecolor="w", frameon=False, mode="expand",
                        bbox_to_anchor=legend_size)
        
        return lgd

    def save_mixed_plot(self, fig_corr, ax_corr, wgd_handles, wgd_labels, super_title):
        legend_size = Kf.define_legend_size(ax_corr)
        chart_box = ax_corr.get_position()
        ax_corr.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
        lgd = Kf.create_legend(ax_corr, legend_size, self.components, wgd_handles, wgd_labels)
        fig_corr.savefig(self.output_file, bbox_extra_artists=(lgd, super_title), transparent=True)
        
    def run(self):
        logger.info("Ks_fitting module init and the following parameters are config information")
        print()
        for key, value in vars(self).items():
            if key not in ["conf", "minus_axis_ratio", "move_distance", "col", "min_block_length", "bins_number"]:
                print(key, "=", value)
        print()
        fig, ax, y_min_ratio, wgd_handles, wgd_labels = self.focal_species_first_run()

        if self.correct_file:
            df = pd.read_csv(self.correct_file, header=0, index_col=None)
            df['Focal_Species'] = df['Focal_Species'].astype(str)
            df['Sister_Species'] = df['Sister_Species'].astype(str)
            divergent_peaks_number = len(df)
            cmap = plt.get_cmap('gist_rainbow')

            color_list = [cmap(i) for i in np.linspace(0, 1, divergent_peaks_number)]
            super_title = fig.suptitle("Rate-adjusted mixed " + r'${K_s}$' + f" distribution for {self.focal_species}", y=0.98)
            y_max = ax.get_ylim()[1]

            ax.spines['bottom'].set_position(('outward', self.move_distance))
            move_tick, gap_scale =self.determine_coordinate_gap(y_max, divergent_peaks_number)
            ax.set_ylim(-move_tick, y_max)
            ax.spines['left'].set_bounds(0, y_max)

            to_plot = df[["Node", "Focal_Species", "Sister_Species", "Adjusted_Mode_Mean", "Adjusted_Mode_Mean_SD", "Original_Mode"]]
            arrow_y = -gap_scale
            z_order = -10
            for index, row in to_plot.iterrows():
                per_node_number = len(to_plot[to_plot["Node"] == row["Node"]])
                z_order -= 2
                node, focal_species, sister_species, mode, sd, original = row.at["Node"], row.at["Focal_Species"], row.at["Sister_Species"], row.at["Adjusted_Mode_Mean"], row.at["Adjusted_Mode_Mean_SD"], row.at["Original_Mode"]
                color = color_list[index]
                if original < mode:
                    means_label = str(round(original, 2)) + r'$\rightarrow$' + str(round(mode, 2))
                else:
                    means_label = str(round(mode, 2)) + r'$\leftarrow$' + str(round(original, 2))
                label = f"({node}) {sister_species} ({means_label})"
                self.plot_divergence_line(ax, mode, sd, color, label, -300, y_min_ratio, height=0.88)
                self.plot_arrow(ax, original, mode, arrow_y, color)
                self.plot_divergence_id(ax, node, mode, y_max, per_node_number, color, z_order)
                arrow_y -= gap_scale
            if self.disable_arrow:
                ax.set_ylim(0, y_max)
            self.save_mixed_plot(fig, ax, wgd_handles, wgd_labels, super_title) # type: ignore
            plt.show()
        logger.info("Plot ks fitting info finished!")
