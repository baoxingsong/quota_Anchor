import matplotlib.pyplot as plt
import logging
import datetime
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from collections import OrderedDict
import numpy as np
from math import pi
from . import base
import sys
import seaborn as sns
from alive_progress import alive_bar

SPECIES_X = -0.1
BEZIER_RATIO = 0.316
DPI = 300
GAP_RATIO = 0.05
CAP_DIAMETER = 0.015
logger = logging.getLogger('main.line')
class Line:
    def __init__(self, config_pra, parameter):
        self.overwrite = False
        self.remove_chromosome_prefix = "chr,CHR,Chr"
        self.chr_font_size = 7
        self.species_name_font_size = 7

        self.input_file = ""
        self.length_file = ""
        self.species_name = ""
        self.output_file_name = ""
        self.figsize = "14,14"

        for i in config_pra.sections():
            if i == 'line':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.chr_font_size = float(self.chr_font_size)
        self.species_name_font_size = float(self.species_name_font_size)
        # 0.6 = 0.3+0.3
        self.ref_height = 0.3
        self.height_gap = 0.3
        self.query_height = 0.6

    @staticmethod
    def circle(radius, x):
        y = np.sqrt(radius ** 2 - x ** 2)
        return y

    @staticmethod
    def line_get_pos_list(df_dup, gap_length):
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                start_list.append(gap_length)
                end_list.append(gap_length + lgh)
                chr_to_start[chr_name[i]] = gap_length
                i += 1
            else:
                start_list.append(end_list[i-1] + gap_length)
                end_list.append(start_list[i] + lgh)
                chr_to_start[chr_name[i]] = start_list[-1]
                i += 1
        return start_list, end_list, chr_to_start

    @staticmethod
    def set_palette():
        colors = ['#FFA07A', '#4b6870', '#4169E1', '#9370DB']
        cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
        return cmap

    @staticmethod
    def set_husl_palette():
        colors = sns.husl_palette(as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("husl_256", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def plot_line_chr(start, end, cap_diameter, height):
        new_start = start + cap_diameter / 2
        new_end = end - cap_diameter / 2

        # right
        t = np.arange(-pi/2, pi/2, pi / 180)
        x = list(cap_diameter / 2 * np.cos(t) + new_end)
        y = list(cap_diameter / 2 * np.sin(t) + height + cap_diameter / 2)

        # top
        x.append(new_end)
        y.append(height + cap_diameter)
        x.append(new_start)
        y.append(height + cap_diameter)

        # left
        t = np.arange(pi/2, 3*pi/2, pi / 180)
        x1 = list(cap_diameter / 2 * np.cos(t) + new_start)
        y1 = list(cap_diameter / 2 * np.sin(t) + height + cap_diameter / 2)
        x += x1
        y += y1

        # bottom
        x.append(new_start)
        y.append(height)
        x.append(new_end)
        y.append(height)

        return x, y

    def plot_collinearity_region(self, pos1, pos2, pos3, pos4, query_height, ref_height, jg_qr_st, jg_qr_ed, jg_rf_st, jg_rf_ed):
        # ref_mar = ref_height + CAP_DIAMETER * 1.1
        # query_mar = query_height - CAP_DIAMETER * 0.1
        ref_mar = ref_height + CAP_DIAMETER * 1.0
        query_mar = query_height
        x, y = [], []
        # p1->p3(ref_block)
        t = np.linspace(pos1, pos3, 1000)
        for i in t:
            if jg_rf_st - CAP_DIAMETER / 2 <= i < jg_rf_st:
                will_x = i
                will_y = ref_height + CAP_DIAMETER / 2 + self.circle(CAP_DIAMETER / 2, jg_rf_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_rf_st <= i <= jg_rf_ed:
                will_x = i
                will_y = ref_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = ref_height + CAP_DIAMETER / 2 + self.circle(CAP_DIAMETER / 2, i - jg_rf_ed)
                x.append(will_x)
                y.append(will_y)
        # x.append(pos1)
        # y.append(ref_mar)
        # x.append(pos3)
        # y.append(ref_mar)

        # p3->p4
        pos3_x = x[-1]
        pos3_y = y[-1]
        if jg_qr_st - CAP_DIAMETER / 2 < pos4 < jg_qr_st:
            will_x = pos4
            will_y = query_height + CAP_DIAMETER / 2 - self.circle(CAP_DIAMETER / 2, jg_qr_st - pos4)
        elif jg_qr_st <= pos4 <= jg_qr_ed:
            will_x = pos4
            will_y = query_mar
        else:
            will_x = pos4
            will_y = query_height + CAP_DIAMETER / 2 - self.circle(CAP_DIAMETER / 2, pos4 - jg_qr_ed)
        t = np.arange(0, 1.01, 0.005)
        dx = will_x - pos3_x
        dy = will_y - pos3_y
        p1, p2, p3, p4 = pos3_x, dx * BEZIER_RATIO + pos3_x, -dx * BEZIER_RATIO + will_x, will_x
        x1 = base.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos3_y, (1 - BEZIER_RATIO) * dy + pos3_y, -(1-BEZIER_RATIO) * dy + will_y, will_y
        y1 = base.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        # p4->p2
        t = np.linspace(pos4, pos2, 1000)
        for i in t:
            if jg_qr_st - CAP_DIAMETER / 2 < i < jg_qr_st:
                will_x = i
                will_y = query_height + CAP_DIAMETER / 2 - self.circle(CAP_DIAMETER / 2, jg_qr_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_qr_st <= i <= jg_qr_ed:
                will_x = i
                will_y = query_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = query_height + CAP_DIAMETER / 2 - self.circle(CAP_DIAMETER / 2, i - jg_qr_ed)
                x.append(will_x)
                y.append(will_y)

        # p2->p1
        pos2_x = x[-1]
        pos2_y = y[-1]
        pos1_x = x[0]
        pos1_y = y[0]
        t = np.arange(0, 1.01, 0.005)
        dx = pos1_x - pos2_x
        dy = pos1_y - pos2_y
        p1, p2, p3, p4 = pos2_x, dx * BEZIER_RATIO + pos2_x, -dx * BEZIER_RATIO + pos1_x, pos1_x
        x1 = base.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos2_y, (1 - BEZIER_RATIO) * dy + pos2_y, -(1-BEZIER_RATIO) * dy + pos1_y, pos1_y
        y1 = base.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        return x, y

    def sub_run(self, collinearity, prefix, ref_length, query_length, chr_plus_gap_length, query_height, ref_height, loop, last_loop, strip_chr_abbr):
        logger.info(f"plot {collinearity} start.")
        ref_chr_list = ref_length['chr'].tolist()
        query_chr_list = query_length['chr'].tolist()
        chr_color_dict = {}
        cmap = self.set_palette()
        i = 1
        for ch in ref_chr_list:
            color_dict = {'chr': cmap(round(i / len(ref_chr_list), 2))}
            chr_color_dict[ch] = color_dict
            i += 1
        # figure geometry info
        ref_gap_length = (chr_plus_gap_length - ref_length['length'].sum()) / (len(ref_chr_list) + 1)
        ref_start_list, ref_end_list, ref_chr_to_start = self.line_get_pos_list(ref_length, ref_gap_length)
        # 20240818 only_for_chr_margin_collinearity_plot
        only_for_chr_margin_collinearity_plot_ref_zip = zip(ref_chr_list, ref_start_list, ref_end_list)
        only_for_chr_margin_collinearity_plot_ref_chr_to_position = {}
        for ch, start, end in only_for_chr_margin_collinearity_plot_ref_zip:
            start_x = start / chr_plus_gap_length + CAP_DIAMETER / 2
            end_x = end / chr_plus_gap_length - CAP_DIAMETER / 2
            only_for_chr_margin_collinearity_plot_ref_chr_to_position[ch] = [start_x, end_x]

        # relative length 1
        species_y = ref_height + CAP_DIAMETER / 2
        plt.text(SPECIES_X, species_y, prefix[0], ha="center", va="center", fontsize=self.species_name_font_size, color='black')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / chr_plus_gap_length
            ref_end_x = ref_end_list[i] / chr_plus_gap_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, CAP_DIAMETER, ref_height)
            plt.fill(x, y, facecolor='white', alpha=.7, edgecolor='black', linewidth=1)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = ref_height + CAP_DIAMETER / 2
            text_chr = ref_chr_list[i][len(prefix[0]):]
            for abbr in strip_chr_abbr:
                if text_chr.startswith(abbr):
                    text_chr = text_chr[len(abbr):]
            plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size, color='black')

        query_gap_length = (chr_plus_gap_length - query_length['length'].sum()) / (len(query_chr_list) + 1)
        query_start_list, query_end_list, query_chr_to_start = self.line_get_pos_list(query_length, query_gap_length)
        # 20240818 only_for_chr_margin_collinearity_plot
        only_for_chr_margin_collinearity_plot_query_zip = zip(query_chr_list, query_start_list, query_end_list)
        only_for_chr_margin_collinearity_plot_query_chr_to_position = {}
        for ch, start, end in only_for_chr_margin_collinearity_plot_query_zip:
            start_x = start / chr_plus_gap_length + CAP_DIAMETER / 2
            end_x = end / chr_plus_gap_length - CAP_DIAMETER / 2
            only_for_chr_margin_collinearity_plot_query_chr_to_position[ch] = [start_x, end_x]
        if loop == last_loop:
            species_y = query_height + CAP_DIAMETER / 2
            plt.text(SPECIES_X, species_y, prefix[1], ha="center", va="center", fontsize=self.species_name_font_size, color='black')
            for i in range(len(query_chr_list)):
                query_start_x = query_start_list[i] / chr_plus_gap_length
                query_end_x = query_end_list[i] / chr_plus_gap_length
                x, y = self.plot_line_chr(query_start_x, query_end_x, CAP_DIAMETER, query_height)
                plt.fill(x, y, facecolor='white', alpha=0.7, edgecolor='black', linewidth=1)
                label_x = (query_start_x + query_end_x) / 2
                label_y = query_height + CAP_DIAMETER / 2
                text_chr = query_chr_list[i][len(prefix[1]):]
                for abbr in strip_chr_abbr:
                    if text_chr.startswith(abbr):
                        text_chr = text_chr[len(abbr):]
                plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size, color='black')
        chr_list = list(OrderedDict.fromkeys(query_chr_list + ref_chr_list))
        chr_to_start = {}
        chr_to_start.update(query_chr_to_start)
        chr_to_start.update(ref_chr_to_start)
        data, gn_to_pos, rf_blk_chr, qry_blk_chr = base.read_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)
        intra = []
        intra_chr_list = []
        i = 0
        dt = datetime.datetime.now()
        time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
        with alive_bar(len(data), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
            for block in data:
                if qry_blk_chr[i] == rf_blk_chr[i]:
                    intra.append(block)
                    intra_chr_list.append(qry_blk_chr[i])
                    i += 1
                else:
                    color = chr_color_dict[rf_blk_chr[i]]['chr']
                    id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                    pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
                    pos1_coord_x = pos1 / chr_plus_gap_length
                    pos2_coord_x = pos2 / chr_plus_gap_length
                    pos3_coord_x = pos3 / chr_plus_gap_length
                    pos4_coord_x = pos4 / chr_plus_gap_length
                    # 20240818 only_for_chr_margin_collinearity_plot
                    query_chr = qry_blk_chr[i]
                    ref_chr = rf_blk_chr[i]
                    judge_fake_query_start_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][0]
                    judge_fake_query_end_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][1]
                    judge_fake_ref_start_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][0]
                    judge_fake_ref_end_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][1]
                    x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x,
                                                         query_height, ref_height,
                                                         judge_fake_query_start_x, judge_fake_query_end_x,
                                                         judge_fake_ref_start_x, judge_fake_ref_end_x)
                    plt.fill(x, y, facecolor=color, alpha=0.7)
                    i += 1
                bar()
        if intra:
            i = 0
            dt = datetime.datetime.now()
            time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
            with alive_bar(len(intra), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
                for block in intra:
                    # color = cmap(round(intra_color_index[i] / len(data), 2))
                    color = chr_color_dict[intra_chr_list[i]]['chr']
                    id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                    pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
                    pos1_coord_x = pos1 / chr_plus_gap_length
                    pos2_coord_x = pos2 / chr_plus_gap_length
                    pos3_coord_x = pos3 / chr_plus_gap_length
                    pos4_coord_x = pos4 / chr_plus_gap_length

                    # 20240818 only_for_chr_margin_collinearity_plot
                    query_chr = ref_chr = intra_chr_list[i]
                    judge_fake_query_start_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][0]
                    judge_fake_query_end_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][1]
                    judge_fake_ref_start_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][0]
                    judge_fake_ref_end_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][1]
                    x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                         judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
                    plt.fill(x, y, facecolor=color, alpha=0.8)
                    i += 1
                    bar()
        logger.info(f"plot {collinearity} finished.")

    def line_init(self):
        print()
        for key, value in vars(self).items():
            if key != "conf" and key not in ["ref_height", "height_gap", "query_height"]:
                print(key, "=", value)
        print()

        if not self.length_file:
            logger.error("Please specify your chromosome length file(Separator: ',')")
            sys.exit(1)

        if not self.species_name:
            logger.error("Please specify your species name(Separator: ',')")
            sys.exit(1)
    @staticmethod
    def split_file(file_string, separator):
        file_list = file_string.split(separator)
        new_file_list = []
        for file in file_list:
            if len(file) == 0:
                continue
            file = file.strip()
            base.file_empty(file)
            new_file_list.append(file)
        return new_file_list

    def split_length(self, strip_prefix):
        length_file_list = self.length_file.split(",")
        strip_length_file_list_df = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            base.file_empty(le)
            length = pd.read_csv(le, sep='\t', header=0, index_col=None)
            length['chr'] = length['chr'].astype(str)
            length['chr'] = strip_prefix[i] + length['chr']
            length['length'] = length['length'].astype(int)
            strip_length_file_list_df.append(length)
            i += 1
        return strip_length_file_list_df

    @staticmethod
    def get_max_chr_plus_gap_length(strip_length_file_list_df):
        length_list = []
        for df in strip_length_file_list_df:
            length_list.append(df['length'].sum())
        total_chr_length = max(length_list)
        chr_plus_gap_length = total_chr_length + total_chr_length * GAP_RATIO
        return chr_plus_gap_length

    def run(self):
        logger.info("Init line and the following parameters are config information.")
        self.line_init()
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)

        strip_chr_abbr = base.split_conf(self.remove_chromosome_prefix, ",")
        strip_collinearity_file_list = self.split_file(self.input_file, ",")
        strip_prefix = base.split_conf(self.species_name, ",")
        strip_length_file_list_df = self.split_length(strip_prefix)

        pair_number = len(strip_collinearity_file_list)
        new_prefix = []
        for i in range(pair_number):
            new_prefix.append([strip_prefix[i], strip_prefix[i + 1]])

        new_length_file_list_df = []
        for i in range(pair_number):
            new_length_file_list_df.append([strip_length_file_list_df[i], strip_length_file_list_df[i+1]])

        chr_plus_gap_length = self.get_max_chr_plus_gap_length(strip_length_file_list_df)
        zipped_three_pair = list(zip(strip_collinearity_file_list, new_prefix, new_length_file_list_df))

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='white')
        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1

        i = 0
        for col, prefix, length in zipped_three_pair:
            self.sub_run(col, prefix, length[0], length[1], chr_plus_gap_length, query_height, ref_height, i, last_loop, strip_chr_abbr)
            query_height = query_height + self.height_gap
            ref_height = ref_height + self.height_gap
            i += 1

        plt.axis('off')
        # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.output_file_name, dpi=DPI, bbox_inches='tight')
        logger.info(f"generate {self.output_file_name} finished.")
        sys.exit(0)
