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
        self.hide_chr = False
        self.italic = False
        self.actual_len = False
        self.color_style = "four_colors"
        self.gap_style = "compact"
        self.sp_chr_color_comma_sep = ""
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

    def determine_fig_par(self, len_number):
        x = 0.7 / (6 * len_number - 5)
        self.height_gap = 5 * x
        self.query_height = 0.3 + 5 * x
        global CAP_DIAMETER
        CAP_DIAMETER = 0.0015

    @staticmethod
    def line_get_pos_list_compact(df_dup, gap_length, total_length, flag):
        tmp_start = (total_length - df_dup['length'].sum() - len(df_dup) * gap_length) / 2
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                if flag:
                    start_list.append(gap_length + tmp_start)
                    end_list.append(gap_length + tmp_start + lgh)
                    chr_to_start[chr_name[i]] = gap_length + tmp_start
                else:
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
    def read_collinearity(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start, ref_qry_ratio):
        ref_ratio = ref_qry_ratio[0]
        query_ratio = ref_qry_ratio[1]
        data = []
        ref_chr_list = []
        query_chr_list = []
        gene_pos_dict = {}
        block_index = 0
        block = []
        direction_list = []
        # print(qry_prefix)
        # print(ref_prefix)
        # print(collinearity)
        # print(chr_list)
        # print(chr_to_start)
        with open(collinearity) as f:
            _ = next(f)
            _ = next(f)
            flag = True
            for line in f:
                if line.startswith("#"):
                    if block:
                        data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                        block = []
                    chr_pair = line.split()[4].split("&")
                    if ref_prefix + chr_pair[0] not in chr_list or qry_prefix + chr_pair[1] not in chr_list:
                        flag = False
                    else:
                        flag = True
                        ref_chr = ref_prefix + chr_pair[0]
                        ref_chr_list.append(ref_chr)
                        query_chr = qry_prefix + chr_pair[1]
                        query_chr_list.append(query_chr)
                        direction = line.split()[5]
                        direction_list.append(direction)
                        block_index += 1
                else:
                    if flag:
                        line_list = line.split()
                        block.append([line_list[0], line_list[5]])
                        gene_pos_dict[line_list[0]] = chr_to_start[ref_chr] + int(line_list[4]) * ref_ratio
                        gene_pos_dict[line_list[5]] = chr_to_start[query_chr] + int(line_list[9]) * query_ratio
                    else:
                        continue
            if block:
                data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
        return data, gene_pos_dict, ref_chr_list, query_chr_list, direction_list

    def get_my_cmap(self):
        if self.color_style == "rainbow":
            cmap = self.set_rainbow()
        elif self.color_style == "husl":
            cmap = self.set_husl_palette()
        elif self.color_style == "four_colors":
            cmap = self.set_palette()
        else:
            cmap = self.set_husl_palette()
        return cmap

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
    def set_rainbow():
        colors = plt.get_cmap("gist_rainbow")
        cmap = LinearSegmentedColormap.from_list("rainbow", colors(np.linspace(0, 1, 256)))
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

    def loose_compact_qry_ref(self, chr_plus_gap_length, length_df, length_df_chr_list, max_gm_gap_length):
        if self.gap_style == "loose":
            gap_length = (chr_plus_gap_length - length_df['length'].sum()) / (len(length_df_chr_list) + 1)
            start_list, end_list, chr_to_start = self.line_get_pos_list_compact(length_df, gap_length, chr_plus_gap_length, flag=False)
        else:
            gap_length = max_gm_gap_length
            start_list, end_list, chr_to_start = self.line_get_pos_list_compact(length_df, gap_length, chr_plus_gap_length, flag=True)
        # 20240818 only_for_chr_margin_collinearity_plot
        cap_judge = zip(length_df_chr_list, start_list, end_list)
        cap_judge_chr = {}
        for ch, start, end in cap_judge:
            start_x = start / chr_plus_gap_length + CAP_DIAMETER / 2
            end_x = end / chr_plus_gap_length - CAP_DIAMETER / 2
            cap_judge_chr[ch] = [start_x, end_x]
        return cap_judge_chr, chr_to_start, start_list, end_list

    @staticmethod
    def delete_prefix(chr_name, sp_name, strip_chr_abbr):
        text_chr = chr_name[len(sp_name):]
        for abbr in strip_chr_abbr:
            if text_chr.startswith(abbr):
                text_chr = text_chr[len(abbr):]
        return text_chr

    def plot_text_sp(self,species_y, sp_name):
        if self.italic:
            plt.text(SPECIES_X, species_y, sp_name, ha="center", va="center", fontsize=self.species_name_font_size, color='black', fontstyle='italic')
        else:
            plt.text(SPECIES_X, species_y, sp_name, ha="center", va="center", fontsize=self.species_name_font_size, color='black')
        return plt

    def get_sp_chr_color_dict(self, sp_list):
        if self.sp_chr_color_comma_sep:
            color_list = base.split_conf(self.sp_chr_color_comma_sep, ",")
            i = 0
            while len(color_list) < len(sp_list):
                dup_color = color_list[i]
                color_list.append(dup_color)
                if i >= len(color_list)-1:
                    i = 0
                else:
                    i += 1
        else:
            color_list = ["black"] * len(sp_list)
        sp_chr_color_dict = dict(list(zip(sp_list, color_list)))
        return sp_chr_color_dict

    def plot_chr(self, chr_list, start_list, end_list, chr_plus_gap_length, height, sp_name, strip_chr_abbr, chr_color_dict):
        for i in range(len(chr_list)):
            start_x = start_list[i] / chr_plus_gap_length
            end_x = end_list[i] / chr_plus_gap_length
            x, y = self.plot_line_chr(start_x, end_x, CAP_DIAMETER, height)
            plt.fill(x, y, facecolor=chr_color_dict[sp_name], alpha=0.8, edgecolor=chr_color_dict[sp_name], linewidth=1, zorder=2)
            label_x = (start_x + end_x) / 2
            label_y = height + CAP_DIAMETER / 2
            text_chr = self.delete_prefix(chr_list[i], sp_name, strip_chr_abbr)
            if not self.hide_chr:
                plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size,
                         color='black', zorder=3)
            else:
                pass

    @staticmethod
    def prepare_coll(query_chr_list, ref_chr_list, query_chr_to_start, ref_chr_to_start):
        chr_list = list(OrderedDict.fromkeys(query_chr_list + ref_chr_list))
        chr_to_start = {}
        chr_to_start.update(query_chr_to_start)
        chr_to_start.update(ref_chr_to_start)
        return chr_list, chr_to_start

    def get_block_corr_color(self, color, direction, block, gn_to_pos, chr_plus_gap_length):
        if self.color_style == "two_colors" and direction == "POSITIVE":
            color = "#F0F0F0"
        if self.color_style == "two_colors" and direction == "NEGATIVE":
            color = "#66AD56"

        id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
        pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
        pos1_coord_x = pos1 / chr_plus_gap_length
        pos2_coord_x = pos2 / chr_plus_gap_length
        pos3_coord_x = pos3 / chr_plus_gap_length
        pos4_coord_x = pos4 / chr_plus_gap_length
        return pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, color

    @staticmethod
    def get_judge_se(cap_judge_chr, ch):
        judge_fake_start_x = cap_judge_chr[ch][0]
        judge_fake_end_x = cap_judge_chr[ch][1]
        return judge_fake_start_x, judge_fake_end_x

    def plot_block(self, cap_judge_query_chr, cap_judge_ref_chr, query_chr, ref_chr, pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height,
                                             ref_height, color, i):
        judge_fake_query_start_x, judge_fake_query_end_x = self.get_judge_se(cap_judge_query_chr, query_chr)
        judge_fake_ref_start_x, judge_fake_ref_end_x = self.get_judge_se(cap_judge_ref_chr, ref_chr)
        x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height, judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
        plt.fill(x, y, facecolor=color, alpha=0.5, zorder=1.5)
        i += 1
        return i

    def sub_run(self, collinearity, prefix, ref_length, query_length, chr_plus_gap_length, query_height, ref_height, loop, last_loop, strip_chr_abbr, max_gm_gap_length, sp_chr_color_dict):
        logger.info(f"Plot {collinearity} start.")
        ref_chr_list = ref_length['chr'].tolist()
        query_chr_list = query_length['chr'].tolist()
        chr_color_dict = {}
        cmap = self.get_my_cmap()

        i = 1
        for ch in ref_chr_list:
            color_dict = {'chr': cmap(round(i / len(ref_chr_list), 2))}
            chr_color_dict[ch] = color_dict
            i += 1

        # figure geometry info
        cap_judge_ref_chr, ref_chr_to_start, ref_start_list, ref_end_list = self.loose_compact_qry_ref(chr_plus_gap_length, ref_length, ref_chr_list, max_gm_gap_length)

        species_y = ref_height + CAP_DIAMETER / 2
        self.plot_text_sp(species_y, prefix[0])
        self.plot_chr(ref_chr_list, ref_start_list, ref_end_list, chr_plus_gap_length, ref_height, prefix[0], strip_chr_abbr, sp_chr_color_dict)

        cap_judge_query_chr, query_chr_to_start, query_start_list, query_end_list = self.loose_compact_qry_ref(chr_plus_gap_length, query_length, query_chr_list, max_gm_gap_length)
        if loop == last_loop:
            species_y = query_height + CAP_DIAMETER / 2
            self.plot_text_sp(species_y, prefix[1])
            self.plot_chr(query_chr_list, query_start_list, query_end_list, chr_plus_gap_length, query_height, prefix[1], strip_chr_abbr, sp_chr_color_dict)
        chr_list, chr_to_start = self.prepare_coll(query_chr_list, ref_chr_list, query_chr_to_start, ref_chr_to_start)
        if self.actual_len:
            data, gn_to_pos, rf_blk_chr, qry_blk_chr, direction_list = base.read_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)
        else:
            ratio_pair = self.ratio_pair[loop]
            data, gn_to_pos, rf_blk_chr, qry_blk_chr, direction_list = self.read_collinearity(prefix[1], prefix[0], collinearity,
                                                                              chr_list, chr_to_start, ratio_pair)
        intra = []
        intra_chr_list = []
        intra_direction = []
        i = 0
        dt = datetime.datetime.now()
        time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
        with alive_bar(len(data), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
            for block in data:
                if qry_blk_chr[i] == rf_blk_chr[i]:
                    intra.append(block)
                    intra_chr_list.append(qry_blk_chr[i])
                    intra_direction.append(direction_list[i])
                    i += 1
                else:
                    color = chr_color_dict[rf_blk_chr[i]]['chr']
                    pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, color = self.get_block_corr_color(color, direction_list[i], block, gn_to_pos, chr_plus_gap_length)

                    # 20240818 only_for_chr_margin_collinearity_plot
                    query_chr = qry_blk_chr[i]
                    ref_chr = rf_blk_chr[i]
                    i = self.plot_block(cap_judge_query_chr, cap_judge_ref_chr, query_chr, ref_chr, pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height, color, i)
                bar()
        if intra:
            i = 0
            dt = datetime.datetime.now()
            time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
            with alive_bar(len(intra), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
                for block in intra:
                    color = chr_color_dict[intra_chr_list[i]]['chr']
                    pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, color = self.get_block_corr_color(color, intra_direction[i], block, gn_to_pos, chr_plus_gap_length)

                    # 20240818 only_for_chr_margin_collinearity_plot
                    query_chr = ref_chr = intra_chr_list[i]
                    i = self.plot_block(cap_judge_query_chr, cap_judge_ref_chr, query_chr, ref_chr, pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height, color, i)
                    bar()
        logger.info(f"Plot {collinearity} finished!")

    def line_init(self):
        print()
        for key, value in vars(self).items():
            if key != "conf" and key not in ["ref_height", "height_gap", "query_height"]:
                print(key, "=", value)
        print()

        if not self.length_file:
            logger.error("Please specify your chromosome length file(Separator: ',')")
            sys.exit(1)
        if not self.input_file:
            logger.error("Please specify your input collinearity file")
            sys.exit(1)
        if not self.output_file_name:
            logger.error("Please specify your output file name")
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

    @staticmethod
    def read_len(file, prefix):
        length = pd.read_csv(file, sep='\t', header=0, index_col=None)
        length['chr'] = length['chr'].astype(str)
        length['chr'] = prefix + length['chr']
        length['length'] = length['length'].astype(int)
        return length

    def split_length(self, strip_prefix):
        length_file_list = self.length_file.split(",")
        strip_length_file_list_df = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            base.file_empty(le)
            length = self.read_len(le, strip_prefix[i])
            strip_length_file_list_df.append(length)
            i += 1
        return strip_length_file_list_df

    def split_length_adjust(self, strip_prefix, ratio_list):
        length_file_list = self.length_file.split(",")
        rm_blank = []
        for idx, le in enumerate(length_file_list):
            le = le.strip()
            if len(le) == 0:
                continue
            rm_blank.append(le)
        strip_length_file_list_df = []

        for idx, le in enumerate(rm_blank):
            base.file_empty(le)
            ratio = ratio_list[idx]
            length = self.read_len(le, strip_prefix[idx])
            length['length'] = length['length'] * ratio
            strip_length_file_list_df.append(length)
        return strip_length_file_list_df

    def get_sp_length(self):
        length_file_list = self.length_file.split(",")
        sp_length = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            base.file_empty(le)
            length = pd.read_csv(le, sep='\t', header=0, index_col=None)
            length['chr'] = length['chr'].astype(str)
            length['length'] = length['length'].astype(int)
            sp_length.append(length['length'].sum())
            i += 1
        return sp_length

    @staticmethod
    def get_expand_ratio(length_list):
        max_length = max(length_list)
        ratio_list = []
        for i in length_list:
            if i < max_length * 0.5:
                ratio_list.append(max_length * 0.5 / i)
            elif max_length * 0.5 <= i < max_length * 0.65:
                ratio_list.append(max_length * 0.65 / i)
            elif max_length * 0.65 <= i < max_length * 0.8:
                ratio_list.append(max_length * 0.8 / i)
            else:
                ratio_list.append(1)
        return ratio_list

    @staticmethod
    def get_max_chr_plus_gap_length(strip_length_file_list_df):
        length_list = []
        chr_number_max = 10
        for df in strip_length_file_list_df:
            length_list.append(df['length'].sum())
        total_chr_length = max(length_list)
        for df in strip_length_file_list_df:
            if df['length'].sum() == total_chr_length:
                chr_number_max = len(df['chr'])
        chr_plus_gap_length = total_chr_length + total_chr_length * GAP_RATIO
        return chr_plus_gap_length, chr_number_max

    def run(self):
        logger.info("Line module init and the following parameters are config information.")
        sp_number = len(self.split_file(self.input_file, ",")) + 1
        self.determine_fig_par(sp_number)
        self.line_init()
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)

        strip_chr_abbr = base.split_conf(self.remove_chromosome_prefix, ",")
        strip_collinearity_file_list = self.split_file(self.input_file, ",")
        strip_prefix = base.split_conf(self.species_name, ",")
        sp_chr_color_dict = self.get_sp_chr_color_dict(strip_prefix)

        if self.actual_len:
            strip_length_file_list_df = self.split_length(strip_prefix)
        else:
            length_list = self.get_sp_length()
            ratio_list = self.get_expand_ratio(length_list)
            strip_length_file_list_df = self.split_length_adjust(strip_prefix, ratio_list)
            self.ratio_pair = []
            for i in range(len(strip_collinearity_file_list)):
                self.ratio_pair.append([ratio_list[i], ratio_list[i + 1]])

        pair_number = len(strip_collinearity_file_list)
        new_prefix = []
        for i in range(pair_number):
            new_prefix.append([strip_prefix[i], strip_prefix[i + 1]])

        new_length_file_list_df = []
        for i in range(pair_number):
            new_length_file_list_df.append([strip_length_file_list_df[i], strip_length_file_list_df[i+1]])

        chr_plus_gap_length, mx_chr_number = self.get_max_chr_plus_gap_length(strip_length_file_list_df)
        max_gm_gap_length = chr_plus_gap_length / (1+GAP_RATIO) * GAP_RATIO / mx_chr_number
        zipped_three_pair = list(zip(strip_collinearity_file_list, new_prefix, new_length_file_list_df))

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        if self.figsize:
            fig_lt = self.figsize.split(",")
            fig, ax = plt.subplots(figsize=(float(fig_lt[0]), float(fig_lt[1])), facecolor='white')
        else:
            fig, ax = plt.subplots(figsize=(14, 14), facecolor='white')
        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1

        i = 0
        for col, prefix, length in zipped_three_pair:
            self.sub_run(col, prefix, length[0], length[1], chr_plus_gap_length, query_height, ref_height, i, last_loop, strip_chr_abbr, max_gm_gap_length, sp_chr_color_dict)
            query_height = query_height + self.height_gap
            ref_height = ref_height + self.height_gap
            i += 1

        plt.axis('off')
        # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        # ax.margins(x=0, y=0)
        plt.savefig(self.output_file_name, dpi=DPI, bbox_inches='tight')
        logger.info(f"Generate {self.output_file_name} finished!")
        sys.exit(0)

