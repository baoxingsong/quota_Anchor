import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.transforms import Bbox
import pandas as pd
from collections import OrderedDict
import numpy as np
from math import pi
from . import base
import sys
import seaborn as sns
# import time


class Line:
    def __init__(self, config_pra, parameter):
        # fig setting refers to famCircle(https://github.com/lkiko/famCircle).
        # gaps between chromosome, chr:gap = 4: 1
        self.overwrite = False
        self.remove_chromosome_prefix = "chr,CHR,Chr"
        self.chr_font_size = 7
        self.species_name_font_size = 7
        for i in config_pra.sections():
            if i == 'line':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)
        print()
        self.gap_ratio = 10
        self.ref_height = 0.3
        self.height_gap = 0.3
        self.query_height = 0.6
        self.width = 0.015
        self.dpi = 1000
        # self.block_gep = 200

    @staticmethod
    def circle(radius, x):
        if x > 0.8:
            print(x)
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
    def plot_line_chr(start, end, width, height):
        new_start = start + width / 2
        new_end = end - width / 2

        # right
        t = np.arange(-pi/2, pi/2, pi / 180)
        x = list(width / 2 * np.cos(t) + new_end)
        y = list(width / 2 * np.sin(t) + height + width / 2)

        # top
        x.append(new_end)
        y.append(height + width)
        x.append(new_start)
        y.append(height + width)

        # left
        t = np.arange(pi/2, 3*pi/2, pi / 180)
        x1 = list(width / 2 * np.cos(t) + new_start)
        y1 = list(width / 2 * np.sin(t) + height + width / 2)
        x += x1
        y += y1

        # bottom
        x.append(new_start)
        y.append(height)
        x.append(new_end)
        y.append(height)

        return x, y

    def plot_collinearity_region(self, pos1, pos2, pos3, pos4, query_height, ref_height, jg_qr_st, jg_qr_ed, jg_rf_st, jg_rf_ed):
        ratio = 0.316
        # ref_mar = ref_height + self.width * 1.1
        # query_mar = query_height - self.width * 0.1
        ref_mar = ref_height + self.width * 1.0
        query_mar = query_height
        x, y = [], []
        # p1->p3(ref_block)
        t = np.linspace(pos1, pos3, 1000)
        for i in t:
            if jg_rf_st - self.width / 2 <= i < jg_rf_st:
                will_x = i
                will_y = ref_height + self.width / 2 + self.circle(self.width / 2, jg_rf_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_rf_st <= i <= jg_rf_ed:
                will_x = i
                will_y = ref_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = ref_height + self.width / 2 + self.circle(self.width / 2, i - jg_rf_ed)
                x.append(will_x)
                y.append(will_y)
        # x.append(pos1)
        # y.append(ref_mar)
        # x.append(pos3)
        # y.append(ref_mar)

        # p3->p4
        pos3_x = x[-1]
        pos3_y = y[-1]
        if jg_qr_st - self.width / 2 < pos4 < jg_qr_st:
            will_x = pos4
            will_y = query_height + self.width / 2 - self.circle(self.width / 2, jg_qr_st - pos4)
        elif jg_qr_st <= pos4 <= jg_qr_ed:
            will_x = pos4
            will_y = query_mar
        else:
            will_x = pos4
            will_y = query_height + self.width / 2 - self.circle(self.width / 2, pos4 - jg_qr_ed)
        t = np.arange(0, 1.01, 0.005)
        dx = will_x - pos3_x
        dy = will_y - pos3_y
        p1, p2, p3, p4 = pos3_x, dx * ratio + pos3_x, -dx * ratio + will_x, will_x
        x1 = base.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos3_y, (1 - ratio) * dy + pos3_y, -(1-ratio) * dy + will_y, will_y
        y1 = base.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        # p4->p2
        t = np.linspace(pos4, pos2, 1000)
        for i in t:
            if jg_qr_st - self.width / 2 < i < jg_qr_st:
                will_x = i
                will_y = query_height + self.width / 2 - self.circle(self.width / 2, jg_qr_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_qr_st <= i <= jg_qr_ed:
                will_x = i
                will_y = query_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = query_height + self.width / 2 - self.circle(self.width / 2, i - jg_qr_ed)
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
        p1, p2, p3, p4 = pos2_x, dx * ratio + pos2_x, -dx * ratio + pos1_x, pos1_x
        x1 = base.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos2_y, (1 - ratio) * dy + pos2_y, -(1-ratio) * dy + pos1_y, pos1_y
        y1 = base.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        return x, y

    def sub_run(self, collinearity, prefix, ref_length, query_length, total_length, query_height, ref_height, loop, last_loop, strip_chr_abbr):
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
        ref_gap_length = (total_length - ref_length['length'].sum()) / (len(ref_chr_list) + 1)
        ref_start_list, ref_end_list, ref_chr_to_start = self.line_get_pos_list(ref_length, ref_gap_length)
        # 20240818 only_for_chr_margin_collinearity_plot
        only_for_chr_margin_collinearity_plot_ref_zip = zip(ref_chr_list, ref_start_list, ref_end_list)
        only_for_chr_margin_collinearity_plot_ref_chr_to_position = {}
        for ch, start, end in only_for_chr_margin_collinearity_plot_ref_zip:
            start_x = start / total_length + self.width / 2
            end_x = end / total_length - self.width / 2
            only_for_chr_margin_collinearity_plot_ref_chr_to_position[ch] = [start_x, end_x]
        # relative length 1
        label_x = -0.1
        label_y = ref_height + self.width / 2
        plt.text(label_x, label_y, prefix[0], ha="center", va="center", fontsize=self.species_name_font_size, color='black')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / total_length
            ref_end_x = ref_end_list[i] / total_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, self.width, ref_height)
            plt.fill(x, y, facecolor='white', alpha=.7, edgecolor='black', linewidth=1)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = ref_height + self.width / 2
            text_chr = ref_chr_list[i][len(prefix[0]):]
            for abbr in strip_chr_abbr:
                if text_chr.startswith(abbr):
                    text_chr = text_chr[len(abbr):]
            plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size, color='black')

        query_gap_length = (total_length - query_length['length'].sum()) / (len(query_chr_list) + 1)
        query_start_list, query_end_list, query_chr_to_start = self.line_get_pos_list(query_length, query_gap_length)
        # 20240818 only_for_chr_margin_collinearity_plot
        only_for_chr_margin_collinearity_plot_query_zip = zip(query_chr_list, query_start_list, query_end_list)
        only_for_chr_margin_collinearity_plot_query_chr_to_position = {}
        for ch, start, end in only_for_chr_margin_collinearity_plot_query_zip:
            start_x = start / total_length + self.width / 2
            end_x = end / total_length - self.width / 2
            only_for_chr_margin_collinearity_plot_query_chr_to_position[ch] = [start_x, end_x]
        if loop == last_loop:
            label_x = -0.1
            label_y = query_height + self.width / 2
            plt.text(label_x, label_y, prefix[1], ha="center", va="center", fontsize=self.species_name_font_size, color='black')
            for i in range(len(query_chr_list)):
                query_start_x = query_start_list[i] / total_length
                query_end_x = query_end_list[i] / total_length
                x, y = self.plot_line_chr(query_start_x, query_end_x, self.width, query_height)
                plt.fill(x, y, facecolor='white', alpha=0.7, edgecolor='black', linewidth=1)
                label_x = (query_start_x + query_end_x) / 2
                label_y = query_height + self.width / 2
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
        for block in data:
            if qry_blk_chr[i] == rf_blk_chr[i]:
                intra.append(block)
                intra_chr_list.append(qry_blk_chr[i])
                i += 1
            else:
                color = chr_color_dict[rf_blk_chr[i]]['chr']
                id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
                pos1_coord_x = pos1 / total_length
                pos2_coord_x = pos2 / total_length
                pos3_coord_x = pos3 / total_length
                pos4_coord_x = pos4 / total_length
                # 20240818 only_for_chr_margin_collinearity_plot
                query_chr = qry_blk_chr[i]
                ref_chr = rf_blk_chr[i]
                judge_fake_query_start_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][0]
                judge_fake_query_end_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][1]
                judge_fake_ref_start_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][0]
                judge_fake_ref_end_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][1]
                x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                     judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
                plt.fill(x, y, facecolor=color, alpha=0.7)
                i += 1
        i = 0
        for block in intra:
            # color = cmap(round(intra_color_index[i] / len(data), 2))
            color = chr_color_dict[intra_chr_list[i]]['chr']
            id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
            pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
            pos1_coord_x = pos1 / total_length
            pos2_coord_x = pos2 / total_length
            pos3_coord_x = pos3 / total_length
            pos4_coord_x = pos4 / total_length

            # 20240818 only_for_chr_margin_collinearity_plot
            query_chr = ref_chr = intra_chr_list[i]
            judge_fake_query_start_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][0]
            judge_fake_query_end_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][1]
            judge_fake_ref_start_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][0]
            judge_fake_ref_end_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][1]
            x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                 judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
            plt.fill(x, y, facecolor=color, alpha=0.7)
            i += 1

    def run(self):
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)
        chr_abbr = self.remove_chromosome_prefix.split(',')
        strip_chr_abbr = []
        for i in chr_abbr:
            if len(i) == 0:
                continue
            i = i.strip()
            strip_chr_abbr.append(i)

        # split collinearity in order to loop
        collinearity_file_list = self.input_file.split(",")
        strip_collinearity_file_list = []
        for col in collinearity_file_list:
            if len(col) == 0:
                continue
            col = col.strip()
            base.file_empty(col)
            strip_collinearity_file_list.append(col)

        # split prefix in order to pair
        prefix = self.species_name.split(",")
        strip_prefix = []
        new_prefix = []
        for prx in prefix:
            if len(prx) == 0:
                continue
            prx = prx.strip()
            strip_prefix.append(prx)
        for i in range(len(strip_collinearity_file_list)):
            new_prefix.append([strip_prefix[i], strip_prefix[i + 1]])

        # split length get dataframe
        length_file_list = self.length_file.split(",")
        strip_length_file_list_df = []
        new_length_file_list_df = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            base.file_empty(le)
            length = pd.read_csv(le, sep='\t', header=0)
            # print(length)
            length['chr'] = length['chr'].astype(str)
            length['chr'] = strip_prefix[i] + length['chr']
            length['length'] = length['length'].astype(int)
            strip_length_file_list_df.append(length)
            i += 1
        # print(strip_length_file_list_df)
        for i in range(len(strip_collinearity_file_list)):
            new_length_file_list_df.append([strip_length_file_list_df[i], strip_length_file_list_df[i+1]])

        # get maximum chr total length
        length_list = []
        for df in strip_length_file_list_df:
            length_list.append(df['length'].sum())
        total_chr_length = max(length_list)
        total_length = total_chr_length + total_chr_length / self.gap_ratio

        zipped_three_pair = list(zip(strip_collinearity_file_list, new_prefix, new_length_file_list_df))
        # print(zipped_three_pair)
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='white')
        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1
        # print(zipped_three_pair)
        i = 0
        for col, prefix, length in zipped_three_pair:
            self.sub_run(col, prefix, length[0], length[1], total_length, query_height, ref_height, i, last_loop, strip_chr_abbr)
            query_height = query_height + self.height_gap

            ref_height = ref_height + self.height_gap
            # sys.exit(1)
            i += 1

        plt.axis('off')
        # # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.output_file_name, dpi=int(self.dpi), bbox_inches='tight')
        print("produce", self.output_file_name, "success")
        sys.exit(0)
