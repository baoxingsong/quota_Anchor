import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.transforms import Bbox
import pandas as pd
from collections import OrderedDict
import numpy as np
from math import pi
from . import base
import sys
# import time


class Line:
    def __init__(self, config_pra):
        # fig setting refers to famCircle(https://github.com/lkiko/famCircle).
        # gaps between chromosome, chr:gap = 4: 1
        self.gap_ratio = 6
        self.ref_height = 0.3
        self.height_gap = 0.08
        self.query_height = 0.38
        self.width = 0.003
        self.dpi = 300
        # self.block_gep = 200
        self.collinearity = config_pra['line']['collinearity']
        self.length_file = config_pra['line']['length_file']
        self.prefix = config_pra['line']['prefix']
        self.remove_chromosome_prefix = config_pra['line']['remove_chromosome_prefix']
        self.font_size = int(config_pra['line']['text_font_size'])
        self.savefig = config_pra['line']['savefig']

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
    def read_proali_collinearity(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start):
        # left -> ref;  right -> query
        data = []
        ref_chr_list = []
        query_chr_list = []
        ref_block_direction = []
        # var = ""
        block_index = 0
        block = []

        with open(collinearity) as f:
            print("read", collinearity, "....")
            _ = next(f)
            _ = next(f)
            flag = True
            flag_number = True
            for line in f:
                if line.startswith("#block begin"):
                    # var = line.split()[2]
                    flag = True
                    flag_number = True
                    if block:
                        data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                        # ref_block_direction.append(var)
                        block = []
                        # var = ""
                elif line.startswith("#block end"):
                    continue
                else:
                    if flag:
                        chr_pair = [line.split()[0], line.split()[3]]
                        if ref_prefix + chr_pair[0] not in chr_list or qry_prefix + chr_pair[1] not in chr_list:
                            # print("block not in length.txt")
                            flag = False
                        else:
                            if flag_number:
                                flag = True
                                ref_chr = ref_prefix + chr_pair[0]
                                ref_chr_list.append(ref_chr)
                                query_chr = qry_prefix + chr_pair[1]
                                query_chr_list.append(query_chr)
                                ref_block_direction.append(line.split()[6])
                                block_index += 1
                                flag_number = False
                        if flag:
                            line_list = line.split()
                            block.append([chr_to_start[ref_chr] + int(line_list[2]), chr_to_start[query_chr] + int(line_list[5])])
                    else:
                        continue
            if block:
                data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                # ref_block_direction.append(var)
            print("parse", collinearity, "success")
        return data, ref_chr_list, query_chr_list, ref_block_direction

    @staticmethod
    def set_husl_palette():
        colors = sns.husl_palette(as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("husl_256", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def set_palette():
        colors = ['#FFA07A', '#4b6870', '#4169E1', '#9370DB']
        color_inversion = ["#66AD56", "#66AD56"]
        color_normal_grey = ["#F0F0F0", "#F0F0F0"]
        cmap_normal = LinearSegmentedColormap.from_list('custom_cmap', color_normal_grey)
        cmap_inversion = LinearSegmentedColormap.from_list('custom_cmap', color_inversion)
        return cmap_normal, cmap_inversion

    @staticmethod
    def get_margin_radius(max_total_length, gap_radio, max_number, margin_radio, width, total_length):
        fake_gap_total_length = max_total_length / gap_radio
        every_margin = fake_gap_total_length / (max_number + 1)
        x = every_margin * margin_radio / total_length
        radius_to_sub_run = (x ** 2 + 1 / 4 * width ** 2) / (2 * x)
        change_length = (radius_to_sub_run - x)
        return radius_to_sub_run, change_length
    
    @staticmethod
    def get_change_radian(margin_radius ,change_length):
        half_radian = np.arccos(change_length / margin_radius)
        radian = 2 * half_radian
        return radian

    @staticmethod
    def plot_line_chr(start, end, width, height, radian, change_length, margin_radius):
        new_start = start + width / 2
        new_end = end - width / 2
        # new_start = start
        # new_end = end

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
        x1 = list(width/2 * np.cos(t) + new_start)
        y1 = list(width/2 * np.sin(t) + height + width / 2)
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
            
            # x.append(i)
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
        # will_x = pos4
        # will_y = query_mar
        
        t = np.arange(0, 1.01, 0.005)
        dx = will_x - pos3_x
        dy = will_y - pos3_y
        p1, p2, p3, p4 = pos3_x, dx * ratio + pos3_x, -dx * ratio + will_x, will_x
        x1 = base.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos3_y, (1 - ratio) * dy + pos3_y, -(1 - ratio) * dy + will_y, will_y
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
            # x.append(i)
            # y.append(query_mar)

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
        p1, p2, p3, p4 = pos2_y, (1 - ratio) * dy + pos2_y, -(1 - ratio) * dy + pos1_y, pos1_y
        y1 = base.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        return x, y

    def sub_run(self, collinearity, prefix, ref_length, query_length, total_length, query_height, ref_height, loop, last_loop, strip_chr_abbr, radian, change_length, margin_radius):
        ref_chr_list = ref_length['chr'].tolist()
        query_chr_list = query_length['chr'].tolist()
        chr_color_dict_normal = {}
        chr_color_dict_inversion = {}
        chr_color_dict_trans = {}
        chr_color_dict_trans_inversion = {}
        cmap_normal, cmap_inversion = self.set_palette()
        # cmap = self.set_husl_palette()
        i = 1
        for ch in ref_chr_list:
            color_dict_normal = {'chr': cmap_normal(round(i / len(ref_chr_list), 2))}
            chr_color_dict_normal[ch] = color_dict_normal
            
            color_dict_inversion = {'chr': cmap_inversion(round(i / len(ref_chr_list), 2))}
            chr_color_dict_inversion[ch] = color_dict_inversion
            
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
        # plt.text(label_x, label_y, prefix[0], ha="center", va="center", fontsize=self.font_size, color='black')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / total_length
            ref_end_x = ref_end_list[i] / total_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, self.width, ref_height, radian, change_length, margin_radius)
            # modify outline?
            if loop in [0, 1, 2, 3 ,4]:
                plt.fill(x, y, facecolor='#E83828', alpha=0.5)
            elif loop in [5,6]:
                plt.fill(x, y, facecolor='#E78A00', alpha=0.5)
            else:
                plt.fill(x, y, facecolor='#3F40C3', alpha=0.5)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = ref_height + self.width / 2
            text_chr = ref_chr_list[i][len(prefix[0]):]
            for abbr in strip_chr_abbr:
                if text_chr.startswith(abbr):
                    text_chr = text_chr[len(abbr):]
            # plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.font_size, color='black')

        query_gap_length = (total_length - query_length['length'].sum()) / (len(query_chr_list) + 1)
        query_start_list, query_end_list, query_chr_to_start = self.line_get_pos_list(query_length, query_gap_length)
        # 20240818 only_for_chr_margin_collinearity_plot
        only_for_chr_margin_collinearity_plot_query_zip = zip(query_chr_list, query_start_list, query_end_list)
        only_for_chr_margin_collinearity_plot_query_chr_to_position = {}
        for ch, start, end in only_for_chr_margin_collinearity_plot_query_zip:
            start_x = start / total_length + self.width / 2
            end_x = end / total_length - self.width / 2
            only_for_chr_margin_collinearity_plot_query_chr_to_position[ch] = [start_x, end_x]
        # print(last_loop)
        # print(loop)
        if loop == last_loop:
            label_x = -0.1
            label_y = query_height + self.width / 2
            # print(prefix[1])
            # print(label_x)
            # print(label_y)
            # plt.text(label_x, label_y, prefix[1], ha="center", va="center", fontsize=self.font_size, color='black')
            for i in range(len(query_chr_list)):
                query_start_x = query_start_list[i] / total_length
                query_end_x = query_end_list[i] / total_length
                x, y = self.plot_line_chr(query_start_x, query_end_x, self.width, query_height, radian, change_length, margin_radius)
                # modify outline ?
                plt.fill(x, y, facecolor='#3F40C3', alpha=0.5)
                label_x = (query_start_x + query_end_x) / 2
                label_y = query_height + self.width / 2
                text_chr = query_chr_list[i][len(prefix[1]):]
                for abbr in strip_chr_abbr:
                    if text_chr.startswith(abbr):
                        text_chr = text_chr[len(abbr):]
                # plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.font_size, color='black')
        # print(ref_chr_list)
        # print(query_chr_list)
        chr_list = list(OrderedDict.fromkeys(query_chr_list + ref_chr_list))
        chr_to_start = {}
        chr_to_start.update(query_chr_to_start)
        chr_to_start.update(ref_chr_to_start)
        # print(prefix)
        # print(collinearity)
        # print(chr_list)
        # print(chr_to_start)
        # data, gn_to_pos, rf_blk_chr, qry_blk_chr = base.read_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)

        data, rf_blk_chr, qry_blk_chr, ref_block_direction = self.read_proali_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)
        intra = []
        intra_chr_list = []
        i = 0
        # intra_color_index = []
        for block in data:
            if qry_blk_chr[i] == rf_blk_chr[i]:
                # print("aaa")
                intra.append(block)
                intra_chr_list.append(qry_blk_chr[i])
                # intra_color_index.append(i)
                i += 1
            else:
                # color = cmap(round(i / len(data), 2))
                if ref_block_direction[i] == '+':
                    color = chr_color_dict_normal[rf_blk_chr[i]]['chr']
                if ref_block_direction[i] == '-':
                    color = chr_color_dict_inversion[rf_blk_chr[i]]['chr']
                # if ref_block_direction[i] == 'Trans':
                #     color = chr_color_dict_trans[rf_blk_chr[i]]['chr']
                # if ref_block_direction[i] == 'Inversion/Trans':
                #     color = chr_color_dict_trans_inversion[rf_blk_chr[i]]['chr']
                # if ref_block_direction[i] == 'Trans/Inversion':
                #     color = chr_color_dict_trans_inversion[rf_blk_chr[i]]['chr']
                # color = chr_color_dict[rf_blk_chr[i]]['chr']
                # id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                # pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
                pos1, pos2, pos3, pos4 = block[0], block[1], block[2], block[3]
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
                plt.fill(x, y, facecolor=color, alpha=0.8)
                i += 1
        i = 0
        for block in intra:
            # color = cmap(round(intra_color_index[i] / len(data), 2))
            color = chr_color_dict_normal[intra_chr_list[i]]['chr']
            # if ref_block_direction[i] == 'Normal':
            #     color = chr_color_dict_normal[rf_blk_chr[i]]['chr']
            # if ref_block_direction[i] == 'Inversion':
            #     color = chr_color_dict_inversion[rf_blk_chr[i]]['chr']
            # if ref_block_direction[i] == 'Trans':
            #     color = chr_color_dict_trans[rf_blk_chr[i]]['chr']
            # if ref_block_direction[i] == 'Inversion/Trans':
            #     color = chr_color_dict_trans_inversion[rf_blk_chr[i]]['chr']
            # id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
            # pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
            pos1, pos2, pos3, pos4 = block[0], block[1], block[2], block[3]
            pos1_coord_x = pos1 / total_length
            pos2_coord_x = pos2 / total_length
            pos3_coord_x = pos3 / total_length
            pos4_coord_x = pos4 / total_length

            query_chr = ref_chr = intra_chr_list[i]
            judge_fake_query_start_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][0]
            judge_fake_query_end_x = only_for_chr_margin_collinearity_plot_query_chr_to_position[query_chr][1]
            judge_fake_ref_start_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][0]
            judge_fake_ref_end_x = only_for_chr_margin_collinearity_plot_ref_chr_to_position[ref_chr][1]
            x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                 judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
            plt.fill(x, y, facecolor=color, alpha=0)
            i += 1

    def run(self):
        chr_abbr = self.remove_chromosome_prefix.split(',')
        strip_chr_abbr = []
        for i in chr_abbr:
            if len(i) == 0:
                continue
            i = i.strip()
            strip_chr_abbr.append(i)

        # split collinearity in order to loop
        collinearity_file_list = self.collinearity.split(",")
        strip_collinearity_file_list = []
        for col in collinearity_file_list:
            if len(col) == 0:
                continue
            col = col.strip()
            strip_collinearity_file_list.append(col)

        # split prefix in order to pair
        prefix = self.prefix.split(",")
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
        # 20241024 to get chr margin avoid overlap
        chr_number  = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            length = pd.read_csv(le, sep='\t', header=0)
            # print(length)
            length['chr'] = length['chr'].astype(str)
            length['chr'] = strip_prefix[i] + length['chr']
            length['length'] = length['length'].astype(int)
            strip_length_file_list_df.append(length)
            chr_number.append(len(length['chr'].tolist()))
            i += 1
        # 20241024 to get chr margin avoid overlap
        max_number = max(chr_number)
        # print(strip_length_file_list_df)
        for i in range(len(strip_collinearity_file_list)):
            new_length_file_list_df.append([strip_length_file_list_df[i], strip_length_file_list_df[i+1]])

        # get maximum chr total length
        length_list = []
        for df in strip_length_file_list_df:
            length_list.append(df['length'].sum())
        total_chr_length = max(length_list)
        total_length = total_chr_length + total_chr_length / self.gap_ratio
        
        # 20241024 to get chr margin avoid overlap
        radito_margin_gap = 0.25
        margin_radius, change_length = self.get_margin_radius(total_chr_length, self.gap_ratio, max_number, radito_margin_gap, self.width, total_length)
        radian = self.get_change_radian(margin_radius, change_length)

        zipped_three_pair = list(zip(strip_collinearity_file_list, new_prefix, new_length_file_list_df))
        # print(zipped_three_pair)
        
        plt.rcParams['font.family'] = "Times New Roman"
        fig, ax = plt.subplots(figsize=(14, 8), facecolor='white')
        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1
        # print(zipped_three_pair)
        i = 0
        for col, prefix, length in zipped_three_pair:
            # print(length)
            # print(last_loop)
            # print(i)
            # print(ref_height)
            # print(query_height)

            self.sub_run(col, prefix, length[0], length[1], total_length, query_height, ref_height, i, last_loop, strip_chr_abbr, radian, change_length, margin_radius)
            query_height = query_height + self.height_gap

            ref_height = ref_height + self.height_gap
            # print(ref_height)
            # print(query_height)
            # sys.exit(1)
            i += 1

        plt.axis('off')
        # # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.savefig, dpi=int(self.dpi))
        print("produce", self.savefig, "success")
        sys.exit(0)
