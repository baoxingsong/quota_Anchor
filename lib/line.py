import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.transforms import Bbox
import pandas as pd
from collections import OrderedDict
import numpy as np
from math import pi
from . import circle
import sys
# import time


class Line:
    def __init__(self, config_pra):
        # fig setting refers to famCircle(https://github.com/lkiko/famCircle).
        # gaps between chromosome, chr:gap = 4: 1
        self.gap_ratio = 6
        self.ref_height = 0.3
        self.height_gap = 0.3
        self.query_height = 0.6
        self.width = 0.015
        self.dpi = 1500
        # self.block_gep = 200
        self.collinearity = config_pra['line']['collinearity']
        self.length_file = config_pra['line']['length_file']
        self.prefix = config_pra['line']['prefix']
        self.font_size = int(config_pra['line']['text_font_size'])
        self.savefig = config_pra['line']['savefig']

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

    def plot_colliearity_region(self, pos1, pos2, pos3, pos4, query_height, ref_height):
        ratio = 0.2
        ref_mar = ref_height + self.width * 1.1
        query_mar = query_height - self.width * 0.1
        x, y = [], []
        # p1->p3(ref_block)
        x.append(pos1)
        y.append(ref_mar)
        x.append(pos3)
        y.append(ref_mar)

        # p3->p4
        t = np.arange(0, 1.01, 0.01)
        dx = pos4 - pos3
        dy = query_mar - ref_mar
        p1, p2, p3, p4 = pos3, dx * ratio + pos3, -dx * ratio + pos4, pos4
        x1 = circle.Circle.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = ref_mar, (1 - ratio) * dy + ref_mar, -(1-ratio) * dy + query_mar, query_mar
        y1 = circle.Circle.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        # p4->p2
        x.append(pos4)
        y.append(query_mar)
        x.append(pos2)
        y.append(query_mar)

        # p2->p1
        t = np.arange(0, 1.01, 0.01)
        dx = pos1 - pos2
        dy = ref_mar - query_mar
        p1, p2, p3, p4 = pos2, dx * ratio + pos2, -dx * ratio + pos1, pos1
        x1 = circle.Circle.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = query_mar, (1 - ratio) * dy + query_mar, -(1-ratio) * dy + ref_mar, ref_mar
        y1 = circle.Circle.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        return x, y

    def sub_run(self, collinearity, prefix, ref_length, query_length, total_length, query_height, ref_height, loop, last_loop):
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
        # relative length 1
        label_x = -0.1
        label_y = ref_height + self.width / 2
        plt.text(label_x, label_y, prefix[0], ha="center", va="center", fontsize=self.font_size, color='white')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / total_length
            ref_end_x = ref_end_list[i] / total_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, self.width, ref_height)
            plt.fill(x, y, facecolor='white', alpha=.7)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = ref_height + self.width / 2
            plt.text(label_x, label_y, ref_chr_list[i][len(prefix[0]):], ha="center", va="center", fontsize=self.font_size, color='black')

        query_gap_length = (total_length - query_length['length'].sum()) / (len(query_chr_list) + 1)
        query_start_list, query_end_list, query_chr_to_start = self.line_get_pos_list(query_length, query_gap_length)
        print(last_loop)
        print(loop)
        if loop == last_loop:
            label_x = -0.1
            label_y = query_height + self.width / 2
            # print(prefix[1])
            # print(label_x)
            # print(label_y)
            plt.text(label_x, label_y, prefix[1], ha="center", va="center", fontsize=self.font_size, color='white')
            for i in range(len(query_chr_list)):
                query_start_x = query_start_list[i] / total_length
                query_end_x = query_end_list[i] / total_length
                x, y = self.plot_line_chr(query_start_x, query_end_x, self.width, query_height)
                plt.fill(x, y, facecolor='white', alpha=0.7)
                label_x = (query_start_x + query_end_x) / 2
                label_y = query_height + self.width / 2
                plt.text(label_x, label_y, query_chr_list[i][len(prefix[1]):], ha="center", va="center", fontsize=self.font_size, color='black')
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
        data, gn_to_pos, rf_blk_chr, qry_blk_chr = circle.Circle.read_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)
        intra = []
        intra_chr_list = []
        i = 0
        # intra_color_index = []
        for block in data:
            if qry_blk_chr[i] == rf_blk_chr[i]:
                intra.append(block)
                intra_chr_list.append(qry_blk_chr[i])
                # intra_color_index.append(i)
                i += 1
            else:
                # color = cmap(round(i / len(data), 2))
                color = chr_color_dict[rf_blk_chr[i]]['chr']
                id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                pos1, pos2, pos3, pos4 = gn_to_pos[id1], gn_to_pos[id2], gn_to_pos[id3], gn_to_pos[id4]
                pos1_coord_x = pos1 / total_length
                pos2_coord_x = pos2 / total_length
                pos3_coord_x = pos3 / total_length
                pos4_coord_x = pos4 / total_length

                x, y = self.plot_colliearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height)
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

            x, y = self.plot_colliearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height)
            plt.fill(x, y, facecolor=color, alpha=0.7)
            i += 1

    def run(self):
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
        plt.rcParams['font.family'] = "Times New Roman"
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1
        print(zipped_three_pair)
        i = 0
        for col, prefix, length in zipped_three_pair:
            # print(length)
            # print(last_loop)
            # print(i)
            # print(ref_height)
            # print(query_height)

            self.sub_run(col, prefix, length[0], length[1], total_length, query_height, ref_height, i, last_loop)
            query_height = query_height + self.height_gap

            ref_height = ref_height + self.height_gap
            # print(ref_height)
            # print(query_height)
            # sys.exit(1)
            i += 1

        plt.axis('off')
        # # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.savefig, dpi=int(self.dpi), bbox_inches='tight')
        sys.exit(0)
