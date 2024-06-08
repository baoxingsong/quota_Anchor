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
        self.query_height = 0.6
        self.width = 0.02
        self.dpi = 1000
        # self.block_gep = 200
        self.collinearity = config_pra['line']['collinearity']
        self.ref_length = config_pra['line']['ref_length']
        self.query_length = config_pra['line']['query_length']
        self.ref_prefix = config_pra['line']['ref_prefix']
        self.qry_prefix = config_pra['line']['query_prefix']
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

    def plot_colliearity_region(self, pos1, pos2, pos3, pos4):
        ratio = 0.2
        ref_mar = self.ref_height + self.width * 1.1
        query_mar = self.query_height - self.width * 0.1
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

    def run(self):
        ref_length = pd.read_csv(self.ref_length, sep='\t', header=0)
        ref_length['chr'] = ref_length['chr'].astype(str)
        ref_length['chr'] = self.ref_prefix + ref_length['chr']
        ref_length['length'] = ref_length['length'].astype(int)

        query_length = pd.read_csv(self.query_length, sep='\t', header=0)
        query_length['chr'] = query_length['chr'].astype(str)
        query_length['chr'] = self.qry_prefix + query_length['chr']
        query_length['length'] = query_length['length'].astype(int)
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
        plt.rcParams['font.family'] = "Times New Roman"
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
        ax.set_aspect('equal')
        total_chr_length = max(query_length['length'].sum(), ref_length['length'].sum())
        total_length = total_chr_length + total_chr_length / self.gap_ratio

        ref_gap_length = (total_length - ref_length['length'].sum()) / (len(ref_chr_list) + 1)
        ref_start_list, ref_end_list, ref_chr_to_start = self.line_get_pos_list(ref_length, ref_gap_length)
        # relative length 1
        label_x = -0.1
        label_y = self.ref_height + self.width / 2
        plt.text(label_x, label_y, self.ref_prefix, ha="center", va="center", fontsize=self.font_size, color='white')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / total_length
            ref_end_x = ref_end_list[i] / total_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, self.width, self.ref_height)
            plt.fill(x, y, facecolor='white', alpha=.7)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = self.ref_height + self.width / 2
            plt.text(label_x, label_y, ref_chr_list[i][len(self.ref_prefix):], ha="center", va="center", fontsize=self.font_size, color='black')

        query_gap_length = (total_length - query_length['length'].sum()) / (len(query_chr_list) + 1)
        query_start_list, query_end_list, query_chr_to_start = self.line_get_pos_list(query_length, query_gap_length)
        label_x = -0.1
        label_y = self.query_height + self.width / 2
        plt.text(label_x, label_y, self.qry_prefix, ha="center", va="center", fontsize=self.font_size, color='white')
        for i in range(len(query_chr_list)):
            query_start_x = query_start_list[i] / total_length
            query_end_x = query_end_list[i] / total_length
            x, y = self.plot_line_chr(query_start_x, query_end_x, self.width, self.query_height)
            plt.fill(x, y, facecolor='white', alpha=.7)
            label_x = (query_start_x + query_end_x) / 2
            label_y = self.query_height + self.width / 2
            plt.text(label_x, label_y, query_chr_list[i][len(self.qry_prefix):], ha="center", va="center", fontsize=self.font_size, color='black')
        chr_list = list(OrderedDict.fromkeys(query_chr_list + ref_chr_list))
        chr_to_start = {}
        chr_to_start.update(query_chr_to_start)
        chr_to_start.update(ref_chr_to_start)

        data, gn_to_pos, rf_blk_chr, qry_blk_chr = circle.Circle.read_collinearity(self.qry_prefix, self.ref_prefix, self.collinearity, chr_list, chr_to_start)
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

                x, y = self.plot_colliearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x)
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

            x, y = self.plot_colliearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x)
            plt.fill(x, y, facecolor=color, alpha=0.5)
            i += 1
        plt.axis('off')
        # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.savefig, dpi=int(self.dpi), bbox_inches='tight')
        sys.exit(0)
