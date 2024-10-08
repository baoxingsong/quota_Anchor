import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
# from matplotlib.transforms import Bbox
import numpy as np
from math import pi
import sys
from . import base


class Circle:
    def __init__(self, config_pra, parameter):
        # fig setting refers to famCircle(https://github.com/lkiko/famCircle).
        # gaps between chromosome circle, chr:gap = 4: 1
        self.overwrite = False
        self.ref_name = "Reference_species"
        self.query_name = "Query_species"
        self.remove_chromosome_prefix = "CHR,Chr,chr"
        self.chr_font_size = 7
        self.species_name_font_size = 7
        for i in config_pra.sections():
            if i == 'circle':
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
        self.gap_ratio = 4
        self.inner_radius = 0.3
        self.outer_radius = 0.31
        self.ring_width = 0.01
        self.dpi = 1000
        # self.block_gep = 200


    @staticmethod
    def set_husl_palette():
        colors = sns.husl_palette(as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("husl_256", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def get_pos_list(df_dup, gap_length):
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                start_list.append(0)
                end_list.append(lgh)
                chr_to_start[chr_name[i]] = 0
                i += 1
            else:
                start_list.append(end_list[i-1] + gap_length)
                end_list.append(start_list[i] + lgh)
                chr_to_start[chr_name[i]] = start_list[-1]
                i += 1
        return start_list, end_list, chr_to_start

    @staticmethod
    def length_to_radian(chr_pos, total_length):
        radian_list = []
        for start, end in chr_pos:
            radian_start = 2 * pi * (start / total_length)
            radian_end = 2 * pi * (end / total_length)
            radian_list.append((radian_start, radian_end))
        return radian_list

    @staticmethod
    def pos_to_radian(pos, total_length):
        pos_radian = 2 * pi * (pos / total_length)
        return pos_radian

    def plot_legend_first(self) -> tuple[list, list]: 
        x = [self.inner_radius, self.inner_radius, self.inner_radius - self.ring_width *2, self.inner_radius - self.ring_width *2]
        y = [self.inner_radius - self.ring_width, self.inner_radius, self.inner_radius, self.inner_radius - self.ring_width]
        return x, y

    def plot_legend_second(self) -> tuple[list, list]: 
        x = [self.inner_radius, self.inner_radius, self.inner_radius - self.ring_width *2, self.inner_radius - self.ring_width *2]
        y = [self.inner_radius - 3 * self.ring_width, self.inner_radius - 2 * self.ring_width, 
             self.inner_radius - 2 * self.ring_width, self.inner_radius - 3 * self.ring_width]
        return x, y
    @staticmethod
    def plot_circle_chr(start_radian, end_radian, radius_a, radius_b, ring_width, ring_radian):
        new_start_radian = start_radian + ring_radian
        new_end_radian = end_radian - ring_radian
        t = np.arange(new_start_radian, new_end_radian, pi / 720)
        x = list(radius_a * np.cos(t))
        y = list(radius_a * np.sin(t))

        t = np.arange(new_end_radian+pi, new_end_radian, -pi / 30)
        x1 = list((radius_a + radius_b) / 2 * np.cos(new_end_radian) + ring_width / 2 * np.cos(t))
        y1 = list((radius_a + radius_b) / 2 * np.sin(new_end_radian) + ring_width / 2 * np.sin(t))
        x += x1
        y += y1

        t = np.arange(new_end_radian, new_start_radian, -pi / 720)
        x += list(radius_b * np.cos(t))
        y += list(radius_b * np.sin(t))

        t = np.arange(new_start_radian, new_start_radian-pi, -pi / 30)
        x1 = list((radius_a + radius_b) / 2 * np.cos(new_start_radian) + ring_width / 2 * np.cos(t))
        y1 = list((radius_a + radius_b) / 2 * np.sin(new_start_radian) + ring_width / 2 * np.sin(t))
        x += x1
        y += y1
        x = list(x)
        y = list(y)
        x.append(x[0])
        y.append(y[0])
        # print(x, y)
        return x, y

    @staticmethod
    def plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, radius):
        ratio = 0.382
        # pos1 -> pos3(ref region)
        t = np.arange(pos1_radian, pos3_radian, pi / 720)
        x = list(radius * np.cos(t))
        y = list(radius * np.sin(t))
        assert pos1_radian < pos3_radian

        # pos3 -> pos4
        pos3_x = radius * np.cos(pos3_radian)
        pos3_y = radius * np.sin(pos3_radian)
        pos4_x = radius * np.cos(pos4_radian)
        pos4_y = radius * np.sin(pos4_radian)
        t = np.arange(0, 1.01, 0.01)
        p0, p1, p2, p3 = pos3_x, pos3_x * ratio, pos4_x * ratio, pos4_x
        x1 = base.bezier3(p0, p1, p2, p3, t)
        p0, p1, p2, p3 = pos3_y, pos3_y * ratio, pos4_y * ratio, pos4_y
        y1 = base.bezier3(p0, p1, p2, p3, t)
        x = x + x1
        y = y + y1

        # pos4 -> pos2
        t = np.arange(min(pos4_radian, pos2_radian), max(pos4_radian, pos2_radian), pi / 720)
        x1 = list(radius * np.cos(t))
        y1 = list(radius * np.sin(t))
        if pos4_radian < pos2_radian:
            pass
        else:
            x1 = x1[::-1]
            y1 = y1[::-1]
        x = x + x1
        y = y + y1

        # pos2 -> pos1
        pos2_x = radius * np.cos(pos2_radian)
        pos2_y = radius * np.sin(pos2_radian)
        pos1_x = radius * np.cos(pos1_radian)
        pos1_y = radius * np.sin(pos1_radian)
        t = np.arange(0, 1.01, 0.01)
        p0, p1, p2, p3 = pos2_x, pos2_x * ratio, pos1_x * ratio, pos1_x
        x1 = base.bezier3(p0, p1, p2, p3, t)
        p0, p1, p2, p3 = pos2_y, pos2_y * ratio, pos1_y * ratio, pos1_y
        y1 = base.bezier3(p0, p1, p2, p3, t)
        x += x1
        y += y1
        x.append(x[0])
        y.append(y[0])
        return x, y

    @staticmethod
    def text_rotation(start, end, total_length):
        angle = (start + end) / 2 / total_length * 360
        if 0 <= angle < 180:
            return angle - 90
        if 180 <= angle <= 360:
            return angle - 270

    def run(self):
        base.file_empty(self.input_file)
        base.file_empty(self.ref_length)
        base.file_empty(self.query_length)
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)
        chr_abbr = self.remove_chromosome_prefix.split(',')
        strip_chr_abbr = []
        for i in chr_abbr:
            if len(i) == 0:
                continue
            i = i.strip()
            strip_chr_abbr.append(i)

        ref_length = pd.read_csv(self.ref_length, sep='\t', header=0)
        ref_length['chr'] = ref_length['chr'].astype(str)
        ref_length['chr'] = self.ref_name + ref_length['chr']
        ref_length['length'] = ref_length['length'].astype(int)

        query_length = pd.read_csv(self.query_length, sep='\t', header=0)
        query_length['chr'] = query_length['chr'].astype(str)
        query_length['chr'] = self.query_name + query_length['chr']
        query_length['length'] = query_length['length'].astype(int)
        df_dup = pd.concat([ref_length, query_length])
        df_dup.reset_index(inplace=True, drop=True)
        df_dup.drop_duplicates(inplace=True)

        # intra 10 inter 20 for zm sb
        chr_list = df_dup['chr'].tolist()

        # husl_palette and set line, chr color
        cmap = self.set_husl_palette()
        chr_color_dict = {}
        i = 1
        for ch in chr_list:
            color_dict = {'chr': cmap(round(i / len(chr_list), 2))}
            chr_color_dict[ch] = color_dict
            i += 1
        # figure geometry info
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='white')
        # fig.patch.set_alpha(0.95)
        ax.set_aspect('equal')
        total_chr_length = df_dup['length'].sum()
        average_length = total_chr_length / len(chr_list)
        gap_length = int(average_length / self.gap_ratio)
        total_length = total_chr_length + len(chr_list) * gap_length

        start_list, end_list, chr_to_start = self.get_pos_list(df_dup, gap_length)
        chr_pos = list(zip(start_list, end_list))
        chr_radian_pos = self.length_to_radian(chr_pos, total_length)
        i = 0
        ring_radian = np.arctan(self.ring_width / 2 / (self.inner_radius + self.ring_width / 2))
        for ch in chr_list:
            x, y = self.plot_circle_chr(chr_radian_pos[i][0], chr_radian_pos[i][1], self.inner_radius, self.outer_radius, self.ring_width, ring_radian)
            # color = chr_color_dict[ch]['chr']
            if not ch.startswith(str(self.ref_name)) and ch.startswith(str(self.query_name)):
                color = "blue"
                ch = ch[len(self.query_name):]
            else:
                color = "red"
                ch = ch[len(self.ref_name):]
            plt.fill(x, y, facecolor=color, alpha=.5)
            # plt.plot(x, y, color='black')
            label_x = self.outer_radius * 1.04 * np.cos((chr_radian_pos[i][0] + chr_radian_pos[i][1]) / 2)
            label_y = self.outer_radius * 1.04 * np.sin((chr_radian_pos[i][0] + chr_radian_pos[i][1]) / 2)
            angle = self.text_rotation(chr_pos[i][0], chr_pos[i][1], total_length)
            for abbr in strip_chr_abbr:
                if ch.startswith(abbr):
                    ch = ch[len(abbr):]
            plt.text(label_x, label_y, ch, ha="center", va="center", fontsize=self.chr_font_size, color='black', rotation=angle)
            i += 1
        data, gene_pos_dict, ref_chr_list, query_chr_list = base.read_collinearity(self.query_name, self.ref_name, self.input_file, chr_list, chr_to_start)

        i = 0
        intra = []
        intra_chr_list = []
        # intra_color_index = []
        for block in data:
            if query_chr_list[i] == ref_chr_list[i]:
                intra.append(block)
                intra_chr_list.append(query_chr_list[i])
                # intra_color_index.append(i)
                i += 1
            else:
                color = chr_color_dict[ref_chr_list[i]]['chr']
                id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                pos1, pos2, pos3, pos4 = gene_pos_dict[id1], gene_pos_dict[id2], gene_pos_dict[id3], gene_pos_dict[id4]
                pos1_radian = self.pos_to_radian(pos1, total_length)
                pos2_radian = self.pos_to_radian(pos2, total_length)
                pos3_radian = self.pos_to_radian(pos3, total_length)
                pos4_radian = self.pos_to_radian(pos4, total_length)

                x, y = self.plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, 0.99*self.inner_radius)
                plt.fill(x, y, facecolor=color, alpha=0.5)
                i += 1
        i = 0
        for block in intra:
            color = chr_color_dict[intra_chr_list[i]]['chr']
            id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
            pos1, pos2, pos3, pos4 = gene_pos_dict[id1], gene_pos_dict[id2], gene_pos_dict[id3], gene_pos_dict[id4]
            pos1_radian = self.pos_to_radian(pos1, total_length)
            pos2_radian = self.pos_to_radian(pos2, total_length)
            pos3_radian = self.pos_to_radian(pos3, total_length)
            pos4_radian = self.pos_to_radian(pos4, total_length)

            x, y = self.plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, 0.99 * self.inner_radius)
            plt.fill(x, y, facecolor=color, alpha=0.3)
            i += 1
        
        # First
        legend_x, legend_y = self.plot_legend_first()
        plt.fill(legend_x, legend_y, facecolor='red', alpha=0.5)
        plt.text(self.inner_radius + 0.5 * self.ring_width, self.inner_radius - 0.5 * self.ring_width, 
                 self.ref_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black')
        
        # Second
        if self.ref_name != self.query_name:
            legend_x, legend_y = self.plot_legend_second()
            plt.fill(legend_x, legend_y, facecolor='blue', alpha=0.5)
            plt.text(self.inner_radius + 0.5 * self.ring_width, self.inner_radius - 2.5 * self.ring_width, 
                     self.query_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black')
        plt.axis('off')
        # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.output_file_name, dpi=int(self.dpi), bbox_inches='tight')
        sys.exit(0)
