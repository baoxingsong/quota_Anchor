import pandas as pd
import seaborn as sns
import logging
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
# from matplotlib.transforms import Bbox
import numpy as np
from math import pi
import sys
import datetime
from alive_progress import alive_bar
from . import base

GAP_RATIO = 0.25
DPI = 500
# 0.31 = 0.3 + 0.01
INNER_RADIUS = 0.3
OUTER_RADIUS = 0.31
CAP_DIAMETER = 0.01

logger = logging.getLogger('main.circle')
class Circle:
    def __init__(self, config_pra, parameter):
        self.overwrite = False
        self.ref_name = "Reference_species"
        self.query_name = "Query_species"
        self.remove_chromosome_prefix = "CHR,Chr,chr"
        self.chr_font_size = 12
        self.figsize = "14,14"
        self.species_name_font_size = 12
        self.italic = False
        self.input_file = ""
        self.ref_length = ""
        self.query_length = ""
        self.output_file_name = ""
        for i in config_pra.sections():
            if i == 'circle':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.species_name_font_size = float(self.species_name_font_size)
        self.chr_font_size = float(self.chr_font_size)

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
        chr_pos = list(zip(start_list, end_list))
        return chr_pos, chr_to_start

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

    @staticmethod
    def plot_legend_first() -> tuple[list, list]:
        x = [INNER_RADIUS, INNER_RADIUS, INNER_RADIUS - CAP_DIAMETER * 1.618, INNER_RADIUS - CAP_DIAMETER *1.618]
        y = [INNER_RADIUS - CAP_DIAMETER, INNER_RADIUS, INNER_RADIUS, INNER_RADIUS - CAP_DIAMETER]
        return x, y

    @staticmethod
    def plot_legend_second() -> tuple[list, list]:
        x = [INNER_RADIUS, INNER_RADIUS, INNER_RADIUS - CAP_DIAMETER *1.618, INNER_RADIUS - CAP_DIAMETER *1.618]
        y = [INNER_RADIUS - 3 * CAP_DIAMETER, INNER_RADIUS - 2 * CAP_DIAMETER, 
             INNER_RADIUS - 2 * CAP_DIAMETER, INNER_RADIUS - 3 * CAP_DIAMETER]
        return x, y

    @staticmethod
    def plot_circle_chr(start_radian, end_radian, radius_a, radius_b, cap_diameter, ring_radian):
        # here need 2 * ring radian < rectangle
        new_start_radian = start_radian + ring_radian
        new_end_radian = end_radian - ring_radian
        t = np.arange(new_start_radian, new_end_radian, pi / 618)
        x = (radius_a * np.cos(t)).tolist()
        y = (radius_a * np.sin(t)).tolist()

        t = np.arange(new_end_radian+pi, new_end_radian, -pi / 30)
        x1 = ((radius_a + radius_b) / 2 * np.cos(new_end_radian) + cap_diameter / 2 * np.cos(t)).tolist()
        y1 = ((radius_a + radius_b) / 2 * np.sin(new_end_radian) + cap_diameter / 2 * np.sin(t)).tolist()
        x += x1
        y += y1

        t = np.arange(new_end_radian, new_start_radian, -pi / 618)
        x += (radius_b * np.cos(t)).tolist()
        y += (radius_b * np.sin(t)).tolist()

        t = np.arange(new_start_radian, new_start_radian-pi, -pi / 30)
        x1 = ((radius_a + radius_b) / 2 * np.cos(new_start_radian) + cap_diameter / 2 * np.cos(t)).tolist()
        y1 = ((radius_a + radius_b) / 2 * np.sin(new_start_radian) + cap_diameter / 2 * np.sin(t)).tolist()
        x += x1
        y += y1
        start_x = x[0]
        start_y = y[0]
        x.append(start_x)
        y.append(start_y)
        return x, y

    @staticmethod
    def plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, radius):
        ratio = 0.382
        # pos1 -> pos3(ref region)
        t = np.arange(pos1_radian, pos3_radian, pi / 618)
        x = list(radius * np.cos(t))
        y = list(radius * np.sin(t))

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
        t = np.arange(min(pos4_radian, pos2_radian), max(pos4_radian, pos2_radian), pi / 618)
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

    def circle_init(self):
        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)
        print()

        base.file_empty(self.input_file)
        base.file_empty(self.ref_length)
        base.file_empty(self.query_length)
        base.output_file_parentdir_exist(self.output_file_name, self.overwrite)

        if not self.query_length:
            logger.error("Please specify your query species chromosome length file")
            sys.exit(1)
        if not self.ref_length:
            logger.error("Please specify your reference species chromosome length file")
            sys.exit(1)
        if not self.input_file:
            logger.error("Please specify your collinearity file as input file")
            sys.exit(1)
        if not self.output_file_name:
            logger.error("Please specify your output file name")
            sys.exit(1)
        chr_abbr = self.remove_chromosome_prefix.split(',')
        strip_chr_abbr = []
        for i in chr_abbr:
            if len(i) == 0:
                continue
            i = i.strip()
            strip_chr_abbr.append(i)

        width, height = self.figsize.split(",")
        return strip_chr_abbr, width, height

    @staticmethod
    def read_length(name, length):
        length = pd.read_csv(length, sep='\t', header=0, index_col=None)
        length['chr'] = length['chr'].astype(str)
        length['chr'] = name + length['chr']
        length['length'] = length['length'].astype(int)
        return length

    def get_chr_info(self):
        ref_length = self.read_length(self.ref_name, self.ref_length)
        query_length = self.read_length(self.query_name, self.query_length)
        df_dup = pd.concat([ref_length, query_length], axis=0)
        df_dup.drop_duplicates(inplace=True)
        df_dup.reset_index(inplace=True, drop=True)
        return df_dup

    def get_block_chr_color(self, chr_list):
        cmap = self.set_husl_palette()
        chr_color_dict = {}
        i = 1
        for ch in chr_list:
            color_dict = {'chr': cmap(round(i / len(chr_list), 2))}
            chr_color_dict[ch] = color_dict
            i += 1
        return chr_color_dict

    @staticmethod
    def init_figure(width, height):
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        fig, ax = plt.subplots(figsize=(float(width), float(height)), facecolor='white')
        # fig.patch.set_alpha(0.5)
        ax.set_aspect('equal')
        return fig, ax

    @staticmethod
    def get_circle_info(df_dup, chr_list):
        total_chr_length = df_dup['length'].sum()
        average_length = total_chr_length / len(chr_list)
        gap_length = int(average_length * GAP_RATIO)
        circle_perimeter = total_chr_length + len(chr_list) * gap_length
        return circle_perimeter, gap_length

    def run(self):
        logger.info("Circle module init and the following parameters are config information.")
        strip_chr_abbr, width, height = self.circle_init()
        df_dup = self.get_chr_info()

        # intra 10 inter 20 for zm sb
        chr_list = df_dup['chr'].tolist()

        # husl_palette and set line, chr color
        chr_color_dict = self.get_block_chr_color(chr_list)

        # figure geometry info
        fig, ax = self.init_figure(width, height)

        circle_perimeter, gap_length = self.get_circle_info(df_dup, chr_list)
        # bp level
        chr_pos, chr_to_start = self.get_pos_list(df_dup, gap_length)
        # radian level
        chr_radian_pos = self.length_to_radian(chr_pos, circle_perimeter)

        # circle cap not bending
        i = 0
        ring_radian = np.arctan(CAP_DIAMETER / 2 / (INNER_RADIUS + CAP_DIAMETER / 2))
        for ch in chr_list:
            x, y = self.plot_circle_chr(chr_radian_pos[i][0], chr_radian_pos[i][1], INNER_RADIUS, OUTER_RADIUS, CAP_DIAMETER, ring_radian)
            # color = chr_color_dict[ch]['chr']
            if not ch.startswith(str(self.ref_name)) and ch.startswith(str(self.query_name)):
                color = "blue"
                ch = ch[len(self.query_name):]
            else:
                color = "red"
                ch = ch[len(self.ref_name):]
            plt.fill(x, y, facecolor=color, alpha=.5)
            # plt.plot(x, y, color='black')
            label_x = OUTER_RADIUS * 1.04 * np.cos((chr_radian_pos[i][0] + chr_radian_pos[i][1]) / 2)
            label_y = OUTER_RADIUS * 1.04 * np.sin((chr_radian_pos[i][0] + chr_radian_pos[i][1]) / 2)
            angle = self.text_rotation(chr_pos[i][0], chr_pos[i][1], circle_perimeter)
            for abbr in strip_chr_abbr:
                if ch.startswith(abbr):
                    ch = ch[len(abbr):]
            plt.text(label_x, label_y, ch, ha="center", va="center", fontsize=self.chr_font_size, color='black', rotation=angle)
            i += 1
        try:
            data, gene_pos_dict, ref_chr_list, query_chr_list, _ = base.read_collinearity(self.query_name, self.ref_name, self.input_file, chr_list, chr_to_start)
            assert len(data) > 0 and len(ref_chr_list) > 0 and len(query_chr_list) > 0 and len(gene_pos_dict) > 0
        except AssertionError:
            logger.error(f'Please check your {self.input_file}, {self.ref_length} and {self.query_length}.')
            sys.exit(1)

        i = 0
        intra = []
        intra_chr_list = []
        # intra_color_index = []
        dt = datetime.datetime.now()
        time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
        with alive_bar(len(data), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
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
                    pos1_radian = self.pos_to_radian(pos1, circle_perimeter)
                    pos2_radian = self.pos_to_radian(pos2, circle_perimeter)
                    pos3_radian = self.pos_to_radian(pos3, circle_perimeter)
                    pos4_radian = self.pos_to_radian(pos4, circle_perimeter)

                    x, y = self.plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, 0.99*INNER_RADIUS)
                    plt.fill(x, y, facecolor=color, alpha=0.5)
                    i += 1
                bar()
        if len(intra) > 0:
            i = 0
            dt = datetime.datetime.now()
            time_now = dt.strftime('%Y/%m/%d %H:%M:%S')
            with alive_bar(len(intra), title=f"[{time_now} INFO]", bar="bubbles", spinner="waves") as bar:
                for block in intra:
                    color = chr_color_dict[intra_chr_list[i]]['chr']
                    id1, id2, id3, id4 = block[0], block[1], block[2], block[3]
                    pos1, pos2, pos3, pos4 = gene_pos_dict[id1], gene_pos_dict[id2], gene_pos_dict[id3], gene_pos_dict[id4]
                    pos1_radian = self.pos_to_radian(pos1, circle_perimeter)
                    pos2_radian = self.pos_to_radian(pos2, circle_perimeter)
                    pos3_radian = self.pos_to_radian(pos3, circle_perimeter)
                    pos4_radian = self.pos_to_radian(pos4, circle_perimeter)

                    x, y = self.plot_closed_region(pos1_radian, pos2_radian, pos3_radian, pos4_radian, 0.99 * INNER_RADIUS)
                    plt.fill(x, y, facecolor=color, alpha=0.8)
                    i += 1
                    bar()
        
        # First
        legend_x, legend_y = self.plot_legend_first()
        plt.fill(legend_x, legend_y, facecolor='red', alpha=0.5)
        if self.italic:
            plt.text(INNER_RADIUS + 0.5 * CAP_DIAMETER, INNER_RADIUS - 0.5 * CAP_DIAMETER,
                     self.ref_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black', fontstyle='italic')
        else:
            plt.text(INNER_RADIUS + 0.5 * CAP_DIAMETER, INNER_RADIUS - 0.5 * CAP_DIAMETER,
                     self.ref_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black')
        
        # Second
        if self.ref_name != self.query_name:
            legend_x, legend_y = self.plot_legend_second()
            plt.fill(legend_x, legend_y, facecolor='blue', alpha=0.5)
            if self.italic:
                plt.text(INNER_RADIUS + 0.5 * CAP_DIAMETER, INNER_RADIUS - 2.5 * CAP_DIAMETER,
                         self.query_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black', fontstyle='italic')
            else:
                plt.text(INNER_RADIUS + 0.5 * CAP_DIAMETER, INNER_RADIUS - 2.5 * CAP_DIAMETER,
                         self.query_name, ha="left", va="center", fontsize=self.species_name_font_size, color='black')
        plt.axis('off')
        # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.output_file_name, dpi=DPI, bbox_inches='tight', transparent=False)
        logger.info(f"Plot {self.output_file_name} finished!")
        sys.exit(0)

