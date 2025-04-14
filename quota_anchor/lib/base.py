import os
import sys
import logging
import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger('main.base')
class FileEmptyError(Exception):
    pass


class ClsVis:
    def __init__(self, stats_file, figure, flag):
        self.stats_file = stats_file
        self.output = figure
        self.ylab = "Number"
        self.flag = flag
        if not self.flag:
            self.title1 = "Different gene types of focal species"
            self.title2 = "Different duplicate pair types of focal species"
        else:
            self.title1 = "Different gene types of focal species(Unique)"
            self.title2 = "Different duplicate pair types of focal species(Unique)"

    @staticmethod
    def get_data(stats_file):
        df = pd.read_csv(stats_file, sep="\t", index_col="Type", header=0)
        df_pair = df.loc[["wgd.pairs", "tandem.pairs", "proximal.pairs", "transposed.pairs", "dispersed.pairs"], :]
        df_gene = df.loc[["wgd.genes", "tandem.genes", "proximal.genes", "transposed.genes", "dispersed.genes", "singleton.genes"], :]
        return df_pair, df_gene

    def plot(self, df, name):
        fig, ax = plt.subplots(figsize=(10, 6.18))
        tp = [ele for ele in df.index]
        length = len(tp)
        counts = [df.loc[tp[i], "Number"] for i in range(length)]
        color_list = ['#FFB6C1', '#C71585', '#FF00FF', '#9932CC', '#00CED1', "#90EE90"]

        bars = ax.bar(tp, counts, color=color_list[:length])
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, height, f'{height}', ha='center', va='bottom')

        ax.set_ylabel(str(self.ylab))

        if len(df) == 6:
            ax.set_title(str(self.title1))
        if len(df) == 5:
            ax.set_title(str(self.title2))
        plt.xticks(rotation=300)
        plt.tight_layout()
        plt.savefig(name,  bbox_inches='tight')
#        plt.show()

    def run(self):
        output_list = self.output
        df_pair, df_gene = self.get_data(self.stats_file)
        self.plot(df_gene, output_list[0])
        self.plot(df_pair, output_list[1])

def split_conf(conf, separator):
    new_conf = []
    split_lt = conf.split(separator)
    for ele in split_lt:
        ele = ele.strip()
        if len(ele) == 0:
            continue
        new_conf.append(ele)
    return new_conf

def read_collinearity(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start):
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
                    gene_pos_dict[line_list[0]] = chr_to_start[ref_chr] + int(line_list[4])
                    gene_pos_dict[line_list[5]] = chr_to_start[query_chr] + int(line_list[9])
                else:
                    continue
        if block:
            data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
    return data, gene_pos_dict, ref_chr_list, query_chr_list, direction_list

     
def bezier3(a, b, c, d, t):
    fvalue = []
    for i in t:
        y = a * (1-i)**3 + 3 * b * i * (1-i)**2 + 3 * c * (1-i) * i**2 + d * i**3
        fvalue.append(y)
    return fvalue


def file_empty(file_path):
    try:
        file_path = os.path.abspath(file_path)
        if os.path.exists(file_path):
            if os.path.getsize(file_path) > 0:
                pass
            else:
                raise FileEmptyError
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file."
        logger.error(exist_error_message)
        sys.exit(1)
    except FileEmptyError:
        empty_error_message = "{0} is empty, please check your parameters and files.".format(file_path)
        logger.error(empty_error_message)
        sys.exit(1)
    except OSError:
        ose_error_message = f"{file_path} does not exist or is inaccessible."
        logger.error(ose_error_message)
        sys.exit(1)


def output_file_parentdir_exist(path, overwrite):
    if os.path.exists(path):
        logger.info(f"Output file {path} already exist.")
        if overwrite:
            os.remove(path)
            logger.info(f"Output file {path} will be overwrote.")
            return
        else:
            logger.info(f"Output file {path} will not be overwrote, and you can set '--overwrite' in the command line to overwrite it.")
            sys.exit(1)
    path = os.path.abspath(path)
    dir_name = os.path.dirname(path)
    if os.path.isdir(dir_name):
        pass
    else:
        logger.info(f"{dir_name} does not exist and the software will recursively create the directory.")
        os.makedirs(dir_name, exist_ok=True)

def get_blank_chr(query_df, ref_df):
    chr_margin = pd.merge(left=ref_df, right=query_df, how="cross")
    chr_margin.columns = ["refChr", "referenceStart",  "refId", "queryChr", "queryStart", "queryId"]
    chr_margin_copy =  chr_margin.copy()
    length = len(chr_margin)
    chr_margin_copy["referenceStart"] = [0] *  length
    chr_margin_copy["queryStart"] = [0] *  length
    chr_margin_copy["refId"] = [0] *  length
    chr_margin_copy["queryId"] = [0] *  length
    chr_margin = pd.concat([chr_margin, chr_margin_copy], axis=0)
    return chr_margin

def output_info(info):
    if info:
        info_list = info.split("\n")
        for i in info_list:
            if i:
                logger.info(f'{i}')
    return

