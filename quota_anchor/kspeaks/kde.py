from ..lib import base
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

logger = logging.getLogger('main.kde')
class Kde:
    def __init__(self, config_pra, parameter):
        self.ref_length = ""
        self.query_length = ""
        self.ks_file = ""
        self.collinearity_file = ""
        self.output_file = ""
        self.overwrite = False

        self.ks_range = "0,3"
        self.figsize = "10,6.18"
        self.bins_number = 100
        for i in config_pra.sections():
            if i == 'kde':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.col = "ks_NG86"
        self.bins_number = int(self.bins_number)
        
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
    def get_block_all_ks_and_median(ks_dict, collinearity_dict):
        block_dict = {}
        median_one = {}
        last_key = -1
        for key, values in collinearity_dict.items():
            block_dict[key] = []
            if last_key != -1:
                if len(block_dict[last_key]) >= 5:
                    median_one[last_key] = np.median(block_dict[last_key])
                else:
                    block_dict.pop(last_key)
            for pair_name in values:
                try:
                    ks = ks_dict[pair_name]
                except KeyError:
                    # non self vs self may be keyError
                    elements = pair_name.split("_")
                    new_pair_name = elements[1] + "_" + elements[0]
                    if new_pair_name in ks_dict:
                        ks = ks_dict[new_pair_name]
                    else:
                        continue
                block_dict[key].append(ks)
            last_key = key
        if len(block_dict[last_key]) >= 5:
            median_one[last_key] = np.median(block_dict[last_key])
        else:
            block_dict.pop(last_key)
        return block_dict, median_one

    def get_all_median_list(self, ks_min, ks_max, ref_chr_list: list[str], query_chr_list: list[str]):
        ks_dict = self.read_ks_get_map(self.ks_file, self.col, ks_min, ks_max)
        collinearity_dict = self.read_collinearity(self.collinearity_file, ref_chr_list, query_chr_list)
        block_dict_ks, median_one = self.get_block_all_ks_and_median(ks_dict, collinearity_dict)
        median_list = list(median_one.values())
        all_list = []
        for _, values in block_dict_ks.items():
            all_list.extend(values)
        
        return all_list, median_list            

    @staticmethod
    def ks_kde(all_list, median_list, modifier):
        all_list_kde = gaussian_kde(all_list, bw_method='scott')
        all_list_kde.set_bandwidth(bw_method=all_list_kde.factor * modifier)
        median_list_kde = gaussian_kde(median_list, bw_method='silverman')
        median_list_kde.set_bandwidth(bw_method=median_list_kde.factor * modifier)
        return all_list_kde, median_list_kde

    def run(self):
        logger.info("Init kde and the following parameters are config information")
        print()
        for key, value in vars(self).items():
            if key != "conf" and key != "col":
                print(key, "=", value)
        print()
        logger.info("Check if the input file exists and is not empty")
        base.file_empty(self.ref_length)
        base.file_empty(self.query_length)
        base.file_empty(self.ks_file)
        base.file_empty(self.collinearity_file)
        base.output_file_parentdir_exist(self.output_file, self.overwrite)
        
        ref_chr_list = pd.read_csv(self.ref_length, sep="\t", header=0, index_col=None)['chr'].astype(str).tolist()
        query_chr_list = pd.read_csv(self.query_length, sep="\t", header=0, index_col=None)['chr'].astype(str).tolist()
        self.figsize = [float(k) for k in self.figsize.split(',')]
        ks_min, ks_max = float(self.ks_range.split(",")[0]), float(self.ks_range.split(",")[1])
        all_list_raw, median_list_raw = self.get_all_median_list(ks_min, ks_max, ref_chr_list, query_chr_list)
        plt.rcParams['ytick.major.pad'] = 0
        plt.rcParams['mathtext.default'] = 'regular'
        plt.rcParams['font.size'] = 18
        fig, ax = plt.subplots(figsize=self.figsize)

        _, bins, _ = ax.hist(all_list_raw, int(self.bins_number), density=True, label='histogram of all pairs', histtype="stepfilled", alpha=0.3, fc='grey')
        _, bins, _ = ax.hist(median_list_raw, int(self.bins_number), density=True, label='histogram of median', histtype="stepfilled", alpha=0.3, fc='red')      
        # modifier parameter: more larger more smooth
        all_list_kde, median_list_kde = self.ks_kde(all_list_raw, median_list_raw, modifier=1)
        x = np.linspace(ks_min, ks_max, self.bins_number)
        # bins_width = (ks_max-ks_min)/self.bins_number
        ax.plot(x, all_list_kde(x), color='#0000FF', label='all pairs', linewidth=1.5)
        ax.plot(x, median_list_kde(x), color='#C71585', label='block median', linewidth=1.5)      
        
        ax.legend(frameon=False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_position(('outward', 10))
        ax.set_xlim([ks_min, ks_max])
        ax.set_xlabel(r'${K_s}$')
        ax.set_ylabel('Density')
        plt.setp(ax.yaxis.get_majorticklabels(), rotation=90, va='center')
        plt.savefig(self.output_file, dpi=500, bbox_inches='tight', transparent=True)
        plt.show()
        logger.info("kde plot finished!")