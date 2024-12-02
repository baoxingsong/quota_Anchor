import sys
from ..lib import base
import random
import logging
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings('ignore', category=SyntaxWarning)

logger = logging.getLogger('main.correct')

class Correct:
    def __init__(self, config_pra, parameter):
        self.species_pair_file = ""
        self.trios_file = ""
        self.species_pair_ks_file = ""
        self.outfile_divergent_peaks = "outfile_divergent_peaks.csv"
        self.overwrite = False

        self.ks_range = "0,3"

        for i in config_pra.sections():
            if i == 'correct':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.col = "ks_NG86"
        self.bins_number = 200
        self.n_replicates = 382

    def get_ks_limit(self):
        ks_range = self.ks_range.split(",")
        ks_min, ks_max = float(ks_range[0]), float(ks_range[1])
        return ks_min, ks_max

    def init_process(self):
        print()
        for key, value in vars(self).items():
            if key != "conf" and key not in ["col", "bins_number", "n_replicates"]:
                print(key, "=", value)
        print()

        if not self.species_pair_file:
            logger.error("Please specify your species pair file path")
            sys.exit(1)
        if not self.trios_file:
            logger.error("Please specify your trios file path")
        if not self.species_pair_ks_file:
            logger.error(r'Please specify your species pair ks file path(Separator ",")')
        if not self.outfile_divergent_peaks:
            logger.error("Please specify your output file path")
            sys.exit(1)

        base.file_empty(self.species_pair_file)
        base.file_empty(self.trios_file)
        base.output_file_parentdir_exist(self.outfile_divergent_peaks, self.overwrite)

        pair_df = pd.read_csv(self.species_pair_file, header=0, index_col=None)
        pair_df["Species_1"] = pair_df["Species_1"].astype(str)
        pair_df["Species_2"] = pair_df["Species_2"].astype(str)

        ks_file_list = base.split_conf(self.species_pair_ks_file, ",")
        try:
            assert len(ks_file_list) == len(pair_df)    
        except AssertionError:
            logger.error("The number of ks file is not equal to the number of species pairs")
            sys.exit(1)
        return pair_df, ks_file_list
    
    @staticmethod
    def get_pair_ks_map(pair_df, ks_file_list):
        pair_ks_dict = {}
        pair_list = pair_df.apply(lambda row: str(row['Species_1']) + '_' + str(row['Species_2']), axis=1).tolist()
        species_pair_to_ks_file = zip(pair_list, ks_file_list)
        for pr, ks in species_pair_to_ks_file:
            pair_ks_dict[pr] = ks
        return pair_ks_dict
    
    def compute_kde(self, ks_list, ks_min, ks_max):
        kde = gaussian_kde(ks_list, bw_method='scott')
        variable_bw = kde.factor * np.random.uniform(0.5, 0.9)
        kde.set_bandwidth(variable_bw)
        kde_x = np.linspace(ks_min, ks_max, int(self.bins_number))
        kde_y = kde(kde_x)
        return kde, kde_x, kde_y

    def get_mode_sd(self, ks_list, ks_min, ks_max ,n_replicates):
        mode_list = []
        median_list = []
        for _ in range(n_replicates):
            sample_list = random.choices(ks_list, k=len(ks_list))
            __, kde_x, kde_y = self.compute_kde(sample_list, ks_min, ks_max)

            index_max_y = np.argmax(kde_y)
            mode_x_value = kde_x[index_max_y]

            mode_list.append(mode_x_value)
            median_list.append(np.median(sample_list))

        mean_mode = np.mean(mode_list, dtype=np.float64)
        std_mode = np.std(mode_list, dtype=np.float64)

        mean_median = np.mean(median_list, dtype=np.float64)
        std_median = np.std(median_list, dtype=np.float64)
        return mean_mode, std_mode, mean_median, std_median
    
    def get_fitting_info(self):
        logger.info(r"The order of species pairs in the species pair file(specify by -s parameter/species_pair_file) must be consistent with the order of the ks file(specify by -k parameter/species_pair_ks_file)")
        pair_df, ks_file_list = self.init_process()
        ks_min, ks_max = self.get_ks_limit()
        col = self.col
        pair_ks_dict = self.get_pair_ks_map(pair_df, ks_file_list)
        
        mode_sd = []
        for pair in pair_ks_dict:
            df = pd.read_csv(pair_ks_dict[pair], header=0, index_col=None, sep="\t")
            df['id1'] = df['id1'].astype(str)
            df['id2'] = df['id2'].astype(str)
            ks_df = df[(df[col] >= max(ks_min, 0)) & (df[col] <= ks_max) & (df[col] != -1)]
            ks_list = ks_df[col].tolist()
            logger.info(f"{pair} kernel density estimation start")
            mode_mean, mode_sd_mean, median_mean, median_sd_mean = self.get_mode_sd(ks_list, ks_min, ks_max, self.n_replicates)
            logger.info(f"{pair} kernel density estimation end")
            mode_sd.append([mode_mean, mode_sd_mean])

        peaks_info_df = pd.DataFrame(mode_sd, columns=["mode", "sd"])
        pair_mode_info_df = pd.concat([pair_df, peaks_info_df], axis=1)
        return pair_mode_info_df

    @staticmethod
    def process_pair_info_df(df):
        df['key'] = df['Species_1'] + '_' + df['Species_2']
        df['reverse_key'] = df['Species_2'] + '_' + df['Species_1']
        # get mode map
        df_keys = pd.concat([df[['key', "mode"]], df[['reverse_key', "mode"]].rename(columns={'reverse_key': 'key'})])
        df_keys.drop_duplicates(subset=["key"] ,inplace=True)
        mode_map = df_keys.set_index('key')["mode"].to_dict()
        
        # get sd map
        df_keys = pd.concat([df[['key', "sd"]], df[['reverse_key', "sd"]].rename(columns={'reverse_key': 'key'})])
        df_keys.drop_duplicates(subset=["key"] ,inplace=True)
        sd_map = df_keys.set_index('key')["sd"].to_dict()

        return mode_map, sd_map

    def correct_trios(self, pair_mode_info_df):
        sisters_per_node = {}
        pre_strategy = []
        mode_map, sd_map = self.process_pair_info_df(pair_mode_info_df)

        trios_df = pd.read_csv(self.trios_file, header=0, index_col=None)
        trios_df['Focal_Species'] = trios_df['Focal_Species'].astype(str)
        trios_df['Sister_Species'] = trios_df['Sister_Species'].astype(str)
        trios_df['Out_Species'] = trios_df['Out_Species'].astype(str)
        
        for index, row in trios_df.iterrows():
            node = row['Node']
            focal_name, sister_name, out_name = row['Focal_Species'], row['Sister_Species'], row['Out_Species']
            focal_sister_name = focal_name + "_" + sister_name
            focal_out_name = focal_name + "_" + out_name
            sister_out_name = sister_name + "_" + out_name
            # raw_mode and sd
            focal_sister_mode, focal_sister_sd = mode_map[focal_sister_name], sd_map[focal_sister_name]
            focal_out_mode, focal_out_sd = mode_map[focal_out_name], sd_map[focal_out_name]
            sister_out_mode, sister_out_sd = mode_map[sister_out_name], sd_map[sister_out_name]

            adjusted_focal = (focal_out_mode + focal_sister_mode - sister_out_mode) / 2.0
            adjusted_sister = (focal_sister_mode + sister_out_mode - focal_out_mode) / 2.0
            # Error propagation rules
            adjusted_focal_sd = np.sqrt(pow(focal_out_sd, 2) + pow(focal_sister_sd, 2) + pow(sister_out_sd, 2)) / 2.0

            adjusted_focal_sister = 2 * adjusted_focal
            adjusted_focal_sister_sd = np.sqrt(2) * adjusted_focal_sd
            # smaller and more accurate
            OC_distance = focal_out_mode - adjusted_focal

            if node not in sisters_per_node:
                sisters_per_node[node] = []
            if sister_name not in sisters_per_node[node]:
                sisters_per_node[node].append(sister_name)
            
            pre_strategy.append([node, focal_name, sister_name, out_name, 
                                 adjusted_focal_sister, adjusted_focal_sister_sd, focal_sister_mode, focal_sister_sd, 
                                 adjusted_focal, adjusted_sister, OC_distance])
        all_trios_correction_df = pd.DataFrame(pre_strategy, columns=["Node", "Focal_Species",
                                "Sister_Species", "Out_Species", "Adjusted_Mode", "Adjusted_Mode_SD", "Original_Mode",
                                "Original_Mode_SD", "Ks_Focal", "Ks_Sister", "Out_Distance"])
        return all_trios_correction_df, sisters_per_node
    
    def mean_or_best_strategy(self, pre_select_trios):
        pair_corrected_list = []
        for name, group in pre_select_trios.groupby(by=['Node', "Sister_Species"]):
            # name (1, 'sorghum') (Node, Sister_Species)
            # strategy one: mean
            adjusted_mode_mean = float(group['Adjusted_Mode'].values.mean())
            focal_mean = float(group['Ks_Focal'].values.mean())
            sister_mean = float(group['Ks_Sister'].values.mean())

            adjusted_sd_list = group['Adjusted_Mode_SD'].to_list()
            sd_err_prop = 0
            for sd in adjusted_sd_list:
                sd_err_prop += pow(sd, 2)
            sd_err_prop = np.sqrt(sd_err_prop) / len(adjusted_sd_list)

            # strategy two: best out
            best_out = group['Out_Distance'].min()
            node_df_sister_best = group[group['Out_Distance'] == best_out]
            mode_mean_best = float(node_df_sister_best['Adjusted_Mode'].values.mean())
            focal_mean_best = float(node_df_sister_best['Ks_Focal'].values.mean())
            sister_mean_best = float(node_df_sister_best['Ks_Sister'].values.mean())

            best_sd_err_prop = 0
            best_adjusted_sd_list = node_df_sister_best['Adjusted_Mode_SD'].values
            for sd in best_adjusted_sd_list:
                best_sd_err_prop += pow(sd, 2)
            sd_mean_best = np.sqrt(best_sd_err_prop) / len(best_adjusted_sd_list)

            pair_corrected_list.append([name[0], group["Focal_Species"].to_list()[0], group["Sister_Species"].to_list()[0], round(adjusted_mode_mean,6), round(sd_err_prop, 6),
                                        round(focal_mean, 6), round(sister_mean, 6), round(mode_mean_best, 6), round(sd_mean_best, 6), round(focal_mean_best, 6), round(sister_mean_best, 6),
                                        round(group["Original_Mode"].to_list()[0], 6), round(group["Original_Mode_SD"].to_list()[0],6)])
        
        df = pd.DataFrame(pair_corrected_list, columns=["Node", "Focal_Species", "Sister_Species",
                                                        "Adjusted_Mode_Mean", "Adjusted_Mode_Mean_SD", "Ks_Focal_Mean", "Ks_Sister_Mean",
                                                        "Adjusted_Mode_Best", "Adjusted_Mode_Best_SD", "Ks_Focal_Best", "Ks_Sister_Best", 
                                                        "Original_Mode", "Original_Mode_SD"])
        df.to_csv(self.outfile_divergent_peaks, header=True, index=False)
        return df

    def run(self):
        logger.info("Trios module init and the following parameters are config information")
        pair_mode_info_df = self.get_fitting_info()

        logger.info("Divergent peak based on ks between focal species and sister species correction start")
        pre_select_trios, sisters_per_node = self.correct_trios(pair_mode_info_df)
        pair_corrected = self.mean_or_best_strategy(pre_select_trios)
        logger.info("Divergent peak based on ks between focal species and sister species correction end")
        logger.info(f"Generate {self.outfile_divergent_peaks} done!")
