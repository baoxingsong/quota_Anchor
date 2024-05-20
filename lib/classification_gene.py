import os
import pandas as pd
from lib import GffFile


def file_exist(path):
    if not os.path.isfile(path):
        error_message = f"{path} don't exist"
        print(error_message)
        raise FileNotFoundError(error_message)


def dir_exist(path):
    if os.path.isdir(path):
        pass
    else:
        os.makedirs(path, exist_ok=True)


class ClassGene:
    def __init__(self, config_pra):
        self.qy_table = str(config_pra['classification']['query_table'])
        self.gff = str(config_pra['classification']['query_gff_file'])
        self.qy_col_fl = str(config_pra["classification"]["query_query_collinearity"])
        self.qy_rf_col = str(config_pra["classification"]["query_ref_collinearity"])
        self.out_dir = str(config_pra["classification"]["out_directory"])
        self.out_prefix = str(config_pra["classification"]["out_prefix"])
        self.pm_threshold = int(config_pra["classification"]["proximal_max_distance"])
        self.seg_anc = int(config_pra["classification"]["seg_anc"])
        # self.type = int(config_pra["classification"]["type"])

        self.wgd_gp_file = self.out_dir + "/" + self.out_prefix + "." + "wgd.pairs"
        self.wgd_gn_file = self.out_dir + "/" + self.out_prefix + "." + "wgd.genes"
        self.tm_gp_file = self.out_dir + "/" + self.out_prefix + "." + "tandem.pairs"
        self.tm_gn_file = self.out_dir + "/" + self.out_prefix + "." + "tandem.genes"
        self.pm_gp_file = self.out_dir + "/" + self.out_prefix + "." + "proximal.pairs"
        self.pm_gn_file = self.out_dir + "/" + self.out_prefix + "." + "proximal.genes"
        self.trp_gp_file = self.out_dir + "/" + self.out_prefix + "." + "transposed.pairs"
        self.trp_gn_file = self.out_dir + "/" + self.out_prefix + "." + "transposed.genes"
        self.dp_gp_file = self.out_dir + "/" + self.out_prefix + "." + "dispersed.pairs"
        self.dp_gn_file = self.out_dir + "/" + self.out_prefix + "." + "dispersed.genes"
        self.si_gn_file = self.out_dir + "/" + self.out_prefix + "." + "singleton.genes"

    def file_jg(self, qy_table, qy_rf_table, qy_col_fl, qy_rf_col, out_dir):
        try:
            file_exist(qy_table)
            file_exist(qy_rf_table)
            file_exist(qy_col_fl)
            qy_rf_lt = qy_rf_col.split(',')
            for file in qy_rf_lt:
                file = file.strip()
                file_exist(file)
        except FileNotFoundError as e:
            print("file not found:", e)
            exit(1)
        dir_exist(out_dir)

    def wgd_pair(self, qy_col_fl, wgd_gp_file, wgd_gn_file, out_prefix):
        wgd_gn_lt = []
        wgd_pr_lt = []
        wgd_md = {}
        wgd_gp_file = open(wgd_gp_file, 'w')
        wgd_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        wgd_gn_file = open(wgd_gn_file, 'w')
        wgd_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        with open(qy_col_fl) as qy_col_fl:
            _ = next(qy_col_fl)
            _ = next(qy_col_fl)
            for line in qy_col_fl:
                if line.startswith("#"):
                    continue
                col_rd_lt = line.split('\t')
                # wgd_pr_lt -> twod   file -> one
                if str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]) not in wgd_pr_lt and str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]) not in wgd_pr_lt:
                    wgd_gp_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3])
                                      + "\t" + str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8]) + "\n")
                    wgd_pr_lt.append(str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]))
                    wgd_pr_lt.append(str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]))
                if str(col_rd_lt[0]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[0]))
                    wgd_md[str(col_rd_lt[0])] = 1
                    wgd_gn_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3] + "\n"))
                if str(col_rd_lt[5]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[5]))
                    wgd_md[str(col_rd_lt[6])] = 1
                    wgd_gn_file.write(str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8] + "\n"))
        wgd_gn_file.close()
        wgd_gp_file.close()
        return wgd_gn_lt, wgd_pr_lt, wgd_md

    def tm_ad_wgd_pr(self, qy_table, tm_gp_file, tm_gn_file, out_prefix, wgd_gn_lt, wgd_pr_lt):
        tm_gn_lt = []
        tm_gp_lt = []
        homo_gn_md = {}
        tm_gp_file = open(tm_gp_file, 'w')
        tm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tm_gn_file = open(tm_gn_file, 'w')
        tm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        with open(qy_table) as qy_table:
            for line in qy_table:
                tm_rd_lt = line.split('\t')
                homo_gn_md[str(tm_rd_lt[0])] = 0
                homo_gn_md[str(tm_rd_lt[6])] = 0
                tandem_condition1 = tm_rd_lt[0] != tm_rd_lt[6]
                tandem_condition2 = tm_rd_lt[1] == tm_rd_lt[7]
                tandem_condition3 = abs(int(tm_rd_lt[2]) - int(tm_rd_lt[8])) == 1
                if tandem_condition3 and tandem_condition2 and tandem_condition1:
                    # tm_gp_lt -> one
                    if int(tm_rd_lt[2]) - int(tm_rd_lt[8]) < 0:
                        if str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in wgd_pr_lt and str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in tm_gp_lt:
                            tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                            tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                            tm_gp_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3])
                                             + "\t" + str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                             + str(tm_rd_lt[9]) + "\n")
                    else:
                        if str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in wgd_pr_lt and str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in tm_gp_lt:
                            tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                            tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                            tm_gp_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                             + str(tm_rd_lt[9]) + "\t" + str(tm_rd_lt[0]) + "\t" + str(out_prefix)
                                             + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3]) + "\n")
                    if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                        tm_gn_lt.append(str(tm_rd_lt[0]))
                        homo_gn_md[str(tm_rd_lt[0])] = 2
                        tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                    if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                        tm_gn_lt.append(str(tm_rd_lt[6]))
                        homo_gn_md[str(tm_rd_lt[6])] = 2
                        tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
            tm_gn_file.close()
            tm_gp_file.close()
        return homo_gn_md, tm_gp_lt

    def pm_ad_wgd_pr(self, qy_table, pm_gp_file, pm_gn_file, out_prefix, wgd_pr_lt, pm_threshold, homo_gn_md):
        pm_gn_lt = []
        pm_gp_lt = []
        pm_pre_df = []
        pm_gp_file = open(pm_gp_file, 'w')
        pm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        pm_gn_file = open(pm_gn_file, 'w')
        pm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        # with open(qy_table) as qy_table:
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            pm_rd_lt = line.split('\t')
            pm_condition1 = pm_rd_lt[0] != pm_rd_lt[6]
            pm_condition2 = pm_rd_lt[1] == pm_rd_lt[7]
            pm_condition3 = 1 < abs(int(pm_rd_lt[2]) - int(pm_rd_lt[8])) <= pm_threshold
            if pm_condition3 and pm_condition2 and pm_condition1:
                if int(pm_rd_lt[2]) - int(pm_rd_lt[8]) < 0:
                    if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in wgd_pr_lt:
                        pm_pre_df.append([str(pm_rd_lt[0]), str(pm_rd_lt[6]), int(pm_rd_lt[8]) - int(pm_rd_lt[2])])
                else:
                    if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in wgd_pr_lt:
                        pm_pre_df.append([str(pm_rd_lt[6]), str(pm_rd_lt[0]), int(pm_rd_lt[2]) - int(pm_rd_lt[8])])
        qy_tb1.close()
        df = pd.DataFrame(pm_pre_df, columns=['Duplicate1', 'Duplicate2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        pm_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate1']) + "\t" + str(row['Duplicate2']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                pm_rd_lt = line.split("\t")
                if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) in pm_gn_pr and str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in pm_gp_lt:
                    pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                    pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                    pm_gp_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3])
                                     + "\t" + str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                     + str(pm_rd_lt[9]) + "\n")
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and str(pm_rd_lt[0]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[0]))
                        homo_gn_md[str(pm_rd_lt[0])] = 3
                        pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    if homo_gn_md[str(pm_rd_lt[6])] not in [1, 2] and str(pm_rd_lt[6]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[6]))
                        homo_gn_md[str(pm_rd_lt[6])] = 3
                        pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) in pm_gn_pr and str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in pm_gp_lt:
                    pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                    pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                    pm_gp_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                     + str(pm_rd_lt[9]) + "\t" + str(pm_rd_lt[0]) + "\t" + str(out_prefix)
                                     + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3]) + "\n")
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and str(pm_rd_lt[0]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[0]))
                        homo_gn_md[str(pm_rd_lt[0])] = 3
                        pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    if homo_gn_md[str(pm_rd_lt[6])] not in [1, 2] and str(pm_rd_lt[6]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[6]))
                        homo_gn_md[str(pm_rd_lt[6])] = 3
                        pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
            pm_gn_file.close()
            pm_gp_file.close()
            return homo_gn_md, pm_gp_lt

    def trp_ad_wgd_pr(self, wgd_gn_lt, qy_rf_col, qy_table, trp_gp_file, trp_gn_file, out_prefix, seg_anc, pm_threshold, homo_gn_md):
        trp_gn_lt = []
        trp_gp_lt = []
        pre_df = []
        df = pd.read_csv(qy_table, header=None, index_col=False, sep='\t')
        df = df.iloc[:, [0, 6, 12]]
        mask = df.iloc[:, 0] != df.iloc[:, 1]
        df = df[mask]
        idt_dict = {str(row[0]) + "\t" + str(row[6]): row[12] for _, row in df.iterrows()}
        # print(idt_dict)
        # max_values = df.groupby(0)[12].transform('max')
        # print(max_values)
        # df_filtered = df[max_values == df[12]]
        # new_list = df.apply(lambda row: str(row[0]) + "\t" + str(row[6]), axis=1).tolist()
        anc = []
        qy_rf_lt = qy_rf_col.split(',')
        for file in qy_rf_lt:
            file = file.strip()
            with open(file) as fi:
                _ = next(fi)
                _ = next(fi)
                for line in fi:
                    if line.startswith("#"):
                        continue
                    col_rf_rd_lt = line.split('\t')
                    anc.append(col_rf_rd_lt[5])
        if seg_anc == 1:
            anc = anc + wgd_gn_lt
            anc = list(set(anc))
        trp_gp_file = open(trp_gp_file, 'w')
        trp_gp_file.write("Transposed" + "\t" + "Location1" + "\t" + "Parental" + "\t" + "Location2" + "\n")
        trp_gn_file = open(trp_gn_file, 'w')
        trp_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            qy_rd_lt = line.split("\t")
            if qy_rd_lt[1] != qy_rd_lt[7]:
                if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [2, 3]:
                    pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
            else:
                if abs(int(qy_rd_lt[2]) - int(qy_rd_lt[8])) > pm_threshold:
                    if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [2, 3]:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                    if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
        qy_tb1.close()
        df = pd.DataFrame(pre_df, columns=['Duplicate', 'Parental', 'identity'])
        df = df.drop_duplicates()
        max_values = df.groupby('Duplicate')['identity'].transform('max')
        df_filtered = df[max_values == df['identity']]
        trp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate']) + "\t" + str(row['Parental']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                qy_rd_lt = line.split("\t")
                if str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) in trp_gn_pr and str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) not in trp_gp_lt:
                    trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                    trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                    trp_gp_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":"
                                      + str(qy_rd_lt[3]) + "\t" + str(qy_rd_lt[6]) + "\t" + str(out_prefix)
                                      + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                    if homo_gn_md[qy_rd_lt[6]] == 0:
                        homo_gn_md[qy_rd_lt[6]] = 4
                    if str(qy_rd_lt[0]) not in trp_gn_lt:
                        trp_gn_lt.append(str(qy_rd_lt[0]))
                        homo_gn_md[str(qy_rd_lt[0])] = 4
                        trp_gn_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in trp_gn_pr and str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) not in trp_gp_lt:
                    trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                    trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                    trp_gp_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":"
                                      + str(qy_rd_lt[9]) + "\t" + str(qy_rd_lt[0]) + "\t" + str(out_prefix)
                                      + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                    if homo_gn_md[qy_rd_lt[0]] == 0:
                        homo_gn_md[qy_rd_lt[0]] = 4
                    if str(qy_rd_lt[6]) not in trp_gn_lt:
                        trp_gn_lt.append(str(qy_rd_lt[6]))
                        homo_gn_md[str(qy_rd_lt[6])] = 4
                        trp_gn_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
        return homo_gn_md, trp_gp_lt

    def dip_ad_pr(self, qy_table, dip_gp_file, dip_gn_file, out_prefix, homo_gn_md, used_gp_lt):
        dp_gp_lt = []
        dp_gn_lt = []
        dp_pre_df = []
        dip_gp_file = open(dip_gp_file, 'w')
        dip_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dip_gn_file = open(dip_gn_file, 'w')
        dip_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            dp_rd_lt = line.split('\t')
            if dp_rd_lt[0] != dp_rd_lt[6]:
                if str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]) not in used_gp_lt:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(dp_rd_lt[12])])
        df = pd.DataFrame(dp_pre_df, columns=['Duplicate1', 'Duplicate2', 'identity'])
        max_values = df.groupby('Duplicate2')['identity'].transform('max')
        df_filtered = df[max_values == df['identity']]
        dp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate2']) + "\t" + str(row['Duplicate1']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                dp_rd_lt = line.split('\t')
                if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in dp_gn_pr and str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) not in dp_gp_lt:
                    dp_gp_lt.append(str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]))
                    dp_gp_lt.append(str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]))
                    dip_gp_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":"
                                      + str(dp_rd_lt[9]) + "\t" + str(dp_rd_lt[0]) + "\t" + str(out_prefix)
                                      + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                if homo_gn_md[str(dp_rd_lt[0])] not in [1, 2, 3, 4] and str(dp_rd_lt[0]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[0]))
                    dip_gn_file.write(str(dp_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                if homo_gn_md[str(dp_rd_lt[6])] not in [1, 2, 3, 4] and str(dp_rd_lt[6]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[6]))
                    dip_gn_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":" + str(dp_rd_lt[9]) + "\n")
        dip_gp_file.close()
        dip_gn_file.close()
        return dp_gn_lt

    def merge_mode(self, homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md

    def singleton(self, gff, homo_gn_md, out_file, out_prefix):
        chr_gn_dict, chr_gn_list, gn_chr_dict, _ = GffFile.readGff(gff)
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        for ch in chr_gn_list:
            for gn in chr_gn_list[ch]:
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + str(out_prefix) + "-" + str(ch) + ":" + str(chr_gn_dict[ch][gn.name].start) + "\n")
        ot_file.close()

    def sort_file(self, gene_file, pair_file):
        for gn_file in gene_file:
            df = pd.read_csv(gn_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(gn_file, header=True, index=False, sep='\t')
        for pr_file in pair_file:
            df = pd.read_csv(pr_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(pr_file, header=True, index=False, sep='\t')

    def run(self):
        self.file_jg(self.qy_table, self.gff, self.qy_col_fl, self.qy_rf_col, self.out_dir)
        wgd_gn_lt, wgd_pr_lt, wgd_md = self.wgd_pair(self.qy_col_fl, self.wgd_gp_file, self.wgd_gn_file, self.out_prefix)
        homo_gn_md, tm_gp_lt = self.tm_ad_wgd_pr(self.qy_table, self.tm_gp_file, self.tm_gn_file, self.out_prefix, wgd_gn_lt, wgd_pr_lt)
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, pm_gp_lt = self.pm_ad_wgd_pr(self.qy_table, self.pm_gp_file, self.pm_gn_file, self.out_prefix,
                                                 wgd_pr_lt, self.pm_threshold, homo_gn_md)
        homo_gn_md, trp_gp_lt = self.trp_ad_wgd_pr(wgd_gn_lt, self.qy_rf_col, self.qy_table, self.trp_gp_file, self.trp_gn_file,
                                                   self.out_prefix, self.seg_anc, self.pm_threshold, homo_gn_md)
        used_gp_lt = wgd_pr_lt + tm_gp_lt + pm_gp_lt + trp_gp_lt
        self.dip_ad_pr(self.qy_table, self.dp_gp_file, self.dp_gn_file, self.out_prefix, homo_gn_md, used_gp_lt)
        self.singleton(self.gff, homo_gn_md, self.si_gn_file, self.out_prefix)
        self.sort_file([self.wgd_gn_file, self.tm_gn_file, self.pm_gn_file, self.trp_gn_file, self.dp_gn_file],
                       [self.wgd_gp_file, self.tm_gp_file, self.pm_gp_file, self.trp_gp_file, self.dp_gp_file])


class ClassGeneUnique:
    def __init__(self, config_pra):
        self.qy_table = str(config_pra['classification']['query_table'])
        self.gff = str(config_pra['classification']['query_gff_file'])
        self.qy_col_fl = str(config_pra["classification"]["query_query_collinearity"])
        self.qy_rf_col = str(config_pra["classification"]["query_ref_collinearity"])
        self.out_dir = str(config_pra["classification"]["out_directory"])
        self.out_prefix = str(config_pra["classification"]["out_prefix"])
        self.pm_threshold = int(config_pra["classification"]["proximal_max_distance"])
        self.seg_anc = int(config_pra["classification"]["seg_anc"])
        # self.type = int(config_pra["classification"]["type"])

        self.wgd_gp_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "wgd.pairs"
        self.wgd_gn_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "wgd.genes"
        self.tm_gp_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "tandem.pairs"
        self.tm_gn_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "tandem.genes"
        self.pm_gp_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "proximal.pairs"
        self.pm_gn_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "proximal.genes"
        self.trp_gp_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "transposed.pairs"
        self.trp_gn_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "transposed.genes"
        self.dp_gp_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "dispersed.pairs"
        self.dp_gn_file = self.out_dir + "/" + self.out_prefix + "-unique" + "." + "dispersed.genes"
        self.si_gn_file = self.out_dir + "/" + self.out_prefix + "." + "singleton.genes"

    def file_jg(self, qy_table, qy_rf_table, qy_col_fl, qy_rf_col, out_dir):
        try:
            file_exist(qy_table)
            file_exist(qy_rf_table)
            file_exist(qy_col_fl)
            qy_rf_lt = qy_rf_col.split(',')
            for file in qy_rf_lt:
                file = file.strip()
                file_exist(file)
        except FileNotFoundError as e:
            print("file not found:", e)
            exit(1)
        dir_exist(out_dir)

    def wgd_pair(self, qy_col_fl, wgd_gp_file, wgd_gn_file, out_prefix):
        wgd_gn_lt = []
        wgd_pr_lt = []
        wgd_md = {}
        wgd_gp_file = open(wgd_gp_file, 'w')
        wgd_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        wgd_gn_file = open(wgd_gn_file, 'w')
        wgd_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        with open(qy_col_fl) as qy_col_fl:
            _ = next(qy_col_fl)
            _ = next(qy_col_fl)
            for line in qy_col_fl:
                if line.startswith("#"):
                    continue
                col_rd_lt = line.split('\t')
                # wgd_pr_lt -> twod   file -> one
                if str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]) not in wgd_pr_lt and str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]) not in wgd_pr_lt:
                    wgd_gp_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3])
                                      + "\t" + str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8]) + "\n")
                    wgd_pr_lt.append(str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]))
                    wgd_pr_lt.append(str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]))
                if str(col_rd_lt[0]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[0]))
                    wgd_md[str(col_rd_lt[0])] = 1
                    wgd_gn_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3] + "\n"))
                if str(col_rd_lt[5]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[5]))
                    wgd_md[str(col_rd_lt[6])] = 1
                    wgd_gn_file.write(str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8] + "\n"))
        wgd_gn_file.close()
        wgd_gp_file.close()
        return wgd_gn_lt, wgd_pr_lt, wgd_md

    def tm_ad_wgd_pr(self, qy_table, tm_gp_file, tm_gn_file, out_prefix, wgd_gn_lt):
        tm_gn_lt = []
        tm_gp_lt = []
        homo_gn_md = {}
        tm_gp_file = open(tm_gp_file, 'w')
        tm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tm_gn_file = open(tm_gn_file, 'w')
        tm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        with (open(qy_table) as qy_table):
            for line in qy_table:
                tm_rd_lt = line.split('\t')
                homo_gn_md[str(tm_rd_lt[0])] = 0
                homo_gn_md[str(tm_rd_lt[6])] = 0
                tandem_condition1 = tm_rd_lt[0] != tm_rd_lt[6]
                tandem_condition2 = tm_rd_lt[1] == tm_rd_lt[7]
                tandem_condition3 = abs(int(tm_rd_lt[2]) - int(tm_rd_lt[8])) == 1
                if tandem_condition3 and tandem_condition2 and tandem_condition1:
                    # tm_gp_lt -> one
                    if int(tm_rd_lt[2]) - int(tm_rd_lt[8]) < 0:
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in wgd_gn_lt and \
                               str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in tm_gp_lt:
                            tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                            tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                            tm_gp_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3])
                                             + "\t" + str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                             + str(tm_rd_lt[9]) + "\n")
                    else:
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in wgd_gn_lt and \
                                str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in tm_gp_lt:
                            tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                            tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                            tm_gp_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                             + str(tm_rd_lt[9]) + "\t" + str(tm_rd_lt[0]) + "\t" + str(out_prefix)
                                             + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3]) + "\n")
                    if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                        tm_gn_lt.append(str(tm_rd_lt[0]))
                        homo_gn_md[str(tm_rd_lt[0])] = 2
                        tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                    if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                        tm_gn_lt.append(str(tm_rd_lt[6]))
                        homo_gn_md[str(tm_rd_lt[6])] = 2
                        tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
            tm_gn_file.close()
            tm_gp_file.close()
        return homo_gn_md, tm_gp_lt

    def pm_ad_wgd_pr(self, qy_table, pm_gp_file, pm_gn_file, out_prefix, pm_threshold, homo_gn_md):
        pm_gn_lt = []
        pm_gp_lt = []
        pm_pre_df = []
        pm_gp_file = open(pm_gp_file, 'w')
        pm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        pm_gn_file = open(pm_gn_file, 'w')
        pm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        # with open(qy_table) as qy_table:
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            pm_rd_lt = line.split('\t')
            pm_condition1 = pm_rd_lt[0] != pm_rd_lt[6]
            pm_condition2 = pm_rd_lt[1] == pm_rd_lt[7]
            pm_condition3 = 1 < abs(int(pm_rd_lt[2]) - int(pm_rd_lt[8])) <= pm_threshold
            if pm_condition3 and pm_condition2 and pm_condition1:
                if int(pm_rd_lt[2]) - int(pm_rd_lt[8]) < 0:
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and homo_gn_md[str(pm_rd_lt[6])] not in [1, 2]:
                        pm_pre_df.append([str(pm_rd_lt[0]), str(pm_rd_lt[6]), int(pm_rd_lt[8]) - int(pm_rd_lt[2])])
                else:
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and homo_gn_md[str(pm_rd_lt[6])] not in [1, 2]:
                        pm_pre_df.append([str(pm_rd_lt[6]), str(pm_rd_lt[0]), int(pm_rd_lt[2]) - int(pm_rd_lt[8])])
        qy_tb1.close()
        df = pd.DataFrame(pm_pre_df, columns=['Duplicate1', 'Duplicate2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        pm_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate1']) + "\t" + str(row['Duplicate2']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                pm_rd_lt = line.split("\t")
                if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) in pm_gn_pr and str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in pm_gp_lt:
                    pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                    pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                    pm_gp_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3])
                                     + "\t" + str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                     + str(pm_rd_lt[9]) + "\n")
                    if str(pm_rd_lt[0]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[0]))
                        homo_gn_md[str(pm_rd_lt[0])] = 3
                        pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    if str(pm_rd_lt[6]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[6]))
                        homo_gn_md[str(pm_rd_lt[6])] = 3
                        pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) in pm_gn_pr and str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in pm_gp_lt:
                    pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                    pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                    pm_gp_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                     + str(pm_rd_lt[9]) + "\t" + str(pm_rd_lt[0]) + "\t" + str(out_prefix)
                                     + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3]) + "\n")
                    if str(pm_rd_lt[0]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[0]))
                        homo_gn_md[str(pm_rd_lt[0])] = 3
                        pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    if str(pm_rd_lt[6]) not in pm_gn_lt:
                        pm_gn_lt.append(str(pm_rd_lt[6]))
                        homo_gn_md[str(pm_rd_lt[6])] = 3
                        pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
            pm_gn_file.close()
            pm_gp_file.close()
            return homo_gn_md, pm_gp_lt

    def trp_ad_wgd_pr(self, wgd_gn_lt, qy_rf_col, qy_table, trp_gp_file, trp_gn_file, out_prefix, seg_anc, pm_threshold, homo_gn_md):
        trp_gn_lt = []
        trp_gp_lt = []
        pre_df = []
        df = pd.read_csv(qy_table, header=None, index_col=False, sep='\t')
        df = df.iloc[:, [0, 6, 12]]
        mask = df.iloc[:, 0] != df.iloc[:, 1]
        df = df[mask]
        idt_dict = {str(row[0]) + "\t" + str(row[6]): row[12] for _, row in df.iterrows()}
        # print(idt_dict)
        # max_values = df.groupby(0)[12].transform('max')
        # print(max_values)
        # df_filtered = df[max_values == df[12]]
        # new_list = df.apply(lambda row: str(row[0]) + "\t" + str(row[6]), axis=1).tolist()
        anc = []
        qy_rf_lt = qy_rf_col.split(',')
        for file in qy_rf_lt:
            file = file.strip()
            with open(file) as fi:
                _ = next(fi)
                _ = next(fi)
                for line in fi:
                    if line.startswith("#"):
                        continue
                    col_rf_rd_lt = line.split('\t')
                    anc.append(col_rf_rd_lt[5])
        if seg_anc == 1:
            anc = anc + wgd_gn_lt
            anc = list(set(anc))
        trp_gp_file = open(trp_gp_file, 'w')
        trp_gp_file.write("Transposed" + "\t" + "Location1" + "\t" + "Parental" + "\t" + "Location2" + "\n")
        trp_gn_file = open(trp_gn_file, 'w')
        trp_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            qy_rd_lt = line.split("\t")
            if qy_rd_lt[1] != qy_rd_lt[7]:
                if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [2, 3] and homo_gn_md[qy_rd_lt[6]] not in [2, 3]:
                    pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [2, 3] and homo_gn_md[qy_rd_lt[0]] not in [2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
            else:
                if abs(int(qy_rd_lt[2]) - int(qy_rd_lt[8])) > pm_threshold:
                    if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [2, 3] and homo_gn_md[qy_rd_lt[6]] not in [2, 3]:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                    if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [2, 3] and homo_gn_md[qy_rd_lt[0]] not in [2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
        qy_tb1.close()
        df = pd.DataFrame(pre_df, columns=['Duplicate', 'Parental', 'identity'])
        df = df.drop_duplicates()
        max_values = df.groupby('Duplicate')['identity'].transform('max')
        df_filtered = df[max_values == df['identity']]
        trp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate']) + "\t" + str(row['Parental']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                qy_rd_lt = line.split("\t")
                if str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) in trp_gn_pr and str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) not in trp_gp_lt:
                    trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                    trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                    trp_gp_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":"
                                      + str(qy_rd_lt[3]) + "\t" + str(qy_rd_lt[6]) + "\t" + str(out_prefix)
                                      + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                    if homo_gn_md[qy_rd_lt[6]] == 0:
                        homo_gn_md[qy_rd_lt[6]] = 4
                    if str(qy_rd_lt[0]) not in trp_gn_lt:
                        trp_gn_lt.append(str(qy_rd_lt[0]))
                        homo_gn_md[str(qy_rd_lt[0])] = 4
                        trp_gn_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in trp_gn_pr and str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) not in trp_gp_lt:
                    trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                    trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                    trp_gp_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":"
                                      + str(qy_rd_lt[9]) + "\t" + str(qy_rd_lt[0]) + "\t" + str(out_prefix)
                                      + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                    if homo_gn_md[qy_rd_lt[0]] == 0:
                        homo_gn_md[qy_rd_lt[0]] = 4
                    if str(qy_rd_lt[6]) not in trp_gn_lt:
                        trp_gn_lt.append(str(qy_rd_lt[6]))
                        homo_gn_md[str(qy_rd_lt[6])] = 4
                        trp_gn_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
        return homo_gn_md, trp_gp_lt

    def dip_ad_pr(self, qy_table, dip_gp_file, dip_gn_file, out_prefix, homo_gn_md):
        dp_gp_lt = []
        dp_gn_lt = []
        dp_pre_df = []
        dip_gp_file = open(dip_gp_file, 'w')
        dip_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dip_gn_file = open(dip_gn_file, 'w')
        dip_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        qy_tb1 = open(qy_table, 'r')
        for line in qy_tb1:
            dp_rd_lt = line.split('\t')
            if dp_rd_lt[0] != dp_rd_lt[6]:
                if homo_gn_md[str(dp_rd_lt[0])] not in [1, 2, 3, 4] and homo_gn_md[str(dp_rd_lt[6])] not in [1, 2, 3, 4]:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(dp_rd_lt[12])])
        df = pd.DataFrame(dp_pre_df, columns=['Duplicate1', 'Duplicate2', 'identity'])
        max_values = df.groupby('Duplicate2')['identity'].transform('max')
        df_filtered = df[max_values == df['identity']]
        dp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate2']) + "\t" + str(row['Duplicate1']), axis=1).tolist()
        with open(qy_table) as qy_tb2:
            for line in qy_tb2:
                dp_rd_lt = line.split('\t')
                if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in dp_gn_pr and str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) not in dp_gp_lt:
                    dp_gp_lt.append(str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]))
                    dp_gp_lt.append(str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]))
                    dip_gp_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":"
                                      + str(dp_rd_lt[9]) + "\t" + str(dp_rd_lt[0]) + "\t" + str(out_prefix)
                                      + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                    if str(dp_rd_lt[0]) not in dp_gn_lt:
                        dp_gn_lt.append(str(dp_rd_lt[0]))
                        dip_gn_file.write(str(dp_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                    if str(dp_rd_lt[6]) not in dp_gn_lt:
                        dp_gn_lt.append(str(dp_rd_lt[6]))
                        dip_gn_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":" + str(dp_rd_lt[9]) + "\n")
        dip_gp_file.close()
        dip_gn_file.close()
        return dp_gn_lt

    def merge_mode(self, homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md

    def singleton(self, gff, homo_gn_md, out_file, out_prefix):
        chr_gn_dict, chr_gn_list, gn_chr_dict, _ = GffFile.readGff(gff)
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        for ch in chr_gn_list:
            for gn in chr_gn_list[ch]:
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + str(out_prefix) + "-" + str(ch) + ":" + str(chr_gn_dict[ch][gn.name].start) + "\n")
        ot_file.close()

    def sort_file(self, gene_file, pair_file):
        for gn_file in gene_file:
            df = pd.read_csv(gn_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(gn_file, header=True, index=False, sep='\t')
        for pr_file in pair_file:
            df = pd.read_csv(pr_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(pr_file, header=True, index=False, sep='\t')

    def run(self):
        self.file_jg(self.qy_table, self.gff, self.qy_col_fl, self.qy_rf_col, self.out_dir)
        wgd_gn_lt, wgd_pr_lt, wgd_md = self.wgd_pair(self.qy_col_fl, self.wgd_gp_file, self.wgd_gn_file, self.out_prefix)
        homo_gn_md, tm_gp_lt = self.tm_ad_wgd_pr(self.qy_table, self.tm_gp_file, self.tm_gn_file, self.out_prefix, wgd_gn_lt)
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, pm_gp_lt = self.pm_ad_wgd_pr(self.qy_table, self.pm_gp_file, self.pm_gn_file, self.out_prefix,
                                                 self.pm_threshold, homo_gn_md)
        homo_gn_md, trp_gp_lt = self.trp_ad_wgd_pr(wgd_gn_lt, self.qy_rf_col, self.qy_table, self.trp_gp_file, self.trp_gn_file,
                                                   self.out_prefix, self.seg_anc, self.pm_threshold, homo_gn_md)
        self.dip_ad_pr(self.qy_table, self.dp_gp_file, self.dp_gn_file, self.out_prefix, homo_gn_md)
        self.singleton(self.gff, homo_gn_md, self.si_gn_file, self.out_prefix)
        self.sort_file([self.wgd_gn_file, self.tm_gn_file, self.pm_gn_file, self.trp_gn_file, self.dp_gn_file],
                       [self.wgd_gp_file, self.tm_gp_file, self.pm_gp_file, self.trp_gp_file, self.dp_gp_file])
