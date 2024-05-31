import os
import pandas as pd
from . import GffFile


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


def read_blast_gff(blast, gff):
    chromosome_gene_dict, chromosome_gene_list, gene_name_chr_dict, _ = GffFile.readGff(gff)
    gene_index = dict()
    for ch in chromosome_gene_list:
        i = 0
        while i < len(chromosome_gene_list[ch]):
            gene_index[chromosome_gene_list[ch][i].name] = i + 1
            i += 1
    data_info = []
    match_pairs = set()
    with open(blast) as f:
        for line in f:
            elements = line.split()
            qseqid = elements[0]
            sseqid = elements[1]
            if qseqid == sseqid:
                continue
            pident = str(float(elements[2]))
            # length = int(elements[3])
            # mismatch = elements[4]
            # gapopen = elements[5]
            # qstart = elements[6]
            # qend = elements[7]
            # sstart = elements[8]
            # send = elements[9]
            # evalue = elements[10]
            bitscore = float(elements[11])
            match_pair = sseqid + "_" + qseqid
            if match_pair not in match_pairs:
                match_pairs.add(match_pair)
                data_info.append([sseqid, gene_name_chr_dict[sseqid],  str(gene_index[sseqid]),
                                  str(chromosome_gene_dict[gene_name_chr_dict[sseqid]][sseqid].start),
                                  str(chromosome_gene_dict[gene_name_chr_dict[sseqid]][sseqid].end),
                                  chromosome_gene_dict[gene_name_chr_dict[sseqid]][sseqid].strand,
                                  qseqid, gene_name_chr_dict[qseqid], str(gene_index[qseqid]),
                                  str(chromosome_gene_dict[gene_name_chr_dict[qseqid]][qseqid].start),
                                  str(chromosome_gene_dict[gene_name_chr_dict[qseqid]][qseqid].end),
                                  chromosome_gene_dict[gene_name_chr_dict[qseqid]][qseqid].strand, bitscore, pident])
    return data_info


class ClassGene:
    def __init__(self, config_pra):
        self.qy_blast = str(config_pra['classification']['query_blast'])
        self.gff = str(config_pra['classification']['query_gff_file'])
        self.qy_col_fl = str(config_pra["classification"]["query_query_collinearity"])
        self.qy_rf_col = str(config_pra["classification"]["query_ref_collinearity"])
        self.out_dir = str(config_pra["classification"]["out_directory"])
        self.out_prefix = str(config_pra["classification"]["out_prefix"])
        self.pm_threshold = int(config_pra["classification"]["proximal_max_distance"])
        self.seg_anc = int(config_pra["classification"]["seg_anc"])
        # self.type = int(config_pra["classification"]["type"])

        self.wgd_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}.wgd.pairs")
        self.wgd_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.wgd.genes")
        self.tm_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}.tandem.pairs")
        self.tm_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.tandem.genes")
        self.pm_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}.proximal.pairs")
        self.pm_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.proximal.genes")
        self.trp_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}.transposed.pairs")
        self.trp_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.transposed.genes")
        self.dp_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}.dispersed.pairs")
        self.dp_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.dispersed.genes")
        self.si_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}.singleton.genes")
        self.stats_file = os.path.join(self.out_dir, f"{self.out_prefix}.stats")

    @staticmethod
    def file_jg(qy_table, gff, qy_col_fl, qy_rf_col, out_dir):
        try:
            file_exist(qy_table)
            file_exist(gff)
            file_exist(qy_col_fl)
            qy_rf_lt = qy_rf_col.split(',')
            for file in qy_rf_lt:
                file = file.strip()
                file_exist(file)
        except FileNotFoundError as e:
            print("file not found:", e)
            exit(1)
        dir_exist(out_dir)

    @staticmethod
    def wgd_pair(qy_col_fl, wgd_gp_file, wgd_gn_file, out_prefix):
        wgd_gn_lt = []
        wgd_pr_lt = []
        wgd_gp_number = 0
        wgd_gn_number = 0
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
                if str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]) not in wgd_pr_lt:
                    wgd_gp_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3])
                                      + "\t" + str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8]) + "\n")
                    wgd_pr_lt.append(str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]))
                    wgd_pr_lt.append(str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]))
                    wgd_gp_number += 1
                if str(col_rd_lt[0]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[0]))
                    wgd_md[str(col_rd_lt[0])] = 1
                    wgd_gn_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3] + "\n"))
                    wgd_gn_number += 1
                if str(col_rd_lt[5]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[5]))
                    wgd_md[str(col_rd_lt[5])] = 1
                    wgd_gn_file.write(str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8] + "\n"))
                    wgd_gn_number += 1
        wgd_gn_file.close()
        wgd_gp_file.close()
        return wgd_gn_lt, wgd_pr_lt, wgd_md, wgd_gp_number, wgd_gn_number

    @staticmethod
    def tm_ad_wgd_pr(data, tm_gp_file, tm_gn_file, out_prefix, wgd_gn_lt, wgd_pr_lt):
        tm_gn_lt = []
        tm_gp_lt = []
        tm_gp_number = 0
        tm_gn_number = 0
        homo_gn_md = {}
        tm_gp_file = open(tm_gp_file, 'w')
        tm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tm_gn_file = open(tm_gn_file, 'w')
        tm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        for tm_rd_lt in data:
            if str(tm_rd_lt[0]) not in homo_gn_md:
                homo_gn_md[str(tm_rd_lt[0])] = 0
            if str(tm_rd_lt[6]) not in homo_gn_md:
                homo_gn_md[str(tm_rd_lt[6])] = 0
            tandem_condition1 = tm_rd_lt[1] == tm_rd_lt[7]
            tandem_condition2 = abs(int(tm_rd_lt[2]) - int(tm_rd_lt[8])) == 1
            if tandem_condition2 and tandem_condition1:
                # tm_gp_lt -> one
                if int(tm_rd_lt[2]) - int(tm_rd_lt[8]) < 0:
                    if str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in wgd_pr_lt and str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in tm_gp_lt:
                        tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                        tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                        tm_gp_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3])
                                         + "\t" + str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                         + str(tm_rd_lt[9]) + "\n")
                        tm_gp_number += 1
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[0]))
                            homo_gn_md[str(tm_rd_lt[0])] = 2
                            tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                            tm_gn_number += 1
                        if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[6]))
                            homo_gn_md[str(tm_rd_lt[6])] = 2
                            tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
                            tm_gn_number += 1
                else:
                    if str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in wgd_pr_lt and str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in tm_gp_lt:
                        tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                        tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                        tm_gp_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                         + str(tm_rd_lt[9]) + "\t" + str(tm_rd_lt[0]) + "\t" + str(out_prefix)
                                         + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3]) + "\n")
                        tm_gp_number += 1
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[0]))
                            homo_gn_md[str(tm_rd_lt[0])] = 2
                            tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                            tm_gn_number += 1
                        if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[6]))
                            homo_gn_md[str(tm_rd_lt[6])] = 2
                            tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
                            tm_gn_number += 1
        tm_gn_file.close()
        tm_gp_file.close()
        return homo_gn_md, tm_gp_lt, tm_gp_number, tm_gn_number

    @staticmethod
    def pm_ad_wgd_pr(data, pm_gp_file, pm_gn_file, out_prefix, wgd_pr_lt, pm_threshold, homo_gn_md):
        pm_gn_lt = []
        pm_gp_lt = []
        pm_gp_number = 0
        pm_gn_number = 0
        pm_pre_df = []
        pm_gp_file = open(pm_gp_file, 'w')
        pm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        pm_gn_file = open(pm_gn_file, 'w')
        pm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        # with open(qy_table) as qy_table:
        for pm_rd_lt in data:
            pm_condition1 = pm_rd_lt[1] == pm_rd_lt[7]
            pm_condition2 = 1 < abs(int(pm_rd_lt[2]) - int(pm_rd_lt[8])) <= pm_threshold
            if pm_condition2 and pm_condition1:
                if int(pm_rd_lt[2]) - int(pm_rd_lt[8]) < 0:
                    if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in wgd_pr_lt:
                        pm_pre_df.append([str(pm_rd_lt[0]), str(pm_rd_lt[6]), int(pm_rd_lt[8]) - int(pm_rd_lt[2])])
                else:
                    if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in wgd_pr_lt:
                        pm_pre_df.append([str(pm_rd_lt[6]), str(pm_rd_lt[0]), int(pm_rd_lt[2]) - int(pm_rd_lt[8])])
        df = pd.DataFrame(pm_pre_df, columns=['Duplicate1', 'Duplicate2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        pm_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate1']) + "\t" + str(row['Duplicate2']), axis=1).tolist()
        for pm_rd_lt in data:
            if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) in pm_gn_pr and str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in pm_gp_lt:
                pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                pm_gp_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3])
                                 + "\t" + str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                 + str(pm_rd_lt[9]) + "\n")
                pm_gp_number += 1
                if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and str(pm_rd_lt[0]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[0]))
                    homo_gn_md[str(pm_rd_lt[0])] = 3
                    pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    pm_gn_number += 1
                if homo_gn_md[str(pm_rd_lt[6])] not in [1, 2] and str(pm_rd_lt[6]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[6]))
                    homo_gn_md[str(pm_rd_lt[6])] = 3
                    pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                    pm_gn_number += 1
            if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) in pm_gn_pr and str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in pm_gp_lt:
                pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                pm_gp_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                 + str(pm_rd_lt[9]) + "\t" + str(pm_rd_lt[0]) + "\t" + str(out_prefix)
                                 + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3]) + "\n")
                pm_gp_number += 1
                if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and str(pm_rd_lt[0]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[0]))
                    homo_gn_md[str(pm_rd_lt[0])] = 3
                    pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    pm_gn_number += 1
                if homo_gn_md[str(pm_rd_lt[6])] not in [1, 2] and str(pm_rd_lt[6]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[6]))
                    homo_gn_md[str(pm_rd_lt[6])] = 3
                    pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                    pm_gn_number += 1
        pm_gn_file.close()
        pm_gp_file.close()
        return homo_gn_md, pm_gp_lt, pm_gp_number, pm_gn_number

    @staticmethod
    def trp_ad_wgd_pr(wgd_gn_lt, qy_rf_col, data, trp_gp_file, trp_gn_file, out_prefix, seg_anc, pm_threshold, homo_gn_md):
        trp_gn_lt = []
        trp_gp_lt = []
        trp_gp_number = 0
        trp_gn_number = 0
        pre_df = []
        df = pd.DataFrame(data)
        df = df.iloc[:, [0, 6, 12]]
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
        for qy_rd_lt in data:
            if qy_rd_lt[1] != qy_rd_lt[7]:
                if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                        avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), avg_bs])
                    elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                        avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), avg_bs])
                    elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
            else:
                if abs(int(qy_rd_lt[2]) - int(qy_rd_lt[8])) > pm_threshold:
                    if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                            avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(
                                idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), avg_bs])
                        elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                    if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                            avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(
                                idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), avg_bs])
                        elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
        df = pd.DataFrame(pre_df, columns=['Duplicate', 'Parental', 'bitscore'])
        df = df.drop_duplicates()
        max_values = df.groupby('Duplicate')['bitscore'].transform('max')
        df_filtered = df[max_values == df['bitscore']]
        trp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate']) + "\t" + str(row['Parental']), axis=1).tolist()
        for qy_rd_lt in data:
            if str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) in trp_gn_pr and str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) not in trp_gp_lt:
                trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                trp_gp_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":"
                                  + str(qy_rd_lt[3]) + "\t" + str(qy_rd_lt[6]) + "\t" + str(out_prefix)
                                  + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                trp_gp_number += 1
                if homo_gn_md[qy_rd_lt[6]] == 0:
                    homo_gn_md[qy_rd_lt[6]] = 4
                if str(qy_rd_lt[0]) not in trp_gn_lt:
                    trp_gn_lt.append(str(qy_rd_lt[0]))
                    homo_gn_md[str(qy_rd_lt[0])] = 4
                    trp_gn_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                    trp_gn_number += 1
            if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in trp_gn_pr and str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) not in trp_gp_lt:
                trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                trp_gp_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":"
                                  + str(qy_rd_lt[9]) + "\t" + str(qy_rd_lt[0]) + "\t" + str(out_prefix)
                                  + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                trp_gp_number += 1
                if homo_gn_md[qy_rd_lt[0]] == 0:
                    homo_gn_md[qy_rd_lt[0]] = 4
                if str(qy_rd_lt[6]) not in trp_gn_lt:
                    trp_gn_lt.append(str(qy_rd_lt[6]))
                    homo_gn_md[str(qy_rd_lt[6])] = 4
                    trp_gn_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                    trp_gn_number += 1
        trp_gn_file.close()
        trp_gp_file.close()
        return homo_gn_md, trp_gp_lt, trp_gp_number, trp_gn_number, idt_dict

    @staticmethod
    def dip_ad_pr(data, dip_gp_file, dip_gn_file, out_prefix, homo_gn_md, used_gp_lt, idt_dict):
        dp_gp_lt = []
        dp_gn_lt = []
        dp_gp_number = 0
        dp_gn_number = 0
        dp_pre_df = []
        dip_gp_file = open(dip_gp_file, 'w')
        dip_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dip_gn_file = open(dip_gn_file, 'w')
        dip_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        for dp_rd_lt in data:
            if str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]) not in used_gp_lt:
                if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in idt_dict and str(dp_rd_lt[0]) + "\t" + str(dp_rd_lt[6]) in idt_dict:
                    avg_bs = (float(idt_dict[str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0])]) + float(idt_dict[str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6])])) / 2
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), avg_bs])
                elif str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in idt_dict:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(idt_dict[str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0])])])
                else:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(dp_rd_lt[12])])
        df = pd.DataFrame(dp_pre_df, columns=['Duplicate1', 'Duplicate2', 'bitscore'])
        max_values = df.groupby('Duplicate2')['bitscore'].transform('max')
        df_filtered = df[max_values == df['bitscore']]
        dp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate2']) + "\t" + str(row['Duplicate1']), axis=1).tolist()
        for dp_rd_lt in data:
            if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in dp_gn_pr and str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) not in dp_gp_lt:
                dp_gp_lt.append(str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]))
                dp_gp_lt.append(str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]))
                dip_gp_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":"
                                  + str(dp_rd_lt[9]) + "\t" + str(dp_rd_lt[0]) + "\t" + str(out_prefix)
                                  + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                dp_gp_number += 1
                if homo_gn_md[str(dp_rd_lt[0])] not in [1, 2, 3, 4] and str(dp_rd_lt[0]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[0]))
                    dip_gn_file.write(str(dp_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                    dp_gn_number += 1
                if homo_gn_md[str(dp_rd_lt[6])] not in [1, 2, 3, 4] and str(dp_rd_lt[6]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[6]))
                    dip_gn_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":" + str(dp_rd_lt[9]) + "\n")
                    dp_gn_number += 1
        dip_gp_file.close()
        dip_gn_file.close()
        return dp_gn_lt, dp_gp_number, dp_gn_number

    @staticmethod
    def merge_mode(homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md

    @staticmethod
    def singleton(gff, homo_gn_md, out_file, out_prefix):
        chr_gn_dict, chr_gn_list, gn_chr_dict, _ = GffFile.readGff(gff)
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        single_number = 0
        index = 0
        for ch in chr_gn_list:
            for gn in chr_gn_list[ch]:
                index += 1
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + str(out_prefix) + "-" + str(ch) + ":" + str(chr_gn_dict[ch][gn.name].start) + "\n")
                    single_number += 1
        ot_file.close()
        return single_number

    @staticmethod
    def sort_file(gene_file, pair_file):
        for gn_file in gene_file:
            df = pd.read_csv(gn_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(gn_file, header=True, index=False, sep='\t')
        for pr_file in pair_file:
            df = pd.read_csv(pr_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(pr_file, header=True, index=False, sep='\t')

    def run(self):
        self.file_jg(self.qy_blast, self.gff, self.qy_col_fl, self.qy_rf_col, self.out_dir)
        wgd_gn_lt, wgd_pr_lt, wgd_md, wgd_gp_number, wgd_gn_number = self.wgd_pair(self.qy_col_fl, self.wgd_gp_file, self.wgd_gn_file, self.out_prefix)
        data = read_blast_gff(self.qy_blast, self.gff)
        homo_gn_md, tm_gp_lt, tm_gp_number, tm_gn_number = self.tm_ad_wgd_pr(data, self.tm_gp_file, self.tm_gn_file, self.out_prefix, wgd_gn_lt, wgd_pr_lt)
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, pm_gp_lt, pm_gp_number, pm_gn_number = self.pm_ad_wgd_pr(data, self.pm_gp_file, self.pm_gn_file, self.out_prefix,
                                                                             wgd_pr_lt, self.pm_threshold, homo_gn_md)
        homo_gn_md, trp_gp_lt, trp_gp_number, trp_gn_number, idt_dict = self.trp_ad_wgd_pr(wgd_gn_lt, self.qy_rf_col, data, self.trp_gp_file, self.trp_gn_file,
                                                                                           self.out_prefix, self.seg_anc, self.pm_threshold, homo_gn_md)
        used_gp_lt = wgd_pr_lt + tm_gp_lt + pm_gp_lt + trp_gp_lt
        _, dp_gp_number, dp_gn_number = self.dip_ad_pr(data, self.dp_gp_file, self.dp_gn_file, self.out_prefix, homo_gn_md, used_gp_lt, idt_dict)
        singleton_number = self.singleton(self.gff, homo_gn_md, self.si_gn_file, self.out_prefix)
        with open(self.stats_file, 'w') as f:
            f.write("Type" + "\t" + "Number" + "\n" +
                    "wgd.pairs" + "\t" + str(wgd_gp_number) + "\n" +
                    "tandem.pairs" + "\t" + str(tm_gp_number) + "\n" +
                    "proximal.pairs" + "\t" + str(pm_gp_number) + "\n" +
                    "transposed.pairs" + "\t" + str(trp_gp_number) + "\n" +
                    "dispersed.pairs" + "\t" + str(dp_gp_number) + "\n" +
                    "wgd.genes" + "\t" + str(wgd_gn_number) + "\n" +
                    "tandem.genes" + "\t" + str(tm_gn_number) + "\n" +
                    "proximal.genes" + "\t" + str(pm_gn_number) + "\n" +
                    "transposed.genes" + "\t" + str(trp_gn_number) + "\n" +
                    "dispersed.genes" + "\t" + str(dp_gn_number) + "\n" +
                    "singleton.genes" + "\t" + str(singleton_number))

        self.sort_file([self.wgd_gn_file, self.tm_gn_file, self.pm_gn_file, self.trp_gn_file, self.dp_gn_file],
                       [self.wgd_gp_file, self.tm_gp_file, self.pm_gp_file, self.trp_gp_file, self.dp_gp_file])


class ClassGeneUnique:
    def __init__(self, config_pra):
        self.qy_blast = str(config_pra['classification']['query_blast'])
        self.gff = str(config_pra['classification']['query_gff_file'])
        self.qy_col_fl = str(config_pra["classification"]["query_query_collinearity"])
        self.qy_rf_col = str(config_pra["classification"]["query_ref_collinearity"])
        self.out_dir = str(config_pra["classification"]["out_directory"])
        self.out_prefix = str(config_pra["classification"]["out_prefix"])
        self.pm_threshold = int(config_pra["classification"]["proximal_max_distance"])
        self.seg_anc = int(config_pra["classification"]["seg_anc"])
        # self.type = int(config_pra["classification"]["type"])

        self.wgd_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.wgd.pairs")
        self.wgd_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.wgd.genes")
        self.tm_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.tandem.pairs")
        self.tm_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.tandem.genes")
        self.pm_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.proximal.pairs")
        self.pm_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.proximal.genes")
        self.trp_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.transposed.pairs")
        self.trp_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.transposed.genes")
        self.dp_gp_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.dispersed.pairs")
        self.dp_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.dispersed.genes")
        self.si_gn_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.singleton.genes")
        self.stats_file = os.path.join(self.out_dir, f"{self.out_prefix}-unique.stats")

    @staticmethod
    def file_jg(qy_blast, gff, qy_col_fl, qy_rf_col, out_dir):
        try:
            file_exist(qy_blast)
            file_exist(gff)
            file_exist(qy_col_fl)
            qy_rf_lt = qy_rf_col.split(',')
            for file in qy_rf_lt:
                file = file.strip()
                file_exist(file)
        except FileNotFoundError as e:
            print("file not found:", e)
            exit(1)
        dir_exist(out_dir)

    @staticmethod
    def wgd_pair(qy_col_fl, wgd_gp_file, wgd_gn_file, out_prefix):
        wgd_gn_lt = []
        wgd_pr_lt = []
        wgd_gp_number = 0
        wgd_gn_number = 0
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
                    wgd_gp_number += 1
                    wgd_pr_lt.append(str(col_rd_lt[0]) + "\t" + str(col_rd_lt[5]))
                    wgd_pr_lt.append(str(col_rd_lt[5]) + "\t" + str(col_rd_lt[0]))
                if str(col_rd_lt[0]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[0]))
                    wgd_md[str(col_rd_lt[0])] = 1
                    wgd_gn_file.write(str(col_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[1]) + ":" + str(col_rd_lt[3] + "\n"))
                    wgd_gn_number += 1
                if str(col_rd_lt[5]) not in wgd_gn_lt:
                    wgd_gn_lt.append(str(col_rd_lt[5]))
                    wgd_md[str(col_rd_lt[5])] = 1
                    wgd_gn_file.write(str(col_rd_lt[5]) + "\t" + str(out_prefix) + "-" + str(col_rd_lt[6]) + ":" + str(col_rd_lt[8] + "\n"))
                    wgd_gn_number += 1
        wgd_gn_file.close()
        wgd_gp_file.close()
        return wgd_gn_lt, wgd_pr_lt, wgd_md, wgd_gp_number, wgd_gn_number

    @staticmethod
    def tm_ad_wgd_pr(data, tm_gp_file, tm_gn_file, out_prefix, wgd_gn_lt):
        tm_gn_lt = []
        tm_gp_lt = []
        tm_gp_number = 0
        tm_gn_number = 0
        homo_gn_md = {}
        tm_gp_file = open(tm_gp_file, 'w')
        tm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tm_gn_file = open(tm_gn_file, 'w')
        tm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        for tm_rd_lt in data:
            if tm_rd_lt[0] not in homo_gn_md:
                homo_gn_md[str(tm_rd_lt[0])] = 0
            if tm_rd_lt[6] not in homo_gn_md:
                homo_gn_md[str(tm_rd_lt[6])] = 0
            tandem_condition1 = tm_rd_lt[1] == tm_rd_lt[7]
            tandem_condition2 = abs(int(tm_rd_lt[2]) - int(tm_rd_lt[8])) == 1
            if tandem_condition2 and tandem_condition1:
                # tm_gp_lt -> one
                if int(tm_rd_lt[2]) - int(tm_rd_lt[8]) < 0:
                    if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in wgd_gn_lt and \
                           str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]) not in tm_gp_lt:
                        tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                        tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                        tm_gp_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3])
                                         + "\t" + str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                         + str(tm_rd_lt[9]) + "\n")
                        tm_gp_number += 1
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[0]))
                            homo_gn_md[str(tm_rd_lt[0])] = 2
                            tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                            tm_gn_number += 1
                        if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[6]))
                            homo_gn_md[str(tm_rd_lt[6])] = 2
                            tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
                            tm_gn_number += 1
                else:
                    if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in wgd_gn_lt and \
                            str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]) not in tm_gp_lt:
                        tm_gp_lt.append(str(tm_rd_lt[0]) + "\t" + str(tm_rd_lt[6]))
                        tm_gp_lt.append(str(tm_rd_lt[6]) + "\t" + str(tm_rd_lt[0]))
                        tm_gp_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":"
                                         + str(tm_rd_lt[9]) + "\t" + str(tm_rd_lt[0]) + "\t" + str(out_prefix)
                                         + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3]) + "\n")
                        tm_gp_number += 1
                        if str(tm_rd_lt[0]) not in wgd_gn_lt and str(tm_rd_lt[0]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[0]))
                            homo_gn_md[str(tm_rd_lt[0])] = 2
                            tm_gn_file.write(str(tm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[1]) + ":" + str(tm_rd_lt[3] + "\n"))
                            tm_gn_number += 1
                        if str(tm_rd_lt[6]) not in wgd_gn_lt and str(tm_rd_lt[6]) not in tm_gn_lt:
                            tm_gn_lt.append(str(tm_rd_lt[6]))
                            homo_gn_md[str(tm_rd_lt[6])] = 2
                            tm_gn_file.write(str(tm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(tm_rd_lt[7]) + ":" + str(tm_rd_lt[9] + "\n"))
                            tm_gn_number += 1
        tm_gn_file.close()
        tm_gp_file.close()
        return homo_gn_md, tm_gp_lt, tm_gp_number, tm_gn_number

    @staticmethod
    def pm_ad_wgd_pr(data, pm_gp_file, pm_gn_file, out_prefix, pm_threshold, homo_gn_md):
        pm_gn_lt = []
        pm_gp_lt = []
        pm_gp_number = 0
        pm_gn_number = 0
        pm_pre_df = []
        pm_gp_file = open(pm_gp_file, 'w')
        pm_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        pm_gn_file = open(pm_gn_file, 'w')
        pm_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        # with open(qy_table) as qy_table:
        for pm_rd_lt in data:
            pm_condition1 = pm_rd_lt[1] == pm_rd_lt[7]
            pm_condition2 = 1 < abs(int(pm_rd_lt[2]) - int(pm_rd_lt[8])) <= pm_threshold
            if pm_condition2 and pm_condition1:
                if int(pm_rd_lt[2]) - int(pm_rd_lt[8]) < 0:
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and homo_gn_md[str(pm_rd_lt[6])] not in [1, 2]:
                        pm_pre_df.append([str(pm_rd_lt[0]), str(pm_rd_lt[6]), int(pm_rd_lt[8]) - int(pm_rd_lt[2])])
                else:
                    if homo_gn_md[str(pm_rd_lt[0])] not in [1, 2] and homo_gn_md[str(pm_rd_lt[6])] not in [1, 2]:
                        pm_pre_df.append([str(pm_rd_lt[6]), str(pm_rd_lt[0]), int(pm_rd_lt[2]) - int(pm_rd_lt[8])])
        df = pd.DataFrame(pm_pre_df, columns=['Duplicate1', 'Duplicate2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        pm_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate1']) + "\t" + str(row['Duplicate2']), axis=1).tolist()
        for pm_rd_lt in data:
            if str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) in pm_gn_pr and str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]) not in pm_gp_lt:
                pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                pm_gp_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3])
                                 + "\t" + str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                 + str(pm_rd_lt[9]) + "\n")
                pm_gp_number += 1
                if str(pm_rd_lt[0]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[0]))
                    homo_gn_md[str(pm_rd_lt[0])] = 3
                    pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    pm_gn_number += 1
                if str(pm_rd_lt[6]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[6]))
                    homo_gn_md[str(pm_rd_lt[6])] = 3
                    pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                    pm_gn_number += 1
            if str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) in pm_gn_pr and str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]) not in pm_gp_lt:
                pm_gp_lt.append(str(pm_rd_lt[6]) + "\t" + str(pm_rd_lt[0]))
                pm_gp_lt.append(str(pm_rd_lt[0]) + "\t" + str(pm_rd_lt[6]))
                pm_gp_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":"
                                 + str(pm_rd_lt[9]) + "\t" + str(pm_rd_lt[0]) + "\t" + str(out_prefix)
                                 + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3]) + "\n")
                pm_gp_number += 1
                if str(pm_rd_lt[0]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[0]))
                    homo_gn_md[str(pm_rd_lt[0])] = 3
                    pm_gn_file.write(str(pm_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[1]) + ":" + str(pm_rd_lt[3] + "\n"))
                    pm_gn_number += 1
                if str(pm_rd_lt[6]) not in pm_gn_lt:
                    pm_gn_lt.append(str(pm_rd_lt[6]))
                    homo_gn_md[str(pm_rd_lt[6])] = 3
                    pm_gn_file.write(str(pm_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(pm_rd_lt[7]) + ":" + str(pm_rd_lt[9] + "\n"))
                    pm_gn_number += 1
        pm_gn_file.close()
        pm_gp_file.close()
        return homo_gn_md, pm_gp_lt, pm_gp_number, pm_gn_number

    @staticmethod
    def trp_ad_wgd_pr(wgd_gn_lt, qy_rf_col, data, trp_gp_file, trp_gn_file, out_prefix, seg_anc, pm_threshold, homo_gn_md):
        trp_gn_lt = []
        trp_gp_lt = []
        trp_gp_number = 0
        trp_gn_number = 0
        pre_df = []
        df = pd.DataFrame(data)
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
        for qy_rd_lt in data:
            if qy_rd_lt[1] != qy_rd_lt[7]:
                if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3] and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                        avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), avg_bs])
                    elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3] and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3]:
                    if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                        avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), avg_bs])
                    elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                    else:
                        pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
            else:
                if abs(int(qy_rd_lt[2]) - int(qy_rd_lt[8])) > pm_threshold:
                    if qy_rd_lt[0] not in anc and qy_rd_lt[6] in anc and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3] and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                            avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])]) + float(
                                idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), avg_bs])
                        elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[0]), str(qy_rd_lt[6]), float(qy_rd_lt[12])])
                    if qy_rd_lt[6] not in anc and qy_rd_lt[0] in anc and homo_gn_md[qy_rd_lt[6]] not in [1, 2, 3] and homo_gn_md[qy_rd_lt[0]] not in [1, 2, 3]:
                        if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict and str(qy_rd_lt[0]) + "\t" + str(qy_rd_lt[6]) in idt_dict:
                            avg_bs = (float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])
                                      + float(idt_dict[str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6])])) / 2
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), avg_bs])
                        elif str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in idt_dict:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(idt_dict[str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0])])])
                        else:
                            pre_df.append([str(qy_rd_lt[6]), str(qy_rd_lt[0]), float(qy_rd_lt[12])])
        df = pd.DataFrame(pre_df, columns=['Duplicate', 'Parental', 'bitscore'])
        df = df.drop_duplicates()
        max_values = df.groupby('Duplicate')['bitscore'].transform('max')
        df_filtered = df[max_values == df['bitscore']]
        trp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate']) + "\t" + str(row['Parental']), axis=1).tolist()
        for qy_rd_lt in data:
            if str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) in trp_gn_pr and str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]) not in trp_gp_lt:
                trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                trp_gp_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":"
                                  + str(qy_rd_lt[3]) + "\t" + str(qy_rd_lt[6]) + "\t" + str(out_prefix)
                                  + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                trp_gp_number += 1
                if homo_gn_md[qy_rd_lt[6]] == 0:
                    homo_gn_md[qy_rd_lt[6]] = 4
                if str(qy_rd_lt[0]) not in trp_gn_lt:
                    trp_gn_lt.append(str(qy_rd_lt[0]))
                    homo_gn_md[str(qy_rd_lt[0])] = 4
                    trp_gn_file.write(str(qy_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                    trp_gn_number += 1
            if str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) in trp_gn_pr and str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]) not in trp_gp_lt:
                trp_gp_lt.append(str(qy_rd_lt[6]) + '\t' + str(qy_rd_lt[0]))
                trp_gp_lt.append(str(qy_rd_lt[0]) + '\t' + str(qy_rd_lt[6]))
                trp_gp_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":"
                                  + str(qy_rd_lt[9]) + "\t" + str(qy_rd_lt[0]) + "\t" + str(out_prefix)
                                  + "-" + str(qy_rd_lt[1]) + ":" + str(qy_rd_lt[3]) + "\n")
                trp_gp_number += 1
                if homo_gn_md[qy_rd_lt[0]] == 0:
                    homo_gn_md[qy_rd_lt[0]] = 4
                if str(qy_rd_lt[6]) not in trp_gn_lt:
                    trp_gn_lt.append(str(qy_rd_lt[6]))
                    homo_gn_md[str(qy_rd_lt[6])] = 4
                    trp_gn_file.write(str(qy_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(qy_rd_lt[7]) + ":" + str(qy_rd_lt[9]) + "\n")
                    trp_gn_number += 1
        trp_gn_file.close()
        trp_gp_file.close()
        return homo_gn_md, trp_gp_lt, trp_gp_number, trp_gn_number, idt_dict

    @staticmethod
    def dip_ad_pr(data, dip_gp_file, dip_gn_file, out_prefix, homo_gn_md, idt_dict):
        dp_gp_lt = []
        dp_gn_lt = []
        dp_gp_number = 0
        dp_gn_number = 0
        dp_pre_df = []
        dip_gp_file = open(dip_gp_file, 'w')
        dip_gp_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dip_gn_file = open(dip_gn_file, 'w')
        dip_gn_file.write("Duplicate" + "\t" + "Location" + "\n")
        for dp_rd_lt in data:
            if homo_gn_md[str(dp_rd_lt[0])] not in [1, 2, 3, 4] and homo_gn_md[str(dp_rd_lt[6])] not in [1, 2, 3, 4]:
                if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in idt_dict and str(dp_rd_lt[0]) + "\t" + str(dp_rd_lt[6]) in idt_dict:
                    avg_bs = (float(idt_dict[str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0])]) + float(idt_dict[str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6])])) / 2
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), avg_bs])
                elif str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in idt_dict:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(idt_dict[str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0])])])
                else:
                    dp_pre_df.append([str(dp_rd_lt[0]), str(dp_rd_lt[6]), float(dp_rd_lt[12])])
        df = pd.DataFrame(dp_pre_df, columns=['Duplicate1', 'Duplicate2', 'bitscore'])
        max_values = df.groupby('Duplicate2')['bitscore'].transform('max')
        df_filtered = df[max_values == df['bitscore']]
        dp_gn_pr = df_filtered.apply(lambda row: str(row['Duplicate2']) + "\t" + str(row['Duplicate1']), axis=1).tolist()
        for dp_rd_lt in data:
            if str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) in dp_gn_pr and str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]) not in dp_gp_lt:
                dp_gp_lt.append(str(dp_rd_lt[0]) + '\t' + str(dp_rd_lt[6]))
                dp_gp_lt.append(str(dp_rd_lt[6]) + '\t' + str(dp_rd_lt[0]))
                dip_gp_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":"
                                  + str(dp_rd_lt[9]) + "\t" + str(dp_rd_lt[0]) + "\t" + str(out_prefix)
                                  + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                dp_gp_number += 1
                if str(dp_rd_lt[0]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[0]))
                    dip_gn_file.write(str(dp_rd_lt[0]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[1]) + ":" + str(dp_rd_lt[3]) + "\n")
                    dp_gn_number += 1
                if str(dp_rd_lt[6]) not in dp_gn_lt:
                    dp_gn_lt.append(str(dp_rd_lt[6]))
                    dip_gn_file.write(str(dp_rd_lt[6]) + "\t" + str(out_prefix) + "-" + str(dp_rd_lt[7]) + ":" + str(dp_rd_lt[9]) + "\n")
                    dp_gn_number += 1
        dip_gp_file.close()
        dip_gn_file.close()
        return dp_gn_lt, dp_gp_number, dp_gn_number

    @staticmethod
    def merge_mode(homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md

    @staticmethod
    def singleton(gff, homo_gn_md, out_file, out_prefix):
        chr_gn_dict, chr_gn_list, gn_chr_dict, _ = GffFile.readGff(gff)
        single_number = 0
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        for ch in chr_gn_list:
            for gn in chr_gn_list[ch]:
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + str(out_prefix) + "-" + str(ch) + ":" + str(chr_gn_dict[ch][gn.name].start) + "\n")
                    single_number += 1
        ot_file.close()
        return single_number

    @staticmethod
    def sort_file(gene_file, pair_file):
        for gn_file in gene_file:
            df = pd.read_csv(gn_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(gn_file, header=True, index=False, sep='\t')
        for pr_file in pair_file:
            df = pd.read_csv(pr_file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(pr_file, header=True, index=False, sep='\t')

    def run(self):
        self.file_jg(self.qy_blast, self.gff, self.qy_col_fl, self.qy_rf_col, self.out_dir)
        wgd_gn_lt, wgd_pr_lt, wgd_md, wgd_gp_number, wgd_gn_number = self.wgd_pair(self.qy_col_fl, self.wgd_gp_file, self.wgd_gn_file, self.out_prefix)
        data = read_blast_gff(self.qy_blast, self.gff)

        homo_gn_md, tm_gp_lt, tm_gp_number, tm_gn_number = self.tm_ad_wgd_pr(data, self.tm_gp_file, self.tm_gn_file, self.out_prefix, wgd_gn_lt)
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, pm_gp_lt, pm_gp_number, pm_gn_number = self.pm_ad_wgd_pr(data, self.pm_gp_file, self.pm_gn_file, self.out_prefix,
                                                                             self.pm_threshold, homo_gn_md)
        homo_gn_md, trp_gp_lt, trp_gp_number, trp_gn_number, idt_dict = self.trp_ad_wgd_pr(wgd_gn_lt, self.qy_rf_col, data, self.trp_gp_file, self.trp_gn_file,
                                                                                           self.out_prefix, self.seg_anc, self.pm_threshold, homo_gn_md)
        dp_gn_lt, dp_gp_number, dp_gn_number = self.dip_ad_pr(data, self.dp_gp_file, self.dp_gn_file, self.out_prefix, homo_gn_md, idt_dict)
        singleton_number = self.singleton(self.gff, homo_gn_md, self.si_gn_file, self.out_prefix)
        with open(self.stats_file, 'w') as f:
            f.write("Type" + "\t" + "Number" + "\n" +
                    "wgd.pairs" + "\t" + str(wgd_gp_number) + "\n" +
                    "tandem.pairs" + "\t" + str(tm_gp_number) + "\n" +
                    "proximal.pairs" + "\t" + str(pm_gp_number) + "\n" +
                    "transposed.pairs" + "\t" + str(trp_gp_number) + "\n" +
                    "dispersed.pairs" + "\t" + str(dp_gp_number) + "\n" +
                    "wgd.genes" + "\t" + str(wgd_gn_number) + "\n" +
                    "tandem.genes" + "\t" + str(tm_gn_number) + "\n" +
                    "proximal.genes" + "\t" + str(pm_gn_number) + "\n" +
                    "transposed.genes" + "\t" + str(trp_gn_number) + "\n" +
                    "dispersed.genes" + "\t" + str(dp_gn_number) + "\n" +
                    "singleton.genes" + "\t" + str(singleton_number))
        self.sort_file([self.wgd_gn_file, self.tm_gn_file, self.pm_gn_file, self.trp_gn_file, self.dp_gn_file],
                       [self.wgd_gp_file, self.tm_gp_file, self.pm_gp_file, self.trp_gp_file, self.dp_gp_file])
