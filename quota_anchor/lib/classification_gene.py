# This script is refers to https://github.com/qiao-xin/DupGen_finder.
# We use our syntenic algorithm rather than MCScanX to identify wgd gene pairs/wgd genes.
# When you specify type=0, duplicates are permitted between different classes of genes to some extent.
import os
import time
from . import base
import pandas as pd
from . import GffFile


def file_judge(query_blast, query_gff, query_query_collinearity, query_ref_collinearity, out_directory):
    base.file_empty(query_blast)
    base.file_empty(query_gff)
    base.file_empty(query_query_collinearity)
    query_ref_collinearity_list = query_ref_collinearity.split(',')
    for file in query_ref_collinearity_list:
        file = file.strip()
        if len(file) == 0:
            continue
        base.file_empty(file)
    os.makedirs(out_directory, exist_ok=True)

def overwrite_judge(files, overwrite): 
    for file in files:
        base.output_file_parentdir_exist(file, overwrite)

def split_data_line(line):
    ref_gene = line[0]
    ref_chr = line[1]
    ref_order = int(line[2])
    ref_start = line[3]

    query_gene = line[6]
    query_chr = line[7]
    query_order = int(line[8])
    query_start = line[9]
    return ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start

def split_collinearity_line(line):
    list_element = line.split('\t')
    ref_gene = list_element[0]
    ref_chr = list_element[1]
    ref_start = list_element[3]
    
    query_gene = list_element[5]
    query_chr = list_element[6]
    query_start = list_element[8]
    return ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start

def get_anc(query_ref_collinearity, seg_anc, wgd_gene_list, anc):
    query_ref_collinearity_list = query_ref_collinearity.split(',')
    for file in query_ref_collinearity_list:
        file = file.strip()
        if len(file) == 0:
            continue
        with open(file) as fi:
            _ = next(fi)
            _ = next(fi)
            for line in fi:
                if line.startswith("#"):
                    continue
                col_rf_rd_lt = line.split('\t')
                anc[col_rf_rd_lt[5]] = True
                anc[col_rf_rd_lt[0]] = True
    if seg_anc != 0:
        for gene in wgd_gene_list:
            anc[gene] = True
    return anc

def write_gene_to_wgd_transposed_file(gene_name, chr, start, wgd_md, file_handler, gene_list, output_prefix, gene_number, mode):
    if gene_name not in gene_list:
        gene_list.append(gene_name)
        wgd_md[gene_name] = mode
        file_handler.write(gene_name + "\t" + output_prefix + "-" + chr + ":" + start + "\n")
        gene_number += 1
    return gene_list, gene_number, wgd_md

def write_gene_to_file(gene_name, chr, start, homo_gn_md, file_handler, gene_list, output_prefix, gene_number, mode):
    if gene_name not in gene_list:
        gene_list.append(gene_name)
        if homo_gn_md[gene_name] == 0:
            homo_gn_md[gene_name] = mode
        file_handler.write(gene_name + "\t" + output_prefix + "-" + chr + ":" + start + "\n")
        gene_number += 1
    return gene_list, gene_number, homo_gn_md

def init_homo_gn_md(ref_gene, ref_chr, ref_order, query_gene, query_chr, query_order, homo_gn_md, anc, used_gp_dict):
    if ref_gene not in homo_gn_md:
        homo_gn_md[ref_gene] = 0
        anc[ref_gene] = False
    if query_gene not in homo_gn_md:
        homo_gn_md[query_gene] = 0
        anc[query_gene] = False
    used_gp_dict[ref_gene + '\t' + query_gene] = False
    used_gp_dict[query_gene + '\t' + ref_gene] = False
    tandem_condition1 = ref_chr == query_chr
    tandem_condition2 = abs(ref_order - query_order) == 1
    condition = tandem_condition1 and tandem_condition2
    return condition, homo_gn_md, anc, used_gp_dict

def transposed_pre_df(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line_bitscore, line_identity, anc, homo_gn_md):
    if not anc[ref_gene] and anc[query_gene] and homo_gn_md[ref_gene] not in [1, 2, 3]:
        try:
            avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(
                bitscore_dict[ref_gene + '\t' + query_gene])) / 2
            avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(
                identity_dict[ref_gene + '\t' + query_gene])) / 2
            pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, avg_bs, avg_id])
        except KeyError:
            if query_gene + '\t' + ref_gene in bitscore_dict:
                pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(bitscore_dict[query_gene + '\t' + ref_gene]), float(identity_dict[query_gene + '\t' + ref_gene])])
            else: 
                pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(line_bitscore), float(line_identity)])
        return pre_df
    if not anc[query_gene] and anc[ref_gene] and homo_gn_md[query_gene] not in [1, 2, 3]:
        try:
            avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(
                bitscore_dict[ref_gene + '\t' + query_gene])) / 2
            avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(
                identity_dict[ref_gene + '\t' + query_gene])) / 2
            pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, avg_bs, avg_id])
        except KeyError:
            if query_gene + '\t' + ref_gene in bitscore_dict:
                pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, float(bitscore_dict[query_gene + '\t' + ref_gene]), float(identity_dict[query_gene + '\t' + ref_gene])])
            else:
                pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, float(line_bitscore), float(line_identity)])
    return pre_df

def transposed_pre_df_unique(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line_bitscore, line_identity, anc, homo_gn_md):
    if not anc[ref_gene] and anc[query_gene] and homo_gn_md[ref_gene] not in [1, 2, 3] and homo_gn_md[query_gene] not in [1, 2, 3]:
        try:
            avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(
                bitscore_dict[ref_gene + '\t' + query_gene])) / 2
            avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(
                identity_dict[ref_gene + '\t' + query_gene])) / 2
            pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, avg_bs, avg_id])
        except KeyError:
            if query_gene + '\t' + ref_gene in bitscore_dict:
                pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(bitscore_dict[query_gene + '\t' + ref_gene]), float(identity_dict[query_gene + '\t' + ref_gene])])
            else: 
                pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(line_bitscore), float(line_identity)])
        return pre_df
    if not anc[query_gene] and anc[ref_gene] and homo_gn_md[query_gene] not in [1, 2, 3] and homo_gn_md[ref_gene] not in [1, 2, 3]:
        try:
            avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(
                bitscore_dict[ref_gene + '\t' + query_gene])) / 2
            avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(
                identity_dict[ref_gene + '\t' + query_gene])) / 2
            pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, avg_bs, avg_id])
        except KeyError:
            if query_gene + '\t' + ref_gene in bitscore_dict:
                pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, float(bitscore_dict[query_gene + '\t' + ref_gene]), float(identity_dict[query_gene + '\t' + ref_gene])])
            else:
                pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, float(line_bitscore), float(line_identity)])
    return pre_df

def get_used_dict(wgd_pair_list, tandem_pair_lt, proximal_pair_list, transposed_pair_list, used_gp_dict):
    used_gp_lt = wgd_pair_list + tandem_pair_lt + proximal_pair_list + transposed_pair_list
    for pair in used_gp_lt:
        used_gp_dict[pair] = True
    return used_gp_dict

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
            pident = float(elements[2])
            bitscore = float(elements[11])
            match_pair = sseqid + "_" + qseqid
            # get the first match
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
    def __init__(self, config_pra, parameter):
        for i in config_pra.sections():
            if i == 'classification':
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
        self.seg_anc = int(self.seg_anc)
        self.proximal_max_distance = int(self.proximal_max_distance)

        self.wgd_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}.wgd.pairs")
        self.wgd_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.wgd.genes")
        self.tandem_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}.tandem.pairs")
        self.tandem_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.tandem.genes")
        self.proximal_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}.proximal.pairs")
        self.proximal_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.proximal.genes")
        self.transposed_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}.transposed.pairs")
        self.transposed_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.transposed.genes")
        self.dispersed_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}.dispersed.pairs")
        self.dispersed_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.dispersed.genes")
        self.singleton_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}.singleton.genes")
        self.stats_file = os.path.join(self.out_directory, f"{self.output_prefix}.stats")


    @staticmethod
    def wgd_process(query_query_collinearity, wgd_pair_file, wgd_gene_file, output_prefix):
        mode = 1
        wgd_gene_list, wgd_pair_list = [], []
        wgd_gene_number, wgd_pair_number = 0, 0
        wgd_md = {}
        
        wgd_pair_file = open(wgd_pair_file, 'w')
        wgd_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        wgd_gene_file = open(wgd_gene_file, 'w')
        wgd_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        
        with open(query_query_collinearity) as wgd_col_file:
            _ = next(wgd_col_file)
            _ = next(wgd_col_file)
            # gene:SORBI_3010G080200	10	828	6837385	6840092	gene:SORBI_3002G223000	2	2281	61442049	61444165	+	0.7043
            for line in wgd_col_file:
                if line.startswith("#"):
                    continue
                # TODO add gene order to the location 
                ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start = split_collinearity_line(line)
                
                if ref_gene + "\t" + query_gene not in wgd_pair_list:
                    wgd_pair_file.write(ref_gene + "\t" + output_prefix + "-" + ref_chr + ":" + ref_start
                                      + "\t" + query_gene + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\n")
                    wgd_pair_list.append(ref_gene + "\t" + query_gene)
                    wgd_pair_list.append(query_gene + "\t" + ref_gene)
                    wgd_pair_number += 1
                
                    wgd_gene_list, wgd_gene_number, wgd_md = write_gene_to_wgd_transposed_file(ref_gene, ref_chr, ref_start, 
                                                                        wgd_md, wgd_gene_file, wgd_gene_list, output_prefix, wgd_gene_number, mode)
                    wgd_gene_list, wgd_gene_number, wgd_md = write_gene_to_wgd_transposed_file(query_gene, query_chr, query_start,
                                                                        wgd_md, wgd_gene_file, wgd_gene_list, output_prefix, wgd_gene_number, mode)
        wgd_gene_file.close()
        wgd_pair_file.close()
        return wgd_gene_list, wgd_pair_list, wgd_md, wgd_pair_number, wgd_gene_number

    @staticmethod
    def tandem_process(data, tandem_pair_file, tandem_gene_file, output_prefix):
        mode = 2
        homo_gn_md = {}
        anc = {}
        used_gp_dict = {}
        tandem_gene_list, tandem_pair_list = [], []
        tandem_gene_number, tandem_pair_number = 0, 0
        tandem_pair_file = open(tandem_pair_file, 'w')
        tandem_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tandem_gene_file = open(tandem_gene_file, 'w')
        tandem_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            tandem_condition, homo_gn_md, anc, used_gp_dict = init_homo_gn_md(ref_gene, ref_chr, ref_order, query_gene, query_chr, query_order, homo_gn_md, anc, used_gp_dict)
            
            if tandem_condition:
                if ref_order - query_order < 0:
                    if ref_gene + "\t" + query_gene not in tandem_pair_list:
                        tandem_pair_list.append(ref_gene + "\t" + query_gene)
                        tandem_pair_list.append(query_gene + "\t" + ref_gene)
                        tandem_pair_file.write(ref_gene + "\t" + output_prefix + "-" + ref_chr + ":" + ref_start
                                         + "\t" + query_gene + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\n")
                        tandem_pair_number += 1

                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(ref_gene, ref_chr, ref_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(query_gene, query_chr, query_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                                                                            
                else:
                    if query_gene + "\t" + ref_gene not in tandem_pair_list:
                        tandem_pair_list.append(ref_gene + "\t" + query_gene)
                        tandem_pair_list.append(query_gene + "\t" + ref_gene)
                        tandem_pair_file.write(query_gene + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\t" + 
                                               ref_gene + "\t" + output_prefix + "-" + ref_chr + ":" + ref_start + "\n")
                        tandem_pair_number += 1
                    
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(ref_gene, ref_chr, ref_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(query_gene, query_chr, query_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
        tandem_gene_file.close()
        tandem_pair_file.close()
        return homo_gn_md, tandem_pair_list, tandem_pair_number, tandem_gene_number, anc, used_gp_dict

    @staticmethod
    def proximal_process(data, proximal_pair_file, proximal_gene_file, output_prefix, proximal_threshold, homo_gn_md):
        mode = 3
        proximal_gene_list, proximal_pair_list = [], []
        proximal_gene_number, proximal_pair_number = 0, 0
        proximal_pre_df = []
        
        proximal_pair_file = open(proximal_pair_file, 'w')
        proximal_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        proximal_gene_file = open(proximal_gene_file, 'w')
        proximal_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            proximal_condition1 = ref_chr == query_chr
            proximal_condition2 = 1 < abs(ref_order - query_order) <= proximal_threshold
            if proximal_condition2 and proximal_condition1:
                if ref_order - query_order < 0:
                    proximal_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start ,query_order - ref_order])
                else:
                    proximal_pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, ref_order - query_order])

        df = pd.DataFrame(proximal_pre_df, columns=['Duplicate1', 'chr1', 'start1', 'Duplicate2', 'chr2', 'start2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        df_filtered.reset_index(drop=True, inplace=True)

        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, distance = list(line)
            proximal_pair_list.append(dup_gene1 + "\t" + dup_gene2)
            proximal_pair_list.append(dup_gene2 + "\t" + dup_gene1)
            proximal_pair_file.write(dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + ref_start
                                + "\t" + dup_gene2 + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\n")

            proximal_gene_list, proximal_gene_number, homo_gn_md = write_gene_to_file(dup_gene1, chr1, start1, homo_gn_md,
                                                                                        proximal_gene_file, proximal_gene_list, output_prefix, proximal_gene_number, mode)
            proximal_gene_list, proximal_gene_number, homo_gn_md = write_gene_to_file(dup_gene2, chr2, start2, homo_gn_md,
                                                                                        proximal_gene_file, proximal_gene_list, output_prefix, proximal_gene_number, mode)

        proximal_gene_file.close()
        proximal_pair_file.close()
        proximal_pair_number = index + 1

        return homo_gn_md, proximal_pair_list, proximal_pair_number, proximal_gene_number

    @staticmethod
    def transposed_process(wgd_gene_list, query_ref_collinearity, data, transposed_pair_file, transposed_gene_file, output_prefix, seg_anc, proximal_threshold, homo_gn_md, anc):
        mode = 4
        transposed_gene_list, transposed_pair_list = [], []
        transposed_gene_number, transposed_pair_number = 0, 0
        pre_df = []
        df = pd.DataFrame(data)
        df = df.iloc[:, [0, 6, 12, 13]]
        bitscore_dict = {str(row[0]) + "\t" + str(row[6]): row[12] for _, row in df.iterrows()}
        identity_dict = {str(row[0]) + "\t" + str(row[6]): row[13] for _, row in df.iterrows()}
        anc = get_anc(query_ref_collinearity, seg_anc, wgd_gene_list, anc)

        transposed_pair_file = open(transposed_pair_file, 'w')
        transposed_pair_file.write("Transposed" + "\t" + "Location1" + "\t" + "Parental" + "\t" + "Location2" + "\n")
        transposed_gene_file = open(transposed_gene_file, 'w')
        transposed_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            if ref_chr != query_chr:
                pre_df = transposed_pre_df(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line[12], line[13], anc, homo_gn_md)
            else:
                if abs(ref_order - query_order) > proximal_threshold:
                    pre_df = transposed_pre_df(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line[12], line[13], anc, homo_gn_md)

        df = pd.DataFrame(pre_df, columns=['Duplicate', 'chr1', 'start1', 'Parental', 'chr2', 'start2', 'bitscore', 'identity'])
        df = df.drop_duplicates()
        max_values = df.groupby('Duplicate')['bitscore'].transform('max')
        df_filtered_sub = df[max_values == df['bitscore']]
        max_values_identity = df_filtered_sub.groupby('Duplicate')['identity'].transform('max')
        df_filtered = df_filtered_sub[max_values_identity == df_filtered_sub['identity']]
        df_filtered.reset_index(drop=True, inplace=True)
        
        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, bitscore, identity = list(line)
            transposed_pair_list.append(dup_gene1 + '\t' + dup_gene2)
            transposed_pair_list.append(dup_gene2 + '\t' + dup_gene1)
            transposed_pair_file.write(dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + start1 + "\t" 
                                        + dup_gene2 + "\t" + output_prefix + "-" + chr2 + ":" + start2 + "\n")
            if homo_gn_md[dup_gene2] == 0:
                homo_gn_md[dup_gene2] = 4
            
            transposed_gene_list, transposed_gene_number, homo_gn_md = write_gene_to_wgd_transposed_file(dup_gene1, chr1, start1, homo_gn_md,
                                                                            transposed_gene_file, transposed_gene_list, output_prefix, transposed_gene_number, mode)

        transposed_gene_file.close()
        transposed_pair_file.close()
        transposed_pair_number = index + 1
        return homo_gn_md, transposed_pair_list, transposed_pair_number, transposed_gene_number, bitscore_dict, identity_dict

    @staticmethod
    def dispersed_process(data, dispersed_pair_file, dispersed_gene_file, output_prefix, homo_gn_md, used_gp_dict, bitscore_dict, identity_dict):
        dispersed_gene_list, dispersed_pair_list = [], []
        dispersed_gene_number, dispersed_pair_number = 0, 0
        dispersed_pre_df = []
        
        dispersed_pair_file = open(dispersed_pair_file, 'w')
        dispersed_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dispersed_gene_file = open(dispersed_gene_file, 'w')
        dispersed_gene_file.write("Duplicate" + "\t" + "Location" + "\n")

        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            if not used_gp_dict[ref_gene + '\t' + query_gene]:
                try:
                    avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(bitscore_dict[ref_gene + '\t' + query_gene])) / 2
                    avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(identity_dict[ref_gene + '\t' + query_gene])) / 2
                    dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, avg_bs, avg_id])
                except KeyError:
                    if ref_gene + '\t' + query_gene in bitscore_dict:
                        dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(bitscore_dict[ref_gene + '\t' + query_gene]), float(identity_dict[ref_gene + '\t' + query_gene])])
                    else:
                        dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(line[12]), float(line[13])])

        df = pd.DataFrame(dispersed_pre_df, columns=['Duplicate1', 'chr1', 'start1', 'Duplicate2', 'chr2', 'start2', 'bitscore', 'identity'])
        
        max_values = df.groupby('Duplicate2')['bitscore'].transform('max')
        df_filtered_sub = df[max_values == df['bitscore']]
        max_values_identity = df_filtered_sub.groupby('Duplicate2')['identity'].transform('max')
        df_filtered = df_filtered_sub[max_values_identity == df_filtered_sub['identity']]
        df_filtered.reset_index(drop=True, inplace=True)

        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, bitscore, identity = list(line)
            dispersed_pair_list.append(dup_gene1 + '\t' + dup_gene2)
            dispersed_pair_list.append(dup_gene2 + '\t' + dup_gene1)
            dispersed_pair_file.write(dup_gene2 + "\t" + output_prefix + "-" + chr2 + ":" + start2 + "\t" + 
                                      dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + start1 + "\n")
            if homo_gn_md[dup_gene1] not in [1, 2, 3, 4] and dup_gene1 not in dispersed_gene_list:
                dispersed_gene_list.append(dup_gene1)
                dispersed_gene_file.write(dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + start1 + "\n")
                dispersed_gene_number += 1
            if homo_gn_md[dup_gene2] not in [1, 2, 3, 4] and dup_gene2 not in dispersed_gene_list:
                dispersed_gene_list.append(dup_gene2)
                dispersed_gene_file.write(dup_gene2 + "\t" + output_prefix + "-" + chr2 + ":" + start2 + "\n")
                dispersed_gene_number += 1
        dispersed_pair_file.close()
        dispersed_gene_file.close()
        dispersed_pair_number = index + 1

        return dispersed_gene_list, dispersed_pair_number, dispersed_gene_number

    @staticmethod
    def merge_mode(homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md

    @staticmethod
    def singleton(gff, homo_gn_md, out_file, output_prefix):
        chr_gene_dict, chr_gene_list, _, _ = GffFile.readGff(gff)
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        singleton_number = 0
        for ch in chr_gene_list:
            for gn in chr_gene_list[ch]:
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + output_prefix + "-" + str(ch) + ":" + str(chr_gene_dict[ch][gn.name].start) + "\n")
                    singleton_number += 1
        ot_file.close()
        return singleton_number

    @staticmethod
    def sort_file(gene_file, pair_file):
        for file in gene_file:
            df = pd.read_csv(file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(file, header=True, index=False, sep='\t')
        for file in pair_file:
            df = pd.read_csv(file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(file, header=True, index=False, sep='\t')

    def run(self):
        #first
        file_judge(self.query_blast, self.query_gff_file, self.query_query_collinearity, self.query_ref_collinearity, self.out_directory)
        overwrite_judge([self.wgd_gene_file, self.tandem_gene_file, self.proximal_gene_file, self.transposed_gene_file, self.dispersed_gene_file,
                       self.wgd_pair_file, self.tandem_pair_file, self.proximal_pair_file, self.transposed_pair_file, self.dispersed_pair_file, 
                       self.singleton_gene_file, self.stats_file], self.overwrite)
        #second
        wgd_gene_list, wgd_pair_list, wgd_md, wgd_pair_number, wgd_gene_number = self.wgd_process(self.query_query_collinearity, self.wgd_pair_file, self.wgd_gene_file, self.output_prefix)
        print()
        print("wgd generation complete!")
        #third
        data = read_blast_gff(self.query_blast, self.query_gff_file)
        homo_gn_md, tandem_pair_lt, tandem_pair_number, tandem_gene_number, anc, used_gp_dict  = self.tandem_process(data, self.tandem_pair_file, self.tandem_gene_file, self.output_prefix)
        print()
        print("tandem generation complete!")
        #fourth
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, proximal_pair_list, proximal_pair_number, proximal_gene_number = self.proximal_process(data, self.proximal_pair_file, self.proximal_gene_file, self.output_prefix,
                                                                            self.proximal_max_distance, homo_gn_md)
        print()
        print("proximal generation complete!")
        #fifth
        homo_gn_md, transposed_pair_list, transposed_pair_number, transposed_gene_number, bitscore_dict, identity_dict = self.transposed_process(
            wgd_gene_list, self.query_ref_collinearity, data, self.transposed_pair_file, self.transposed_gene_file, self.output_prefix,
              self.seg_anc, self.proximal_max_distance, homo_gn_md, anc)
        print()
        print("transposed generation complete!")
        #sixth
        used_gp_dict = get_used_dict(wgd_pair_list, tandem_pair_lt, proximal_pair_list, transposed_pair_list, used_gp_dict)
        _, dispersed_pair_number, dispersed_gene_number = self.dispersed_process(data, self.dispersed_pair_file, self.dispersed_gene_file, self.output_prefix, homo_gn_md, used_gp_dict, bitscore_dict, identity_dict)
        print()
        print("dispersed generation complete!")
        #seventh
        singleton_number = self.singleton(self.query_gff_file, homo_gn_md, self.singleton_gene_file, self.output_prefix)
        print()
        print("singleton generation complete!")
        #eighth
        with open(self.stats_file, 'w') as f:
            f.write("Type" + "\t" + "Number" + "\n" +
                    "wgd.pairs" + "\t" + str(wgd_pair_number) + "\n" +
                    "tandem.pairs" + "\t" + str(tandem_pair_number) + "\n" +
                    "proximal.pairs" + "\t" + str(proximal_pair_number) + "\n" +
                    "transposed.pairs" + "\t" + str(transposed_pair_number) + "\n" +
                    "dispersed.pairs" + "\t" + str(dispersed_pair_number) + "\n" +
                    "wgd.genes" + "\t" + str(wgd_gene_number) + "\n" +
                    "tandem.genes" + "\t" + str(tandem_gene_number) + "\n" +
                    "proximal.genes" + "\t" + str(proximal_gene_number) + "\n" +
                    "transposed.genes" + "\t" + str(transposed_gene_number) + "\n" +
                    "dispersed.genes" + "\t" + str(dispersed_gene_number) + "\n" +
                    "singleton.genes" + "\t" + str(singleton_number))
        #nineth
        self.sort_file([self.wgd_gene_file, self.tandem_gene_file, self.proximal_gene_file, self.transposed_gene_file, self.dispersed_gene_file, self.singleton_gene_file],
                       [self.wgd_pair_file, self.tandem_pair_file, self.proximal_pair_file, self.transposed_pair_file, self.dispersed_pair_file])
        print()
        print("stats generation complete!")
        #plot
        base.ClsVis(self.stats_file, [self.stats_file + ".gene.png", self.stats_file + ".pair.png"]).run()



class ClassGeneUnique:
    def __init__(self, config_pra, parameter):
        for i in config_pra.sections():
            if i == 'classification':
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
        self.seg_anc = int(self.seg_anc)
        self.proximal_max_distance = int(self.proximal_max_distance)


        self.wgd_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.wgd.pairs")
        self.wgd_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.wgd.genes")
        self.tandem_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.tandem.pairs")
        self.tandem_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.tandem.genes")
        self.proximal_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.proximal.pairs")
        self.proximal_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.proximal.genes")
        self.transposed_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.transposed.pairs")
        self.transposed_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.transposed.genes")
        self.dispersed_pair_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.dispersed.pairs")
        self.dispersed_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.dispersed.genes")
        self.singleton_gene_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.singleton.genes")
        self.stats_file = os.path.join(self.out_directory, f"{self.output_prefix}-unique.stats")

    @staticmethod
    def wgd_process(query_query_collinearity, wgd_pair_file, wgd_gene_file, output_prefix):
        mode = 1
        wgd_gene_list, wgd_pair_list = [], []
        wgd_gene_number, wgd_pair_number = 0, 0
        wgd_md = {}
        wgd_pair_file = open(wgd_pair_file, 'w')
        wgd_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        wgd_gene_file = open(wgd_gene_file, 'w')
        wgd_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        with open(query_query_collinearity) as query_query_collinearity:
            _ = next(query_query_collinearity)
            _ = next(query_query_collinearity)
            for line in query_query_collinearity:
                # gene:SORBI_3010G080200	10	828	6837385	6840092	gene:SORBI_3002G223000	2	2281	61442049	61444165	+	0.7043
                if line.startswith("#"):
                    continue
                ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start = split_collinearity_line(line)

                if ref_gene + "\t" + query_gene not in wgd_pair_list:
                    wgd_pair_file.write(ref_gene + "\t" + output_prefix + "-" + ref_chr + ":" + ref_start
                                      + "\t" + query_gene + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\n")
                    wgd_pair_list.append(ref_gene + "\t" + query_gene)
                    wgd_pair_list.append(query_gene + "\t" + ref_gene)
                    wgd_pair_number += 1
                
                    wgd_gene_list, wgd_gene_number, wgd_md = write_gene_to_wgd_transposed_file(ref_gene, ref_chr, ref_start, 
                                                                        wgd_md, wgd_gene_file, wgd_gene_list, output_prefix, wgd_gene_number, mode)
                    wgd_gene_list, wgd_gene_number, wgd_md = write_gene_to_wgd_transposed_file(query_gene, query_chr, query_start,
                                                                        wgd_md, wgd_gene_file, wgd_gene_list, output_prefix, wgd_gene_number, mode)
        wgd_gene_file.close()
        wgd_pair_file.close()
        return wgd_gene_list, wgd_pair_list, wgd_md, wgd_pair_number, wgd_gene_number

    @staticmethod
    def tandem_process(data, tandem_pair_file, tandem_gene_file, output_prefix, wgd_gene_list):
        mode = 2
        homo_gn_md = {}
        anc = {}
        used_gp_dict = {}
        tandem_gene_list, tandem_pair_list = [], []
        tandem_gene_number, tandem_pair_number = 0, 0
        tandem_pair_file = open(tandem_pair_file, 'w')
        tandem_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        tandem_gene_file = open(tandem_gene_file, 'w')
        tandem_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            tandem_condition, homo_gn_md, anc, used_gp_dict = init_homo_gn_md(ref_gene, ref_chr, ref_order, query_gene, query_chr, query_order, homo_gn_md, anc, used_gp_dict)
            if tandem_condition:
                if ref_order - query_order < 0:
                    # here is more strict unique condition
                    if ref_gene not in wgd_gene_list and query_gene not in wgd_gene_list and \
                           ref_gene + "\t" + query_gene not in tandem_pair_list:
                        tandem_pair_list.append(ref_gene + "\t" + query_gene)
                        tandem_pair_list.append(query_gene + "\t" + ref_gene)
                        tandem_pair_file.write(ref_gene + "\t" + output_prefix + "-" + ref_chr + ":" + ref_start
                                         + "\t" + query_gene + "\t" + output_prefix + "-" + query_chr + ":"
                                         + query_start + "\n")
                        tandem_pair_number += 1
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(ref_gene, ref_chr, ref_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(query_gene, query_chr, query_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                else:
                    if ref_gene not in wgd_gene_list and query_gene not in wgd_gene_list and \
                            query_gene + "\t" + ref_gene not in tandem_pair_list:
                        tandem_pair_list.append(ref_gene + "\t" + query_gene)
                        tandem_pair_list.append(query_gene + "\t" + ref_gene)
                        tandem_pair_file.write(query_gene + "\t" + output_prefix + "-" + query_chr + ":"
                                         + query_start + "\t" + ref_gene + "\t" + output_prefix
                                         + "-" + ref_chr + ":" + ref_start + "\n")
                        tandem_pair_number += 1
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(ref_gene, ref_chr, ref_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
                        tandem_gene_list, tandem_gene_number, homo_gn_md = write_gene_to_file(query_gene, query_chr, query_start, homo_gn_md,
                                                                                            tandem_gene_file, tandem_gene_list, output_prefix, tandem_gene_number, mode)
        tandem_gene_file.close()
        tandem_pair_file.close()
        return homo_gn_md, tandem_pair_list, tandem_pair_number, tandem_gene_number, anc, used_gp_dict

    @staticmethod
    def proximal_process(data, proximal_pair_file, proximal_gene_file, output_prefix, proximal_threshold, homo_gn_md):
        mode = 3
        proximal_gene_list, proximal_pair_list = [], []
        proximal_gene_number, proximal_pair_number = 0, 0
        proximal_pre_df = []

        proximal_pair_file = open(proximal_pair_file, 'w')
        proximal_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        proximal_gene_file = open(proximal_gene_file, 'w')
        proximal_gene_file.write("Duplicate" + "\t" + "Location" + "\n")

        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            proximal_condition1 = ref_chr == query_chr
            proximal_condition2 = 1 < abs(ref_order - query_order) <= proximal_threshold
            if proximal_condition2 and proximal_condition1:
                if ref_order - query_order < 0:
                    if homo_gn_md[ref_gene] not in [1, 2] and homo_gn_md[query_gene] not in [1, 2]:
                        proximal_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start ,query_order - ref_order])
                else:
                    if homo_gn_md[ref_gene] not in [1, 2] and homo_gn_md[query_gene] not in [1, 2]:
                        proximal_pre_df.append([query_gene, query_chr, query_start, ref_gene, ref_chr, ref_start, ref_order - query_order])
        df = pd.DataFrame(proximal_pre_df, columns=['Duplicate1', 'chr1', 'start1', 'Duplicate2', 'chr2', 'start2', 'Distance'])
        df = df.drop_duplicates()
        min_values = df.groupby('Duplicate1')['Distance'].transform('min')
        df_filtered = df[min_values == df['Distance']]
        df_filtered.reset_index(drop=True, inplace=True)
        
        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, distance = list(line)
            proximal_pair_list.append(dup_gene1 + "\t" + dup_gene2)
            proximal_pair_list.append(dup_gene2 + "\t" + dup_gene1)
            proximal_pair_file.write(dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + ref_start
                                + "\t" + dup_gene2 + "\t" + output_prefix + "-" + query_chr + ":" + query_start + "\n")

            proximal_gene_list, proximal_gene_number, homo_gn_md = write_gene_to_file(dup_gene1, chr1, start1, homo_gn_md,
                                                                                        proximal_gene_file, proximal_gene_list, output_prefix, proximal_gene_number, mode)
            proximal_gene_list, proximal_gene_number, homo_gn_md = write_gene_to_file(dup_gene2, chr2, start2, homo_gn_md,
                                                                                        proximal_gene_file, proximal_gene_list, output_prefix, proximal_gene_number, mode)

        proximal_gene_file.close()
        proximal_pair_file.close()
        proximal_pair_number = index + 1

        return homo_gn_md, proximal_pair_list, proximal_pair_number, proximal_gene_number


    @staticmethod
    def transposed_process(wgd_gene_list, query_ref_collinearity, data, transposed_pair_file, transposed_gene_file, output_prefix, seg_anc, proximal_threshold, homo_gn_md, anc):
        mode = 4
        transposed_gene_list, transposed_pair_list = [], []
        transposed_gene_number, transposed_pair_number = 0, 0
        pre_df = []
        df = pd.DataFrame(data)
        df = df.iloc[:, [0, 6, 12, 13]]
        bitscore_dict = {str(row[0]) + "\t" + str(row[6]): row[12] for _, row in df.iterrows()}
        identity_dict = {str(row[0]) + "\t" + str(row[6]): row[13] for _, row in df.iterrows()}

        anc = get_anc(query_ref_collinearity, seg_anc, wgd_gene_list, anc)
        transposed_pair_file = open(transposed_pair_file, 'w')
        transposed_pair_file.write("Transposed" + "\t" + "Location1" + "\t" + "Parental" + "\t" + "Location2" + "\n")
        transposed_gene_file = open(transposed_gene_file, 'w')
        transposed_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            if ref_chr != query_chr:
                # TODO transposed_pre_df_unique: homo_gn_md[PARENTAL] not in [1,2,3]
                pre_df = transposed_pre_df_unique(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line[12], line[13], anc, homo_gn_md)
            else:
                if abs(ref_order - query_order) > proximal_threshold:
                    pre_df = transposed_pre_df_unique(ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, bitscore_dict, identity_dict, pre_df, line[12], line[13], anc, homo_gn_md)
        df = pd.DataFrame(pre_df, columns=['Duplicate', 'chr1', 'start1', 'Parental', 'chr2', 'start2', 'bitscore', 'identity'])
        df = df.drop_duplicates()
    
        max_values = df.groupby('Duplicate')['bitscore'].transform('max')
        df_filtered_sub = df[max_values == df['bitscore']]
        max_values_identity = df_filtered_sub.groupby('Duplicate')['identity'].transform('max')
        df_filtered = df_filtered_sub[max_values_identity == df_filtered_sub['identity']]
        df_filtered.reset_index(drop=True, inplace=True)

        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, bitscore, identity = list(line)
            transposed_pair_list.append(dup_gene1 + '\t' + dup_gene2)
            transposed_pair_list.append(dup_gene2 + '\t' + dup_gene1)
            transposed_pair_file.write(dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + start1 + "\t" 
                                        + dup_gene2 + "\t" + output_prefix + "-" + chr2 + ":" + start2 + "\n")
            if homo_gn_md[dup_gene2] == 0:
                homo_gn_md[dup_gene2] = 4
            
            transposed_gene_list, transposed_gene_number, homo_gn_md = write_gene_to_wgd_transposed_file(dup_gene1, chr1, start1, homo_gn_md,
                                                                            transposed_gene_file, transposed_gene_list, output_prefix, transposed_gene_number, mode)

        transposed_gene_file.close()
        transposed_pair_file.close()
        transposed_pair_number = index + 1
        return homo_gn_md, transposed_pair_list, transposed_pair_number, transposed_gene_number, bitscore_dict, identity_dict


    @staticmethod
    def dispersed_process(data, dispersed_pair_file, dispersed_gene_file, output_prefix, homo_gn_md, bitscore_dict, identity_dict):
        mode = 5
        dispersed_gene_list, dispersed_pair_list = [], []
        dispersed_gene_number, dispersed_pair_number = 0, 0
        dispersed_pre_df = []

        dispersed_pair_file = open(dispersed_pair_file, 'w')
        dispersed_pair_file.write("Duplicate1" + "\t" + "Location1" + "\t" + "Duplicate2" + "\t" + "Location2" + "\n")
        dispersed_gene_file = open(dispersed_gene_file, 'w')
        dispersed_gene_file.write("Duplicate" + "\t" + "Location" + "\n")
        
        for line in data:
            # refGene refChr refId refStart refEnd refStrand// queryGene queryChr queryId queryStart queryEnd  queryStrand// bitscore identity
            ref_gene, ref_chr, ref_order, ref_start, query_gene, query_chr, query_order, query_start = split_data_line(line)
            if homo_gn_md[ref_gene] not in [1, 2, 3, 4] and homo_gn_md[query_gene] not in [1, 2, 3, 4]:
                try:
                    avg_bs = (float(bitscore_dict[query_gene + '\t' + ref_gene]) + float(bitscore_dict[ref_gene + '\t' + query_gene])) / 2
                    avg_id = (float(identity_dict[query_gene + '\t' + ref_gene]) + float(identity_dict[ref_gene + '\t' + query_gene])) / 2
                    dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, avg_bs, avg_id])
                except KeyError:
                    if ref_gene + '\t' + query_gene in bitscore_dict:
                        dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(bitscore_dict[ref_gene + '\t' + query_gene]), float(identity_dict[ref_gene + '\t' + query_gene])])
                    else:
                        dispersed_pre_df.append([ref_gene, ref_chr, ref_start, query_gene, query_chr, query_start, float(line[12]), float(line[13])])
        df = pd.DataFrame(dispersed_pre_df, columns=['Duplicate1', 'chr1', 'start1', 'Duplicate2', 'chr2', 'start2', 'bitscore', 'identity'])
        max_values = df.groupby('Duplicate2')['bitscore'].transform('max')
        df_filtered_sub = df[max_values == df['bitscore']]
        max_values_identity = df_filtered_sub.groupby('Duplicate2')['identity'].transform('max')
        df_filtered = df_filtered_sub[max_values_identity == df_filtered_sub['identity']]
        df_filtered.reset_index(drop=True, inplace=True)

        for index, line in df_filtered.iterrows():
            dup_gene1, chr1, start1, dup_gene2, chr2, start2, bitscore, identity = list(line)
            dispersed_pair_list.append(dup_gene1 + '\t' + dup_gene2)
            dispersed_pair_list.append(dup_gene2 + '\t' + dup_gene1)
            dispersed_pair_file.write(dup_gene2 + "\t" + output_prefix + "-" + chr2 + ":" + start2 + "\t" + 
                                      dup_gene1 + "\t" + output_prefix + "-" + chr1 + ":" + start1 + "\n")
            
            
            dispersed_gene_list, dispersed_gene_number, homo_gn_md = write_gene_to_file(dup_gene1, chr1, start1, homo_gn_md,
                                                                                        dispersed_gene_file, dispersed_gene_list, output_prefix, dispersed_gene_number, mode)
            dispersed_gene_list, dispersed_gene_number, homo_gn_md = write_gene_to_file(dup_gene2, chr2, start2, homo_gn_md,
                                                                                        dispersed_gene_file, dispersed_gene_list, output_prefix, dispersed_gene_number, mode)
        dispersed_pair_file.close()
        dispersed_gene_file.close()
        dispersed_pair_number = index + 1

        return dispersed_gene_list, dispersed_pair_number, dispersed_gene_number

    @staticmethod
    def merge_mode(homo_gn_md, wgd_md):
        homo_gn_md.update(wgd_md)
        return homo_gn_md
 

    @staticmethod
    def singleton(gff, homo_gn_md, out_file, output_prefix):
        chr_gene_dict, chr_gene_list, _, _ = GffFile.readGff(gff)
        ot_file = open(out_file, 'w')
        ot_file.write("GeneId" + "\t" + "Location" + "\n")
        singleton_number = 0
        for ch in chr_gene_list:
            for gn in chr_gene_list[ch]:
                if gn.name not in homo_gn_md:
                    ot_file.write(str(gn.name) + "\t" + output_prefix + "-" + str(ch) + ":" + str(chr_gene_dict[ch][gn.name].start) + "\n")
                    singleton_number += 1
        ot_file.close()
        return singleton_number

    @staticmethod
    def sort_file(gene_file, pair_file):
        for file in gene_file:
            df = pd.read_csv(file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(file, header=True, index=False, sep='\t')
        for file in pair_file:
            df = pd.read_csv(file, sep="\t", header=0, index_col=None)
            df.sort_values(by=df.columns[1], inplace=True)
            df.to_csv(file, header=True, index=False, sep='\t')

    def run(self):
        #first
        file_judge(self.query_blast, self.query_gff_file, self.query_query_collinearity, self.query_ref_collinearity, self.out_directory)
        overwrite_judge([self.wgd_gene_file, self.tandem_gene_file, self.proximal_gene_file, self.transposed_gene_file, self.dispersed_gene_file,
                self.wgd_pair_file, self.tandem_pair_file, self.proximal_pair_file, self.transposed_pair_file, self.dispersed_pair_file, 
                self.singleton_gene_file, self.stats_file], self.overwrite)
        #second
        wgd_gene_list, _, wgd_md, wgd_pair_number, wgd_gene_number = self.wgd_process(self.query_query_collinearity, self.wgd_pair_file, self.wgd_gene_file, self.output_prefix)
        print()
        print("wgd generation complete!")
        #third
        data = read_blast_gff(self.query_blast, self.query_gff_file)
        homo_gn_md, _, tandem_pair_number, tandem_gene_number, anc, _  = self.tandem_process(data, self.tandem_pair_file, self.tandem_gene_file, self.output_prefix, wgd_gene_list)
        print()
        print("tandem generation complete!")
        #fourth
        homo_gn_md = self.merge_mode(homo_gn_md, wgd_md)
        homo_gn_md, _, proximal_pair_number, proximal_gene_number = self.proximal_process(data, self.proximal_pair_file, self.proximal_gene_file, self.output_prefix,
                                                                                        self.proximal_max_distance, homo_gn_md)
        print()
        print("proximal generation complete!")  
        #fifth
        homo_gn_md, _, transposed_pair_number, transposed_gene_number, bitscore_dict, identity_dict = self.transposed_process(
            wgd_gene_list, self.query_ref_collinearity, data, self.transposed_pair_file, self.transposed_gene_file, self.output_prefix,
              self.seg_anc, self.proximal_max_distance, homo_gn_md, anc)
        _, dispersed_pair_number, dispersed_gene_number = self.dispersed_process(data, self.dispersed_pair_file, self.dispersed_gene_file, self.output_prefix, homo_gn_md, bitscore_dict, identity_dict)
        print()
        print("transposed generation complete!")
        #sixth
        singleton_number = self.singleton(self.query_gff_file, homo_gn_md, self.singleton_gene_file, self.output_prefix)
        print()
        print("singleton generation complete!")
        #seventh
        with open(self.stats_file, 'w') as f:
            f.write("Type" + "\t" + "Number" + "\n" +
                    "wgd.pairs" + "\t" + str(wgd_pair_number) + "\n" +
                    "tandem.pairs" + "\t" + str(tandem_pair_number) + "\n" +
                    "proximal.pairs" + "\t" + str(proximal_pair_number) + "\n" +
                    "transposed.pairs" + "\t" + str(transposed_pair_number) + "\n" +
                    "dispersed.pairs" + "\t" + str(dispersed_pair_number) + "\n" +
                    "wgd.genes" + "\t" + str(wgd_gene_number) + "\n" +
                    "tandem.genes" + "\t" + str(tandem_gene_number) + "\n" +
                    "proximal.genes" + "\t" + str(proximal_gene_number) + "\n" +
                    "transposed.genes" + "\t" + str(transposed_gene_number) + "\n" +
                    "dispersed.genes" + "\t" + str(dispersed_gene_number) + "\n" +
                    "singleton.genes" + "\t" + str(singleton_number))
        #eighth
        self.sort_file([self.wgd_gene_file, self.tandem_gene_file, self.proximal_gene_file, self.transposed_gene_file, self.dispersed_gene_file],
                       [self.wgd_pair_file, self.tandem_pair_file, self.proximal_pair_file, self.transposed_pair_file, self.dispersed_pair_file])
        print()
        print("stats generation complete!")
        #plot
        base.ClsVis(self.stats_file, [self.stats_file + ".gene.png", self.stats_file + ".pair.png"]).run()

