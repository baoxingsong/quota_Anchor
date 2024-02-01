# This script is modified from WGDI's ks.py

import subprocess
import pandas as pd
import os
import numpy as np
from Bio import SeqIO
from Bio.Phylo.PAML import yn00


def read_collinearity(collinearity):
    block_len = []
    with open(collinearity) as f:
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[2].split(sep="N=")[1]))
    # read anchorWave output collinearity file as a dataframe, every line of the result occurs once different from wgdi.
    collinearity_df = pd.read_table(collinearity, comment="#", header=0)

    # get two gene lists
    ref_gene_list = list(collinearity_df.loc[:, "refGene"])  # zea mays (ref) collinearity gene list
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  # sorghum (query) collinearity gene list
    assert (len(ref_gene_list) == len(query_gene_list))
    # ref_to_query = dict(zip(ref_gene_list, query_gene_list))
    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    df = pd.DataFrame(query_to_ref)
    df[0] = df[0].astype(str)
    df[1] = df[1].astype(str)
    df.index = df[0] + ',' + df[1]
    return df, block_len


class Ks:
    def __init__(self, config_pra, config_soft):
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'

        self.muscle = config_soft['software']['muscle']
        self.mafft = config_soft['software']['mafft']
        self.align_software = config_pra['ks']['align_software']
        self.yn00 = config_soft['software']['yn00']
        self.pal2nal = config_soft['software']['pal2nal']

        self.output_query_cds_seq = config_pra['gffread']['output_query_cds_seq']
        self.output_ref_cds_seq = config_pra['gffread']['output_ref_cds_seq']

        self.collinearity = config_pra['ks']['collinearity']
        self.cds_file = config_pra['ks']['cds_file']
        self.pep_file = config_pra['ks']['pep_file']
        self.ks_file = config_pra['ks']['ks_file']
        self.type = config_pra['ks']['type']

    def pair_kaks(self, k):
        self.align()
        pal = self.run_pal2nal()
        if not pal:
            return []
        kaks = self.run_yn00()
        if kaks is None:
            return []
        kaks_new = [kaks[k[0]][k[1]]['NG86']['dN'], kaks[k[0]][k[1]]['NG86']
                    ['dS'], kaks[k[0]][k[1]]['YN00']['dN'], kaks[k[0]][k[1]]['YN00']['dS']]
        return kaks_new

    def align(self):
        if self.align_software == 'mafft':
            command_line = [self.mafft, "--quiet", "--auto", self.pair_pep_file]
            try:
                with open(self.prot_align_file, 'w') as output_file:
                    subprocess.run(command_line, stdout=output_file, check=True)
            except subprocess.CalledProcessError as e:
                print("class_Ks's function align failed: ", e)
        if self.align_software == 'muscle':
            command_line = [self.muscle, '-align', self.pair_pep_file, '-output', self.prot_align_file]
            try:
                subprocess.run(command_line, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print("class_Ks's function align failed: ", e)

    def run_pal2nal(self):
        command_line = [self.pal2nal, self.prot_align_file, self.pair_cds_file, '-output', 'paml', '-nogap']
        try:
            with open(self.mrtrans, 'w') as output_file:
                subprocess.run(command_line, stdout=output_file, check=True)
        except subprocess.CalledProcessError as e:
            print("class_Ks's function run_pal2nal failed: ", e)
            return None
        return True

    def run_yn00(self):
        yn = yn00.Yn00()
        yn.alignment = self.mrtrans
        yn.out_file = self.pair_yn
        yn.set_options(icode=0, commonf3x4=0, weighting=0, verbose=1)
        try:
            run_result = yn.run(command=self.yn00)
        except Exception as e:
            print("class ks function run_yn00 failed", e)
            run_result = None
        return run_result

    def sub_run(self):
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        # query_to_ref index = "query+ref"
        df_pairs, _ = read_collinearity(self.collinearity)
        if os.path.exists(self.ks_file):
            ks = pd.read_csv(self.ks_file, sep='\t')
            ks = ks.drop_duplicates()
            kscopy = ks.copy()
            names = ks.columns.tolist()
            names[0], names[1] = names[1], names[0]
            kscopy.columns = names
            ks = pd.concat([ks, kscopy])
            ks['id'] = ks['id1']+','+ks['id2']
            # df_pairs.index length equal to len(df_pairs)
            # query_ref intersect1d query_ref and ref_query in the interspecific, ks = pd.concat([ks, kscopy]) is useless
            # when len(df_pairs.index) "fake_less than" len(ks['id'].to_numpy()) this is useless
            # when coll have five duplicate and ks have one duplicate, don't add
            # when coll have more pair than ks , this maybe useful
            # intra specific, add A B or B A
            df_pairs.drop(np.intersect1d(df_pairs.index,
                                         ks['id'].to_numpy()), inplace=True)
            ks_file = open(self.ks_file, 'a+')
        else:
            # maybe have duplicate and length equal to len(df_pairs) in the interspecific
            # maybe don't have duplicate and length don't equal to len(df_pairs) in the intraspecific
            # only retain A B or B A
            ks_file = open(self.ks_file, 'w')
            ks_file.write(
                '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00'])+'\n')
        df_pairs = df_pairs[(df_pairs[0].isin(cds.keys())) & (df_pairs[1].isin(
            cds.keys())) & (df_pairs[0].isin(pep.keys())) & (df_pairs[1].isin(pep.keys()))]
        pairs = df_pairs[[0, 1]].to_numpy()
        if len(pairs) > 0 and self.type == 'intra':
            allpairs = []
            pair_hash = {}
            for k in pairs:
                if k[0] == k[1]:
                    continue
                elif k[0]+','+k[1] in pair_hash or k[1]+','+k[0] in pair_hash:
                    continue
                else:
                    pair_hash[k[0]+','+k[1]] = 1
                    pair_hash[k[1]+','+k[0]] = 1
                    allpairs.append(k)
            pairs = allpairs
        for k in pairs:
            cds_gene1, cds_gene2 = cds[k[0]], cds[k[1]]
            cds_gene1.id, cds_gene2.id = 'gene1', 'gene2'
            pep_gene1, pep_gene2 = pep[k[0]], pep[k[1]]
            pep_gene1.id, pep_gene2.id = 'gene1', 'gene2'
            SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
            SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
            kaks = self.pair_kaks(['gene1', 'gene2'])
            if not kaks:
                continue
            ks_file.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
        ks_file.close()
        for file in (self.pair_pep_file, self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file,
                     '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            try:
                os.remove(file)
            except OSError:
                pass
