import pandas as pd
import numpy as np
import subprocess
import os
import shutil
from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from multiprocessing import Pool, cpu_count
from tqdm import trange


class DuplicateKs:
    def __init__(self, config_pra, config_soft):
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'

        self.align_software = str(config_pra['ks']['align_software'])
        self.cds_file = str(config_pra['ks']['cds_file'])
        self.pep_file = str(config_pra['ks']['pep_file'])
        self.inuput_pair_dir = str(config_pra['ks']['inuput_pair_dir'])
        self.species_prefix = str(config_pra['ks']['species_prefix'])
        self.ouput_dir = str(config_pra['ks']['ouput_dir'])
        self.work_dir = os.getcwd()

        if cpu_count() > 32:
            self.process = 12
        else:
            self.process = cpu_count()-1

        if self.align_software == "muscle":
            self.muscle = config_soft['software']['muscle']
        if self.align_software == "mafft":
            self.mafft = config_soft['software']['mafft']
        self.yn00 = config_soft['software']['yn00']
        self.pal2nal = config_soft['software']['pal2nal']

    def get_input_pair_file_path(self):
        dup_pair_file_path = []
        suffix = [".wgd.pairs", ".tandem.pairs", ".proximal.pairs", ".transposed.pairs", ".dispersed.pairs"]
        for suf in suffix:
            dup_pair_file_path.append(os.path.join(self.inuput_pair_dir, self.species_prefix + suf))
        return dup_pair_file_path

    def get_ouput_pair_file_path(self):
        dup_pair_file_path = []
        suffix = [".wgd.pairs.ks", ".tandem.pairs.ks", ".proximal.pairs.ks", ".transposed.pairs.ks", ".dispersed.pairs.ks"]
        for suf in suffix:
            dup_pair_file_path.append(os.path.join(self.ouput_dir, self.species_prefix + suf))
        return dup_pair_file_path

    @staticmethod
    def read_pair_file(file):
        df = pd.read_csv(file, header=0, index_col=None, sep="\t")
        df = df.iloc[:, [0, 2]]
        df.iloc[:, :2] = df.iloc[:, :2].astype(str)
        df.columns = range(df.shape[1])
        df.index = df[0] + '\t' + df[1]
        return df

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

    @staticmethod
    def print_error(value):
        print("Process pool error, the cause of the error is :", value)

    def first_layer_run(self, dup_pair_file, ks_file):
        os.chdir(self.work_dir)
        if os.path.exists("tmp"):
            shutil.rmtree('tmp')
            os.mkdir('tmp')
        else:
            os.mkdir("tmp")
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        # query_to_ref index = "query+ref"
        df_pairs = self.read_pair_file(dup_pair_file)
        df_pairs = df_pairs[(df_pairs[0].isin(cds.keys())) & (df_pairs[1].isin(
            cds.keys())) & (df_pairs[0].isin(pep.keys())) & (df_pairs[1].isin(pep.keys()))]
        # all pairs
        pairs = df_pairs[[0, 1]].to_numpy()
        n = int(np.ceil(len(pairs) / float(self.process)))
        os.chdir(os.path.join(self.work_dir, 'tmp'))
        pool = Pool(int(self.process))
        for i in range(int(self.process)):
            if os.path.exists('process_' + str(i)):
                pass
            else:
                os.mkdir('process_' + str(i))
            if i < int(self.process)-1:
                sub_pr = pairs[i*n:i*n+n]
            else:
                sub_pr = pairs[i*n:]
            pool.apply_async(self.secondary_layer_run, args=(sub_pr, ks_file, i, cds, pep), error_callback=self.print_error)
        pool.close()
        pool.join()
        shutil.rmtree(os.path.join(self.work_dir, 'tmp'))

    def secondary_layer_run(self, pairs, ks_file, i, cds, pep):
        os.chdir(os.path.join(self.work_dir, 'tmp', 'process_' + str(i)))
        # ks_File is absolute path
        if os.path.exists(ks_file):
            # ks = pd.read_csv(ks_file, sep='\t', low_memory=False)
            # ks.iloc[:, :2] = ks.iloc[:, :2].astype(str)
            # combined_list1 = (ks.iloc[:, 0] + '\t' + ks.iloc[:, 1]).tolist()
            # combined_list2 = (ks.iloc[:, 1] + '\t' + ks.iloc[:, 0]).tolist()
            # jd_lt = set(combined_list1+combined_list2)
            # df_pairs = df_pairs[~df_pairs.index.isin(jd_lt)]
            ks_file_handle = open(ks_file, 'a+')
        else:
            ks_file_handle = open(ks_file, 'w')
            ks_file_handle.write(
                '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00'])+'\n')
        if i == int(self.process) - 1:
            for i in trange(len(pairs)):
                k = pairs[i]
                SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
                SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
                kaks = self.pair_kaks(k)
                if not kaks:
                    continue
                ks_file_handle.write('\t'.join([str(i) for i in list(k) + list(kaks)]) + '\n')
            ks_file_handle.close()
        else:
            for k in pairs:
                SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
                SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
                kaks = self.pair_kaks(k)
                if not kaks:
                    continue
                ks_file_handle.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
            ks_file_handle.close()
        for file in (self.pair_pep_file, self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file,
                     '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            try:
                os.remove(file)
            except OSError:
                pass

    def run(self):
        path_list = self.get_input_pair_file_path()
        out_path_list = self.get_ouput_pair_file_path()
        for i in range(len(path_list)):
            self.first_layer_run(path_list[i], out_path_list[i])
