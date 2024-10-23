# This script is refers to wgdi 's ks.py and famcircle 's ks.py and biopython 's paml module.
# https://github.com/lkiko/famCircle
# https://github.com/SunPengChuan/wgdi/blob/master/wgdi/ks.py
# https://github.com/biopython/biopython/blob/master/Bio/Phylo/PAML

import subprocess
import pandas as pd
import os, re
import traceback
from . import base
import shutil
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from alive_progress import alive_bar


def read_collinearity(collinearity):
    block_len = []
    with open(collinearity) as f:
        _ = next(f)
        for line in f:
            if line.startswith('#'):
                block_len.append(int(line.split()[2].split(sep="N=")[1]))
    collinearity_df = pd.read_table(collinearity, comment="#", header=0, low_memory=False)

    ref_gene_list = list(collinearity_df.loc[:, "refGene"])  
    query_gene_list = list(collinearity_df.loc[:, "queryGene"])  
    assert (len(ref_gene_list) == len(query_gene_list))

    query_to_ref = list(zip(query_gene_list, ref_gene_list))
    df = pd.DataFrame(query_to_ref)
    df[0] = df[0].astype(str)
    df[1] = df[1].astype(str)
    df.index = df[0] + ',' + df[1]
    return df, block_len


class Ks:
    def __init__(self, config_pra, config_soft,parameter):
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'
        self.yn00_ctl = "yn00.ctl"
        self.overwrite = False
        for i in config_pra.sections():
            if i == 'ks':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "conf" and key not in ["prot_align_file", "mrtrans", "pair_yn", "pair_pep_file", "pair_cds_file"]:
                print(key, "=", value)
        print()

        self.work_dir = os.getcwd()
        for key in config_soft['software']:
            setattr(self, key, config_soft["software"][key])

        if cpu_count() > 32:
            self.process = 12
        else:
            self.process = cpu_count() - 1

    def pair_kaks(self, k):
        self.align()
        pal = self.run_pal2nal()
        if not pal:
            return []
        self.run_yn00()
        NG86, YN00 = self.read_yn00_result()
        kaks = [NG86, YN00]
        return kaks

    def align(self):
        if self.align_software == 'mafft':
            command_line = [self.mafft, "--quiet", "--auto", self.pair_pep_file]
            try:
                with open(self.prot_align_file, 'w') as output_file:
                    subprocess.run(command_line, stdout=output_file, stderr=subprocess.DEVNULL, check=True)
            except subprocess.CalledProcessError as e:
                print("class Ks's function align failed: ", e)
        if self.align_software == 'muscle':
            command_line = [self.muscle, '-align', self.pair_pep_file, '-output', self.prot_align_file]
            try:
                subprocess.run(command_line, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print("class Ks's function align failed: ", e)

    def run_pal2nal(self):
        command_line = [self.pal2nal, self.prot_align_file, self.pair_cds_file, '-output', 'paml', '-nogap']
        try:
            with open(self.mrtrans, 'w') as output_file:
                subprocess.run(command_line, stdout=output_file, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError as e:
            print("class Ks's function run_pal2nal failed: ", e)
            return None
        return True

    def run_yn00(self):
        command_line = [self.yn00]
        try:
            result = subprocess.run(command_line, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError as e:
            print("class Ks's function run_yn00 failed: ", e)
        return

    def read_yn00_result(self):
        with open(self.pair_yn) as f:
            content = f.readlines()
            sub_content = []
            for line in content:
                if "(A) Nei-Gojobori (1986) method" in line:
                    a_index = content.index(line)
                if "(C) LWL85, LPB93 & LWLm methods" in line:
                    c_index = content.index(line)
            sub_content = content[a_index:c_index]
            for line in sub_content:
                matrix_row_res = re.match(r"^([^\s]+?)(\s+-?\d+\.\d+.*$|\s*$|-1.0000\s*\(.*$)", line)
                if matrix_row_res is not None:
                    line_floats_res = re.findall(r"-*\d+\.\d+", matrix_row_res.group(2))
                    line_floats = [float(val) for val in line_floats_res]
                    for i in range(0, len(line_floats), 3):
                        NG86 = {}
                        NG86["omega"] = line_floats[i]
                        NG86["dN"] = line_floats[i + 1]
                        NG86["dS"] = line_floats[i + 2]

                line_floats_res = re.findall(r"-*\d+\.\d+", line)
                line_floats = [float(val) for val in line_floats_res]
                row_res = re.match(r"\s+(\d+)\s+(\d+)", line)
                if row_res is not None:
                    YN00 = {}
                    YN00["t"] = line_floats[2]
                    YN00["kappa"] = line_floats[3]
                    YN00["omega"] = line_floats[4]
                    YN00["dN"] = line_floats[5]
                    YN00["dS"] = line_floats[7]
        return NG86, YN00

    def first_layer_run(self):
        base.file_empty(self.collinearity)
        base.file_empty(self.cds_file)
        base.file_empty(self.pep_file)
        base.output_file_parentdir_exist(self.ks_file, self.overwrite)
        if os.path.exists("tmp"):
            shutil.rmtree('tmp')
            os.mkdir('tmp')
        else:
            os.mkdir("tmp")
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        # query_to_ref index = "query+ref"
        df_pairs, _ = read_collinearity(self.collinearity)
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
        n = int(np.ceil(len(pairs) / float(self.process)))
        os.chdir(os.path.join(self.work_dir, 'tmp'))
        pool = Pool(int(self.process))
        for i in range(int(self.process)):
            if os.path.exists('process_' + str(i)):
                pass
            else:
                os.mkdir('process_' + str(i))
            if i < int(self.process) - 1:
                sub_pr = pairs[i * n:i * n + n]
            else:
                sub_pr = pairs[i * n:]
            pool.apply_async(self.secondary_layer_run, args=(sub_pr, i, cds, pep))
        pool.close()
        pool.join()
        shutil.rmtree(os.path.join(self.work_dir, 'tmp'))

    def secondary_layer_run(self, pairs, i, cds, pep):
        os.chdir(self.work_dir)
        if os.path.exists(self.ks_file):
            ks_file_handle = open(self.ks_file, 'a+')
        else:
            ks_file_handle = open(self.ks_file, 'w')
            ks_file_handle.write(
                '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'omega_NG86', 'ka_YN00', 'ks_YN00', 'omega_YN00', 't'])+'\n')
        os.chdir(os.path.join(self.work_dir, 'tmp', 'process_' + str(i)))
        with open(self.yn00_ctl, 'w') as ctl:
            ctl.write("seqfile = " + self.mrtrans + "\n" + 
                      "outfile = " + self.pair_yn + "\n" +
                      "verbose = 1" + "\n" +
                      "icode = 0" + "\n" +
                      "weighting = 0" + "\n" +
                      "commonf3x4 = 0")
        if i == int(self.process) - 1:
            with alive_bar(len(pairs), title="quota_Anchor", bar="bubbles", spinner="waves") as bar:
                for i in range(len(pairs)):
                    k = pairs[i]
                    SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
                    SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
                    kaks = self.pair_kaks(k)
                    if not kaks:
                        continue
                    kaks_list = [kaks[0]["dN"], kaks[0]["dS"], kaks[0]["omega"],
                                  kaks[1]["dN"], kaks[1]["dS"], kaks[1]["omega"], kaks[1]["t"]]
                    ks_file_handle.write('\t'.join([str(i) for i in list(k)+kaks_list])+'\n')
                    bar()
                ks_file_handle.close()
        else:
            for k in pairs:
                SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
                SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
                kaks = self.pair_kaks(k)
                if not kaks:
                    continue
                kaks_list = [kaks[0]["dN"], kaks[0]["dS"], kaks[0]["omega"],
                                  kaks[1]["dN"], kaks[1]["dS"], kaks[1]["omega"], kaks[1]["t"]]
                ks_file_handle.write('\t'.join([str(i) for i in list(k)+kaks_list])+'\n')
            ks_file_handle.close()
        for file in (self.pair_pep_file, self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file,
                     '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            try:
                os.remove(file)
            except FileNotFoundError as e1:
                print("ks class remove file error, file not be found:", e1)
            except OSError as e2:
                print("ks class remove file error, the path is a directory:", e2)
                
