import os
import sys
from . import base
import subprocess    
from . import combineBlastAndStrandInformation


class Prepare:
    def __init__(self, config_pra, config_soft, parameter):
        self.overwrite = False
        self.strand = "plus"
        self.skip_blast = parameter.skip_blast
        self.bitscore = 100
        self.align_length = 0
        self.outfmt = 6
        self.thread = 6
        self.dtype = "prot"
        # config
        for i in config_pra.sections():
            if i == "combineBlastAndStrand":
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])

        if not self.skip_blast:
            if "align" in config_pra.sections():
                self.align = config_pra["align"]["align"]
            if hasattr(self, 'align') and  self.align == "diamond" and "diamond" in config_pra.sections():
                for key in config_pra['diamond']:
                    setattr(self, key, config_pra["diamond"][key])
            if hasattr(self, 'align') and (self.align == "blastp" or self.align == "blastn") and "blast" in config_pra.sections():
                for key in config_pra['blast']:
                    setattr(self, key, config_pra["blast"][key])
        # command line
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        # software
        softwawre_list = []
        for key in config_soft['software']:
            softwawre_list.append(key)
            setattr(self, key, config_soft["software"][key])
        
        print()
        for key, value in vars(self).items():
            if key not in softwawre_list and key != "conf":
                if not self.skip_blast:
                    if self.align == "diamond":
                        key_list = ["thread", "outfmt", "dtype", "strand"]
                        if key not in key_list:
                            print(key, "=", value)
                    if self.align == "blastp" or self.align == "blastn":
                        if key != "strand":
                            print(key, "=", value)
                else:
                    blast_key_list = ["align", "database_name", "output_blast_result", "strand",
                                      "ref_seq", "query_seq", "max_target_seqs", "evalue", "thread", "outfmt", "dtype"]
                    if key not in blast_key_list:
                        print(key, "=", value)
        print()
        
    def diamond_makedb(self, ref_protein, database):
        command_line = [self.diamond, 'makedb', '--in', ref_protein, '--db', database]
        try:
            diamond_db_log = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print(f"{diamond_db_log}, class Prepare's function diamond_makedb failed")

    def run_diamond_blastp(self, database, query_file, blast_file, max_target, e_value):
        command_line = [self.diamond, 'blastp', '--db', database, '-q', query_file, '-o', blast_file, '-k', max_target, '-e', e_value]
        try:
            diamond_blastp_log = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'{diamond_blastp_log}, running diamond blastp occurs error')
    
    def mkblastdb(self, ref_seq, database, dtype):
        command_line = [self.makeblastdb, '-in', ref_seq, '-dbtype', dtype, '-out', database]
        try:
            blast_db_log = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print(f"{blast_db_log}, class Prepare's function makeblastdb failed")

    def run_blastn(self, blast_database, query_file, blast_file, e_value, thread, outfmt, max_target_seqs, strand):
        command_line = [self.align, '-query', query_file, '-out', blast_file, '-evalue', e_value,
                        '-db', blast_database, '-num_threads', thread, '-max_target_seqs', max_target_seqs, '-outfmt', outfmt, '-strand', strand]
        try:
            blast_log = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'{blast_log}, running blastn occurs error')

    def run_blastp(self, blast_database, query_file, blast_file, e_value, thread, outfmt, max_target_seqs):
        command_line = [self.align, '-query', query_file, '-out', blast_file, '-evalue', e_value,
                        '-db', blast_database, '-num_threads', thread, '-max_target_seqs', max_target_seqs, '-outfmt', outfmt]
        try:
            blast_log = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'{blast_log}, running blastp occurs error')
    def run_all_process(self):
        flag = True
        if not self.skip_blast:
            flag = False
            base.file_empty(self.ref_seq)
            base.file_empty(self.query_seq)
            base.output_file_parentdir_exist(self.output_file, self.overwrite)
            
            if self.align == "diamond":
                self.diamond_makedb(self.ref_seq, self.database_name)
                self.run_diamond_blastp(self.database_name, self.query_seq, self.output_blast_result,
                                        self.max_target_seqs, self.evalue)
            if self.align == "blastn":
                
                self.mkblastdb(self.ref_seq, self.database_name, self.dtype)
                self.run_blastn(self.database_name, self.query_seq, self.output_blast_result,
                                self.evalue, str(self.thread), str(self.outfmt), self.max_target_seqs, self.strand)
            if self.align == "blastp":
                self.mkblastdb(self.ref_seq, self.database_name, self.dtype)
                self.run_blastp(self.database_name, self.query_seq, self.output_blast_result,
                                self.evalue, str(self.thread), str(self.outfmt), self.max_target_seqs)
        base.file_empty(self.ref_gff_file)
        base.file_empty(self.query_gff_file)
        base.file_empty(self.blast_file)
        if flag:
            base.output_file_parentdir_exist(self.output_file, self.overwrite)
        combineBlastAndStrandInformation.anchorwave_quota(self.ref_gff_file, self.query_gff_file, self.blast_file,
                                                          self.output_file, self.bitscore, self.align_length)
        base.file_empty(self.output_file)
