import subprocess
from . import longestPeps, combineBlastAndStrandInformation
from multiprocessing import Pool, cpu_count
import os
import sys
from . import base


class Prepare:
    def __init__(self, config_pra, config_soft):
        self.config_pra = config_pra

        self.gffread = config_soft['software']['gffread']

        self.query_genome_seq = config_pra['gffread']['query_genome_seq']
        self.query_gff_file = config_pra['gffread']['query_gff_file']
        self.output_query_pep_seq = config_pra['gffread']['output_query_pep_seq']
        self.ref_genome_seq = config_pra['gffread']['ref_genome_seq']
        self.ref_gff_file = config_pra['gffread']['ref_gff_file']
        self.output_ref_pep_seq = config_pra['gffread']['output_ref_pep_seq']
        if config_pra['gffread']['use_S_parameter']:
            self.S = '-S'
        else:
            self.S = ""

        self.out_query_longest_pep_name = config_pra['longest_pep']['out_query_longest_pep_name']
        self.out_ref_longest_pep_name = config_pra['longest_pep']['out_ref_longest_pep_name']

        if config_pra['align']['align'] == "diamond":
            self.diamond = config_soft['software']['diamond']
            self.database_name = config_pra['diamond']['database_name']
            self.blast_result = config_pra['diamond']['output_blast_result']
            self.max_target_seqs = config_pra['diamond']['max_target_seqs']
            self.evalue = config_pra['diamond']['evalue']

        if config_pra['align']['align'] == "blastp":
            self.blastp = config_soft['software']['blastp']
            self.makeblastdb = config_soft['software']['makeblastdb']
            self.blast_database = config_pra['blastp']['database_name']
            self.dtype = config_pra['blastp']['dtype']
            self.max_target_seqs = config_pra['blastp']['max_target_seqs']
            self.e_value = config_pra['blastp']['evalue']
            self.num_thread = config_pra['blastp']['thread']
            self.blast_result = config_pra['blastp']['output_blast_result']
            self.outfmt = config_pra['blastp']['outfmt']

        self.out_file = config_pra['combineBlastAndStrand']['out_file']
        self.bitscore = config_pra['combineBlastAndStrand']['bitscore']
        self.align_length = config_pra['combineBlastAndStrand']['align_length']

    @staticmethod
    def file_empty(file_path):
        if os.path.isfile(file_path):
            if os.path.getsize(file_path) > 0:
                pass
            else:
                sys.exit(1)
        else:
            error_message = f"{file_path} don't exist"
            print(error_message)
            raise FileNotFoundError(error_message)

    def run_gff_read_get_protein(self, fasta, gff, output_protein_file):
        command_line = [self.gffread, '-g', fasta, '-y', output_protein_file, gff, self.S]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print("class_Collinearity's function run_gffread failed: ", e)

    def diamond_make_db(self, ref_protein, database):
        command_line = [self.diamond, 'makedb', '--in', ref_protein, '--db', database]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print("class_Collinearity's function diamond_make_db failed")

    def mkblastdb(self, ref_protein, database, dtype):
        command_line = [self.makeblastdb, '-in', ref_protein, '-dbtype', dtype, '-out', database]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print("class_Collinearity's function makeblastdb failed")

    def run_diamond_blastp(self, database, query_file, blast_file, max_target, e_value):
        command_line = [self.diamond, 'blastp', '--db', database, '-q', query_file, '-o', blast_file, '-k', max_target, '-e', e_value]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'Error running Diamond blastp: {e}')

    def run_blastp(self, blast_database, query_file, blast_file, e_value, thread, outfmt, max_target_seqs):
        command_line = [self.blastp, '-query', query_file, '-out', blast_file, '-evalue', e_value,
                        '-db', blast_database, '-num_threads', thread, '-max_target_seqs', max_target_seqs, '-outfmt', outfmt]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'Error running blastp: {e}')

    def run_all_process(self):
        self.file_empty(self.query_genome_seq)
        self.file_empty(self.query_gff_file)
        self.run_gff_read_get_protein(self.query_genome_seq, self.query_gff_file, self.output_query_pep_seq)
        self.file_empty(self.output_query_pep_seq)

        self.file_empty(self.ref_genome_seq)
        self.file_empty(self.ref_gff_file)
        self.run_gff_read_get_protein(self.ref_genome_seq, self.ref_gff_file, self.output_ref_pep_seq)
        self.file_empty(self.output_ref_pep_seq)

        longestPeps.longestPeps(self.query_gff_file, self.query_genome_seq, self.output_query_pep_seq,
                                self.out_query_longest_pep_name)
        self.file_empty(self.out_query_longest_pep_name)
        longestPeps.longestPeps(self.ref_gff_file, self.ref_genome_seq, self.output_ref_pep_seq, self.out_ref_longest_pep_name)
        self.file_empty(self.out_ref_longest_pep_name)

        if self.config_pra['align']['align'] == "diamond":
            self.diamond_make_db(self.out_ref_longest_pep_name, self.database_name)
            self.run_diamond_blastp(self.database_name, self.out_query_longest_pep_name, self.blast_result,
                                    self.max_target_seqs, self.evalue)
        if self.config_pra['align']['align'] == "blastp":
            self.mkblastdb(self.out_ref_longest_pep_name, self.blast_database, self.dtype)
            self.run_blastp(self.blast_database, self.out_query_longest_pep_name, self.blast_result,
                            self.e_value, self.num_thread, self.outfmt, self.max_target_seqs)

        combineBlastAndStrandInformation.anchorwave_quota(self.ref_gff_file, self.query_gff_file, self.blast_result,
                                                          self.out_file, self.bitscore, self.align_length)
        self.file_empty(self.out_file)


class Longest:
    def __init__(self, config_pra, config_soft):
        self.config_pra = config_pra

        self.gffread = config_soft['software']['gffread']

        self.genome_seq = config_pra['gffread']['genome_seq']
        self.gff_file = config_pra['gffread']['gff_file']
        self.out_pep_seq = config_pra['gffread']['out_pep_seq']

        if config_pra['gffread']['use_S_parameter']:
            self.S = '-S'
        else:
            self.S = ""

        self.out_longest_pep_name = config_pra['longest_pep']['out_longest_pep_name']
        self.process_number = int(config_pra['longest_pep']['thread'])

    def run_gff_read_get_protein(self, fasta, gff, output_protein_file):
        command_line = [self.gffread, '-g', fasta, '-y', output_protein_file, gff, self.S]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print("class_Collinearity's function run_gffread failed: ", e)
    
    def sub_run(self, sub_run_number, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx):

        pool = Pool(sub_run_number)
        for i in range(idx, idx + sub_run_number):
            base.file_empty(genome_list[i])
            base.file_empty(gff_list[i])
            self.run_gff_read_get_protein(genome_list[i], gff_list[i], out_pep_list[i])
            base.file_empty(out_pep_list[i])
            pool.apply_async(longestPeps.longestPeps, args=(gff_list[i], genome_list[i], out_pep_list[i], out_longest_pep_name_list[i]))
        pool.close()
        pool.join()
    
    def run_all_process(self):
        genome_list = base.split_conf(self.genome_seq, ",")
        gff_list = base.split_conf(self.gff_file, ",")
        out_pep_list = base.split_conf(self.out_pep_seq, ",")
        out_longest_pep_name_list = base.split_conf(self.out_longest_pep_name, ",")
        try:
            assert len(genome_list) == len(gff_list) == len(out_pep_list) == len(out_longest_pep_name_list)
        except AssertionError as e:
            print(f"AssertionError: {e} please check your separator for config file!")

        ideal_process_number = len(gff_list)
        sub_run_number = min(ideal_process_number, cpu_count(), self.process_number)
        interger_pool = int(ideal_process_number / sub_run_number) 
        final_pool = ideal_process_number % sub_run_number
        
        idx = 0
        for _ in range(interger_pool):
            self.sub_run(sub_run_number, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx)
            idx += sub_run_number
        if not final_pool:
            self.sub_run(final_pool, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx)

        for i in range(ideal_process_number):
            base.file_empty(out_longest_pep_name_list[i])

