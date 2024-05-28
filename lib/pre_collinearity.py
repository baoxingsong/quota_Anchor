import subprocess
from . import longestPeps, combineBlastAndStrandInformation


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
        self.run_gff_read_get_protein(self.query_genome_seq, self.query_gff_file, self.output_query_pep_seq)
        self.run_gff_read_get_protein(self.ref_genome_seq, self.ref_gff_file, self.output_ref_pep_seq)

        longestPeps.longestPeps(self.query_gff_file, self.query_genome_seq, self.output_query_pep_seq,
                                self.out_query_longest_pep_name)
        longestPeps.longestPeps(self.ref_gff_file, self.ref_genome_seq, self.output_ref_pep_seq, self.out_ref_longest_pep_name)

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
