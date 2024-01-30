import subprocess
import longestPeps
import combineBlastAndStrandInformation


class Collinearity:
    def __init__(self, config_pra, config_soft):
        self.gffread = config_soft['software']['gffread']
        self.AnchorWave = config_soft['software']['AnchorWave']
        self.diamond = config_soft['software']['diamond']
        self.blastp = config_soft['software']['blastp']

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

        self.database_name = config_pra['diamond']['database_name']
        self.output_blast_result = config_pra['diamond']['output_blast_result']
        self.max_target_seqs = config_pra['diamond']['max_target_seqs']
        self.evalue = config_pra['diamond']['evalue']

        self.input_file_name = config_pra['AnchorWave']['input_file_name']
        self.R = config_pra["AnchorWave"]["R"]
        self.Q = config_pra["AnchorWave"]["Q"]
        self.maximum_gap_size = config_pra["AnchorWave"]["maximum_gap_size"]
        self.delete_tandem = config_pra["AnchorWave"]["delete_tandem"]
        self.tandem_dis = config_pra["AnchorWave"]["tandem_dis"]
        self.output_coll_name = config_pra["AnchorWave"]["output_coll_name"]

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

    def run_diamond_blastp(self, database, query_file, blast_file, max_target, e_value):
        command_line = [self.diamond, 'blastp', '--db', database, '-q', query_file, '-o', blast_file, '-k', max_target, '-e', e_value]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f'Error running Diamond blastp: {e}')

    def run_anchorwave_pro(self, input_file, output_file, r_value, q_value, maximum_gap_size, delete_tandem, tandem_dis):
        command_line = [self.AnchorWave, 'pro', '-i', input_file, '-o', output_file, '-R', r_value, '-Q', q_value, '-D', maximum_gap_size,
                        '-m', delete_tandem, '-W', tandem_dis]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError:
            print("class_Collinearity's function run_AnchorWave_pro failed:")

    def run_all_processes(self):
        self.run_gff_read_get_protein(self.query_genome_seq, self.query_gff_file, self.output_query_pep_seq)
        self.run_gff_read_get_protein(self.ref_genome_seq, self.ref_gff_file, self.output_ref_pep_seq)
        longestPeps.longestPeps(self.query_gff_file, self.query_genome_seq, self.output_query_pep_seq, self.out_query_longest_pep_name)
        longestPeps.longestPeps(self.ref_gff_file, self.ref_genome_seq, self.output_ref_pep_seq, self.out_ref_longest_pep_name)

        self.diamond_make_db(self.out_ref_longest_pep_name, self.database_name)
        self.run_diamond_blastp(self.database_name, self.out_query_longest_pep_name, self.output_blast_result, self.max_target_seqs, self.evalue)

        combineBlastAndStrandInformation.anchorwave_quota(self.ref_gff_file, self.query_gff_file, self.output_blast_result, self.input_file_name)
        self.run_anchorwave_pro(self.input_file_name, self.output_coll_name, self.R, self.Q, self.maximum_gap_size, self.delete_tandem, self.tandem_dis)
        return self.output_coll_name
