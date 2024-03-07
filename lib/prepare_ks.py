import subprocess
from lib import longestCds


# def merge_cds_pep(cds_pep1, cds_pep2, merge_file):
#
#     try:
#         with open(merge_file, 'w') as output_file:
#             subprocess.run(["cat", cds_pep1, cds_pep2], stdout=output_file, check=True)
#     except subprocess.CalledProcessError as e:
#         print("script ks' s function cat_cds_pep failed: ", e)


def merge_cds_pep(cds_pep1, cds_pep2, merge_file):
    with open(cds_pep1, 'r') as f1, open(cds_pep2, 'r') as f2:
        content1 = f1.read()
        content2 = f2.read()
    with open(merge_file, 'w') as output_file:
        output_file.write(content1 + content2)


class Prepare:
    def __init__(self, config_pra, config_soft):

        self.gffread = config_soft['software']['gffread']

        self.query_genome_seq = config_pra['gffread']['query_genome_seq']
        self.query_gff_file = config_pra['gffread']['query_gff_file']
        self.output_query_cds_seq = config_pra['gffread']['output_query_cds_seq']
        self.ref_genome_seq = config_pra['gffread']['ref_genome_seq']
        self.ref_gff_file = config_pra['gffread']['ref_gff_file']
        self.output_ref_cds_seq = config_pra['gffread']['output_ref_cds_seq']

        self.raw_query_prot = config_pra['longestcds']['raw_query_prot']
        self.raw_ref_prot = config_pra['longestcds']['raw_ref_prot']
        self.out_query_cds = config_pra['longestcds']['out_query_cds']
        self.out_ref_cds = config_pra['longestcds']['out_ref_cds']

        self.cds1 = config_pra['ks']['cds1']
        self.cds2 = config_pra['ks']['cds2']
        self.pep1 = config_pra['ks']['pep1']
        self.pep2 = config_pra['ks']['pep2']
        self.cds_file = config_pra['ks']['cds_file']
        self.pep_file = config_pra['ks']['pep_file']

    def run_gff_read_get_cds(self, fasta, gff, output_cds_file):
        command_line = [self.gffread, '-g', fasta, '-x', output_cds_file, gff]
        try:
            subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print("class_Ks's function run_gff_read_get_cds failed: ", e)

    def run(self):
        self.run_gff_read_get_cds(self.query_genome_seq, self.query_gff_file, self.output_query_cds_seq)
        self.run_gff_read_get_cds(self.ref_genome_seq, self.ref_gff_file, self.output_ref_cds_seq)
        longestCds.longest_cds(self.query_gff_file, self.query_genome_seq, self.raw_query_prot, self.output_query_cds_seq, self.out_query_cds)
        longestCds.longest_cds(self.ref_gff_file, self.ref_genome_seq, self.raw_ref_prot, self.output_ref_cds_seq, self.out_ref_cds)
        merge_cds_pep(self.cds1, self.cds2, self.cds_file)
        merge_cds_pep(self.pep1, self.pep2, self.pep_file)
