import logging
import sys
from . import base
import subprocess    
from . import combineBlastAndStrandInformation

logger = logging.getLogger('main.pre_collinearity')

class Prepare:
    def __init__(self, config_pra, config_soft, parameter):
        self.overwrite = False
        self.strand = "plus"
        self.skip_blast = parameter.skip_blast
        self.bitscore = 100
        self.align_length = 0
        self.outfmt = 6
        self.evalue = "1e-10"
        self.thread = 6
        self.dtype = "prot"
        # blank string
        # attributes = [
        #     'ref_gff_file', 'query_gff_file', 'ref_seq', 'query_seq', 'blast_file',
        #     'database_name', 'output_blast_result', 'max_target_seqs', 'output_file'
        #     'diamond', 'blastn', 'blastp', 'makeblastdb'
        # ]
        # for attr in attributes:
        #     setattr(self, attr, "")
        self.ref_gff_file, self.query_gff_file, self.ref_seq, self.query_seq, self.blast_file = "", "", "", "", ""
        self.database_name, self.output_blast_result, self.max_target_seqs, self.output_file = "", "", "", ""
        self.diamond, self.blastp, self.blastn, self.makeblastdb = "", "", "", ""

        # config file
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
        software_list = []
        for key in config_soft['software']:
            software_list.append(key)
            setattr(self, key, config_soft["software"][key])
        self.software_list = software_list
        # ref query length file
        if not hasattr(self, "ref_length"):
            self.ref_length = ""
        if not hasattr(self, "query_length"):
            self.query_length = ""
        # convert str
        self.max_target_seqs = str(self.max_target_seqs)
        self.outfmt = str(self.outfmt)
        self.thread = str(self.thread)

    def diamond_makedb(self, ref_protein, database):
        command_line = [self.diamond, 'makedb', '--in', ref_protein, '--db', database]
        try:
            logger.info(f"Run diamond makedb and generate {database} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,  text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run diamond makedb and generate {database} end.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Run diamond makedb and generate {database} failed!")

            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)

    def run_diamond_blastp(self, database, query_file, blast_file, max_target, e_value):
        command_line = [self.diamond, 'blastp', '--db', database, '-q', query_file, '-o', blast_file, '-k', max_target, '-e', e_value]
        try:
            logger.info(f"Run diamond blastp and generate {blast_file} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,  text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run diamond blastp and generate {blast_file} end.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Run diamond blastp and generate {blast_file} failed!")

            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)
    
    def mkblastdb(self, ref_seq, database, dtype):
        command_line = [self.makeblastdb, '-in', ref_seq, '-dbtype', dtype, '-out', database]
        try:
            logger.info(f"Run makeblastdb and generate {database} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,  text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run makeblastdb and generate {database} end.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Run makeblastdb and generate {database} failed!")

            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)

    def run_blastn(self, blast_database, query_file, blast_file, e_value, thread, outfmt, max_target_seqs, strand):
        command_line = [self.blastn, '-query', query_file, '-out', blast_file, '-evalue', e_value,
                        '-db', blast_database, '-num_threads', thread, '-max_target_seqs', max_target_seqs, '-outfmt', outfmt, '-strand', strand]
        try:
            logger.info(f"Run blastn and generate {blast_file} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,  text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run blastn and generate {blast_file} end.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Run blastn and generate {blast_file} failed!")

            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)

    def run_blastp(self, blast_database, query_file, blast_file, e_value, thread, outfmt, max_target_seqs):
        command_line = [self.blastp, '-query', query_file, '-out', blast_file, '-evalue', e_value,
                        '-db', blast_database, '-num_threads', thread, '-max_target_seqs', max_target_seqs, '-outfmt', outfmt]
        try:
            logger.info(f"Run blastp and generate {blast_file} start.")
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,  text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
            logger.info(f"Run blastp and generate {blast_file} end.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Run blastp and generate {blast_file} failed!")

            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)

    def print_config(self):
        print()
        for key, value in vars(self).items():
            if key not in self.software_list and key != "conf" and key != "software_list" and key != "output_blast_result":
                if not self.skip_blast:
                    try:
                        if self.align == "diamond":
                            key_list = ["thread", "outfmt", "dtype", "strand", "skip_blast"]
                            if key not in key_list:
                                print(key, "=", value)
                        if self.align == "blastp" or self.align == "blastn":
                            if key != "strand" and key != "skip_blast":
                                print(key, "=", value)
                    except AttributeError:
                        logger.error(r'You need to specify the --skip_blast option if you want to skip the blast step. Otherwise, specify --align diamond/blastp/blastn')
                        sys.exit(1)
                else:
                    blast_key_list = ["align", "database_name", "output_blast_result", "strand",
                                      "ref_seq", "query_seq", "max_target_seqs", "evalue", "thread", "outfmt", "dtype"]
                    if key not in blast_key_list:
                        print(key, "=", value)
        print()
    def pre_col_init(self):
        if not self.skip_blast:
            if not self.ref_seq:
                logger.error("Please specify your reference species protein/cds file")
                sys.exit(1)
            if not self.query_seq:
                logger.error("Please specify your query species protein/cds file")
                sys.exit(1)
        if not self.ref_gff_file:
            logger.error("Please specify your reference species gff file")
            sys.exit(1)
        if not self.query_gff_file:
            logger.error("Please specify your query species gff file")
            sys.exit(1)
        if not self.output_file:
            logger.error("Please specify your output table file name (called blast file containing gene position information)")
            sys.exit(1)
        if not self.blast_file:
            logger.error("Please specify your blast file name")
            sys.exit(1)

    def run_all_process(self):
        logger.info("Pre_collinearity module init and the following parameters are config information.")
        self.print_config()
        self.pre_col_init()
        if self.ref_length:
            base.file_empty(self.ref_length)
        if self.query_length:
            base.file_empty(self.query_length)

        base.file_empty(self.ref_gff_file)
        base.file_empty(self.query_gff_file)
        base.output_file_parentdir_exist(self.output_file, self.overwrite)

        self.output_blast_result = self.blast_file
        if not self.skip_blast:
            base.file_empty(self.ref_seq)
            base.file_empty(self.query_seq)
            
            if self.align == "diamond":
                self.diamond_makedb(self.ref_seq, self.database_name)
                self.run_diamond_blastp(self.database_name, self.query_seq, self.output_blast_result,
                                        self.max_target_seqs, self.evalue)
            if self.align == "blastn":
                self.mkblastdb(self.ref_seq, self.database_name, self.dtype)
                self.run_blastn(self.database_name, self.query_seq, self.output_blast_result,
                                self.evalue, self.thread, self.outfmt, self.max_target_seqs, self.strand)
            if self.align == "blastp":
                self.mkblastdb(self.ref_seq, self.database_name, self.dtype)
                self.run_blastp(self.database_name, self.query_seq, self.output_blast_result,
                                self.evalue, self.thread, self.outfmt, self.max_target_seqs)
        base.file_empty(self.blast_file)
        combineBlastAndStrandInformation.anchorwave_quota(self.ref_gff_file, self.query_gff_file, self.blast_file,
                                                          self.output_file, self.bitscore, self.align_length, self.query_length, self.ref_length)
        base.file_empty(self.output_file)
        logger.info(f'Generate {self.output_file} finished!')
