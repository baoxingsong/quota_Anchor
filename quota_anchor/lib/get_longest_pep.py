import subprocess
import logging
from . import longestPeps
from multiprocessing import Pool, cpu_count
import os
import sys
from . import base

logger = logging.getLogger('main.get_longest_pep')

class Longest:
    def __init__(self, config_pra, config_soft, parameter):
        
        self.thread = 1
        self.genome_file = ""
        self.gff_file = ""
        self.out_pep_file = ""
        self.out_longest_pep_file = ""
        self.merge = None
        self.overwrite = False

        self.gffread = config_soft['software']['gffread']
        for i in config_pra.sections():
            if i == 'longest_pep':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])

        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.thread = int(self.thread)


    @staticmethod
    def pep_file_empty(file_path):
        try:
            file_path = os.path.abspath(file_path)
            if os.path.exists(file_path):
                if os.path.getsize(file_path) > 0:
                    pass
                else:
                    raise base.FileEmptyError
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file."
            logger.error(exist_error_message)
            sys.exit(1)
        except base.FileEmptyError:
            empty_error_message = "{0} is empty. And gff file and genome file may be don't match.".format(file_path)
            logger.error(empty_error_message)
            sys.exit(1)
        except OSError:
            ose_error_message = f"{file_path} does not exist or is not accessible."
            logger.error(ose_error_message)
            sys.exit(1)

    def run_gffread_get_protein(self, fasta, gff, output_protein_file):
        command_line = [self.gffread, '-g', fasta, '-y', output_protein_file, gff, '-S']
        try:
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            stderr_gff_read = result.stderr
            stdout_gff_read = result.stdout
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
        except subprocess.CalledProcessError as e:
            logger.error(f"Generate {output_protein_file} failed by gffread!")
            error_message = e.stderr
            base.output_info(error_message)

            output_message = e.stdout
            base.output_info(output_message)
            sys.exit(1)
    
    @staticmethod
    def merge_cds_pep(cds_pep_list, merge_file):
        output_file = open(merge_file, 'w')
        for i in range(len(cds_pep_list)):
            file_handler = open(cds_pep_list[i], 'r')
            content = file_handler.read()
            output_file.write(content)
            file_handler.close()
        output_file.close()
    # Todo: merge gffread and longest
    # def gff_read_and_longest(self, genome_file, gff_file, pep_file, longest_pep_file):
    #     self.run_gffread_get_protein(genome_file, gff_file, pep_file)
    #     self.pep_file_empty(pep_file)
    #     longestPeps.longestPeps(gff_file, genome_file, pep_file, longest_pep_file)
    #     self.pep_file_empty(longest_pep_file)

    def sub_run(self, sub_run_include_species, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx):
        pool = Pool(sub_run_include_species)
        for i in range(idx, idx + sub_run_include_species):
            logger.info(f"Generate {out_pep_list[i]} start.")
            self.run_gffread_get_protein(genome_list[i], gff_list[i], out_pep_list[i])
            self.pep_file_empty(out_pep_list[i])
            logger.info(f"Generate {out_pep_list[i]} done!")
            pool.apply_async(longestPeps.longestPeps, args=(gff_list[i], genome_list[i], out_pep_list[i], out_longest_pep_name_list[i]))
        pool.close()
        pool.join()

    def split_para(self):
        genome_list = base.split_conf(self.genome_file, ",")
        for i in genome_list:
            base.file_empty(i)
        gff_list = base.split_conf(self.gff_file, ",")
        for i in gff_list:
            base.file_empty(i)
        out_pep_list = base.split_conf(self.out_pep_file, ",")
        for i in out_pep_list:
            base.output_file_parentdir_exist(i, self.overwrite)
        out_longest_pep_name_list = base.split_conf(self.out_longest_pep_file, ",")
        for i in out_longest_pep_name_list:
            base.output_file_parentdir_exist(i, self.overwrite)
        try:
            assert len(genome_list) == len(gff_list) == len(out_pep_list) == len(out_longest_pep_name_list)
        except AssertionError as e:
            logger.info(f"AssertionError: {e}, please check your separator of variables!")
        return genome_list, gff_list, out_pep_list, out_longest_pep_name_list

    def loop_set(self, gff_list):
        ideal_process_number = len(gff_list)
        cpu_number = cpu_count()
        sub_run_include_species = min(ideal_process_number, cpu_number, self.thread)
        loop_times = int(ideal_process_number / sub_run_include_species)
        final_pool = ideal_process_number % sub_run_include_species
        return loop_times, sub_run_include_species, final_pool, ideal_process_number

    def run_all_process(self):
        logger.info("Longest_pep module init and the following parameters are config information.")
        print()
        for key, value in vars(self).items():
            if key != "gffread" and key != "conf" and value is not None:
                print(key, "=", value)
        print()
        if not self.genome_file:
            logger.error("Please specify your genome file path")
            sys.exit(1)
        if not self.gff_file:
            logger.error("Please specify your gff file path")
            sys.exit(1)
        if not self.out_pep_file:
            logger.error("Please specify your output pep file path")
            sys.exit(1)
        if not self.out_longest_pep_file:
            logger.error("Please specify your longest pep file path")
            sys.exit(1)

        genome_list, gff_list, out_pep_list, out_longest_pep_name_list = self.split_para()
        if hasattr(self, 'merge') and self.merge is not None:
            base.output_file_parentdir_exist(self.merge, self.overwrite)
        loop_times, sub_run_include_species, final_pool, ideal_process_number = self.loop_set(gff_list)

        idx = 0
        for _ in range(loop_times):
            self.sub_run(sub_run_include_species, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx)
            idx += sub_run_include_species
        if final_pool:
            self.sub_run(final_pool, genome_list, gff_list, out_pep_list, out_longest_pep_name_list, idx)

        if hasattr(self, 'merge') and self.merge is not None:
            self.merge_cds_pep(out_longest_pep_name_list, self.merge)
            base.file_empty(self.merge)
        logger.info("Generate species longest protein sequence file finished!")
