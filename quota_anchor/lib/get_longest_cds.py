import subprocess
import logging
from . import longestCds
from multiprocessing import Pool, cpu_count
import os
import sys
from . import base

logger = logging.getLogger('main.get_longest_cds')

class Longest:
    def __init__(self, config_pra, config_soft, parameter):
        
        self.thread = 1
        self.overwrite = False
        self.genome_file = ""
        self.gff_file = ""
        self.out_cds_file = ""
        self.out_longest_cds_file = ""
        self.merge = None
        self.gffread = config_soft['software']['gffread']
        for i in config_pra.sections():
            if i == 'longest_cds':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

        self.thread = int(self.thread)

    @staticmethod
    def merge_cds_cds(cds_cds_list, merge_file):
        output_file = open(merge_file, 'w')
        for i in range(len(cds_cds_list)):
            file_handler = open(cds_cds_list[i], 'r')
            content = file_handler.read()
            output_file.write(content)
            file_handler.close()
        output_file.close()

    @staticmethod
    def cds_file_empty(file_path):
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
            empty_error_message = "{0} is empty. gff file and genome file don't match.".format(file_path)
            logger.error(empty_error_message)
            sys.exit(1)
        except OSError:
            ose_error_message = f"{file_path} does not exist or is not accessible."
            logger.error(ose_error_message)
            sys.exit(1)

    def run_gffread_get_cds(self, fasta, gff, output_cds_file):
        command_line = [self.gffread, '-g', fasta, '-x', output_cds_file, gff]
        try:
            result = subprocess.run(command_line, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_gff_read = result.stderr.decode()
            stdout_gff_read = result.stdout.decode()
            base.output_info(stderr_gff_read)
            base.output_info(stdout_gff_read)
        except subprocess.CalledProcessError:
            logger.error(f"generate {output_cds_file} failed by gffread!")
            sys.exit(1)

    def sub_run(self, sub_run_number, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx):
        pool = Pool(sub_run_number)
        for i in range(idx, idx + sub_run_number):
            logger.info(f"generate {out_cds_list[i]} start.")
            self.run_gffread_get_cds(genome_list[i], gff_list[i], out_cds_list[i])
            self.cds_file_empty(out_cds_list[i])
            logger.info(f"generate {out_cds_list[i]} done!")
            pool.apply_async(longestCds.longest_cds, args=(gff_list[i], genome_list[i], out_cds_list[i], out_longest_cds_name_list[i]))
        pool.close()
        pool.join()

    def split_para(self):
        genome_list = base.split_conf(self.genome_file, ",")
        for i in genome_list:
            base.file_empty(i)
        gff_list = base.split_conf(self.gff_file, ",")
        for i in gff_list:
            base.file_empty(i)
        out_cds_list = base.split_conf(self.out_cds_file, ",")
        for i in out_cds_list:
            base.output_file_parentdir_exist(i, self.overwrite)
        out_longest_cds_name_list = base.split_conf(self.out_longest_cds_file, ",")
        for i in out_longest_cds_name_list:
            base.output_file_parentdir_exist(i, self.overwrite)
        try:
            assert len(genome_list) == len(gff_list) == len(out_cds_list) == len(out_longest_cds_name_list)
        except AssertionError as e:
            logger.info(f"AssertionError: {e}, please check your separator of variables!")
        return genome_list, gff_list, out_cds_list, out_longest_cds_name_list

    def loop_set(self, gff_list):
        ideal_process_number = len(gff_list)
        cpu_number = cpu_count()
        sub_run_include_species = min(ideal_process_number, cpu_number, self.thread)
        loop_times = int(ideal_process_number / sub_run_include_species)
        final_pool = ideal_process_number % sub_run_include_species
        return loop_times, sub_run_include_species, final_pool, ideal_process_number

    def run_all_process(self):
        logger.info("Init longest_cds and the following parameters are config information.")
        print()
        for key, value in vars(self).items():
            if key != "gffread" and key != "conf" and value is not None:
                print(key, "=", value)
        print()

        genome_list, gff_list, out_cds_list, out_longest_cds_name_list = self.split_para()
        if hasattr(self, 'merge') and self.merge is not None:
            base.output_file_parentdir_exist(self.merge, self.overwrite)
        loop_times, sub_run_include_species, final_pool, ideal_process_number = self.loop_set(gff_list)
        
        idx = 0
        for _ in range(loop_times):
            self.sub_run(sub_run_include_species, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx)
            idx += sub_run_include_species
        if final_pool:
            self.sub_run(final_pool, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx)

        if hasattr(self, 'merge') and self.merge is not None:
            self.merge_cds_cds(out_longest_cds_name_list, self.merge)
            base.file_empty(self.merge)
        logger.info("get longest cds finished!")