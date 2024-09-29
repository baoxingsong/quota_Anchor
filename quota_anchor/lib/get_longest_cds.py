import subprocess
from . import longestCds
from multiprocessing import Pool, cpu_count
import os
import sys
from . import base


class Longest:
    def __init__(self, config_pra, config_soft, parameter):
        
        self.thread = 1
        self.overwrite = False
        self.gffread = config_soft['software']['gffread']
        for i in config_pra.sections():
            if i == 'longest_cds':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "gffread" and key != "conf":
                print(key, "=", value)
        print()
        self.thread = int(self.thread)

    @staticmethod
    def merge_cds_pep(cds_pep_list, merge_file):
        output_file = open(merge_file, 'w')
        for i in range(len(cds_pep_list)):
            file_handler = open(cds_pep_list[i], 'r')
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
        except FileNotFoundError as e1:
            exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file"
            print(exist_error_message)  
            sys.exit(1)
        except base.FileEmptyError as e2:
            empty_error_message = "{0} is empty. gff file and genome file don't match.".format(file_path)
            print(empty_error_message)
            sys.exit(1)
        except OSError as e3:
            OSE_error_message = f"{file_path} does not exist or is inaccessible."
            print(OSE_error_message)
            sys.exit(1)

    def run_gffread_get_cds(self, fasta, gff, output_cds_file):
        command_line = [self.gffread, '-g', fasta, '-x', output_cds_file, gff]
        try:
            result = subprocess.run(command_line, check=True)
        except subprocess.CalledProcessError as e:
            print(f"{result}, class Longest's function run_gffread_get_cds failed: ")
    
    def sub_run(self, sub_run_number, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx):

        pool = Pool(sub_run_number)
        for i in range(idx, idx + sub_run_number):
            self.run_gffread_get_cds(genome_list[i], gff_list[i], out_cds_list[i])
            self.cds_file_empty(out_cds_list[i])
            pool.apply_async(longestCds.longest_cds, args=(gff_list[i], genome_list[i], out_cds_list[i], out_longest_cds_name_list[i]))
        pool.close()
        pool.join()
    
    def run_all_process(self):
        
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
            print(f"AssertionError: {e}, please check your separator of variables!")

        ideal_process_number = len(gff_list)
        sub_run_number = min(ideal_process_number, cpu_count(), self.thread)
        interger_pool = int(ideal_process_number / sub_run_number) 
        final_pool = ideal_process_number % sub_run_number
        
        idx = 0
        for _ in range(interger_pool):
            self.sub_run(sub_run_number, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx)
            idx += sub_run_number
        if final_pool:
            self.sub_run(final_pool, genome_list, gff_list, out_cds_list, out_longest_cds_name_list, idx)

        for i in range(ideal_process_number):
            base.file_empty(out_longest_cds_name_list[i])
        if hasattr(self, 'merge') and self.merge != None:
            self.merge_cds_pep(out_longest_cds_name_list, self.merge)
            base.file_empty(self.merge)