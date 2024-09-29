import pandas as pd
from . import base
from . import GffFile
import os
import sys
from . import base

class Lens:
    def __init__(self, config_pra, parameter):
        self.overwrite = False
        for i in config_pra.sections():
            if i == 'length':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)
        print()
    
    @staticmethod
    def file_empty(file_path):
        try:
            file_path = os.path.abspath(file_path)
            if os.path.exists(file_path):
                if os.path.getsize(file_path) > 0:
                    return True
                else:
                    raise base.FileEmptyError
            else:
                raise FileNotFoundError
        except FileNotFoundError as e1:
            exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file"
            print(exist_error_message)  
            return False
        except base.FileEmptyError as e2:
            empty_error_message = "{0} is empty".format(file_path)
            print(empty_error_message)
            return False
        except OSError as e3:
            OSE_error_message = f"{file_path} does not exist or is inaccessible."
            print(OSE_error_message)
            return False


    @staticmethod
    def output_file_parentdir_exist(path, overwite):
        if os.path.exists(path):
            print(f"{path} already exist.")
            if overwite:
                print(f"{path} will be overwrited.")
                return True
            else:
                print(f"{path} will not be overwrited.")
                return False
        path = os.path.abspath(path)
        dir_name = os.path.dirname(path)
        if os.path.isdir(dir_name):
            pass
        else:
            print(f"{dir_name} does not exist and software will make directorys.")
            os.makedirs(dir_name, exist_ok=True)
        return True

    @staticmethod
    def read_fai_gff(selected_prefix, fai_file, output, gff):
        regex = ""
        df = pd.read_csv(fai_file, sep="\t", header=None, index_col=None)
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(int)
        df = df.iloc[:, :2]
        start_list = [i.strip() for i in selected_prefix.split(',')]
        for i in start_list:
            if len(regex) == 0:
                if i == "0-9":
                    regex = r"^\d"
                    continue
                regex = "^{}".format(i)
            else:
                if i == "number":
                    regex += r"|^\d"
                    continue
                regex += "|^{}".format(i)
        lens = df[df[0].str.match(regex)]
        lens.reset_index(drop=True, inplace=True)
        chr_list = lens.iloc[:, 0].copy()
        _, ref_chr_gene_list, _, _ = GffFile.readGff(gff)
        new_column = []
        for ch in chr_list:
            total_gene = len(ref_chr_gene_list[ch])
            new_column.append(total_gene)

        total_gene_series = pd.DataFrame(new_column)

        new_lens = pd.concat([lens, total_gene_series], axis=1)
        new_lens.columns = ["chr", "length", "total_gene"]
        new_lens.to_csv(output, sep='\t', index=False, header=True)

    def run(self):
        new_fai_file_list = base.split_conf(self.fai_file, ",")
        new_start_with = base.split_conf(self.select_fai_chr_startswith, ":")
        new_length_file = base.split_conf(self.output_length_file, ",")
        gff_file_list = base.split_conf(self.gff_file, ",")

        zipped_three_pair = list(zip(new_fai_file_list, new_start_with, new_length_file, gff_file_list))
        for fai, start, length, gff in zipped_three_pair:
            flag1 = self.file_empty(fai)
            if not flag1:
                continue
            flag2 = self.file_empty(gff)
            if not flag2:
                continue
            flag3 = self.output_file_parentdir_exist(length, self.overwrite)
            if not flag3:
                continue
            self.read_fai_gff(start, fai, length, gff)
