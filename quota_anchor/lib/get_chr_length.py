import sys
import pandas as pd
from . import GffFile
import os
from . import base
import logging

logger = logging.getLogger('main.get_chr_length')
class Lens:
    def __init__(self, config_pra, parameter):
        self.fai_file = ""
        self.select_fai_chr_startswith = "chr,Chr,CHR,0-9"
        self.output_length_file = ""
        self.gff_file = ""
        self.overwrite = False
        for i in config_pra.sections():
            if i == 'length':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])
        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)
    
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
        except FileNotFoundError:
            exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file."
            logger.error(exist_error_message)
            return False
        except base.FileEmptyError:
            empty_error_message = "{0} is empty.".format(file_path)
            logger.error(empty_error_message)
            return False
        except OSError:
            ose_error_message = f"{file_path} does not exist or is not accessible."
            logger.error(ose_error_message)
            return False


    @staticmethod
    def output_file_parentdir_exist(path, overwrite):
        if os.path.exists(path):
            logger.info(f"{path} already exist.")
            if overwrite:
                logger.info(f"{path} will be overwrote.")
                return True
            else:
                logger.info(f"{path} will not be overwrote.")
                return False
        path = os.path.abspath(path)
        dir_name = os.path.dirname(path)
        if os.path.isdir(dir_name):
            pass
        else:
            logger.info(f"{dir_name} does not exist and software will make directories.")
            os.makedirs(dir_name, exist_ok=True)
        return True

    @staticmethod
    def read_fai_gff(selected_prefix, fai_file, output, gff):
        logger.info(f'generate {output} start.')
        regex = ""
        df = pd.read_csv(fai_file, sep="\t", header=None, index_col=None)
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(int)
        df = df.iloc[:, :2]
        start_list = [i.strip() for i in selected_prefix.split(',')]
        for i in start_list:
            if len(i) == 0:
                continue   
            else: 
                if len(regex) == 0:
                    if i == "0-9":
                        regex = r"^\d"
                        continue
                    regex = "^{}".format(i)
                else:
                    if i == "0-9":
                        regex += r"|^\d"
                        continue
                    regex += "|^{}".format(i)
        lens = df[df[0].str.match(regex)]
        if len(lens) == 0:
            logger.error(f'{fai_file} don\'t have a chromosome name starts with regex:{regex}.')
            logger.info(f'skip generate {output}')
            return
        lens.reset_index(drop=True, inplace=True)
        chr_list = lens.iloc[:, 0].copy()
        _, ref_chr_gene_list, _, _ = GffFile.readGff(gff)
        new_column = []
        for ch in chr_list:
            try:
                total_gene = len(ref_chr_gene_list[ch])
            except KeyError:
                total_gene = 0
            new_column.append(total_gene)

        total_gene_series = pd.DataFrame(new_column)
        new_lens = pd.concat([lens, total_gene_series], axis=1)
        new_lens.columns = ["chr", "length", "total_gene"]
        new_lens.to_csv(output, sep='\t', index=False, header=True)
        logger.info(f'generate {output} end!')

    def run(self):
        logger.info("Init get_chr_length and the following parameters are config information.")

        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)
        print()

        new_fai_file_list = base.split_conf(self.fai_file, ",")
        new_start_with = base.split_conf(self.select_fai_chr_startswith, ":")
        new_length_file = base.split_conf(self.output_length_file, ",")
        gff_file_list = base.split_conf(self.gff_file, ",")
        try:
            assert len(new_fai_file_list) == len(new_start_with) == len(new_length_file) == len(gff_file_list)
        except AttributeError:
            logger.error('please check your separator for four variables.')
            sys.exit(1)
        zipped_three_pair = list(zip(new_fai_file_list, new_start_with, new_length_file, gff_file_list))
        # skip trios pair which has some problems
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
        logger.info("get chromosome information finished!")

