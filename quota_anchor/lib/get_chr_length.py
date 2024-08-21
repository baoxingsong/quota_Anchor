import pandas as pd
from . import base
from . import GffFile


class Lens:
    def __init__(self, config_pra):
        self.fai_file = config_pra['length']['fai_file']
        self.select_fai_chr_startswith = config_pra['length']['select_fai_chr_startswith']
        self.gff_file = config_pra['length']['gff_file']
        self.length_file = config_pra['length']['length_file']

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
                if i == "number":
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
        new_length_file = base.split_conf(self.length_file, ",")
        gff_file_list = base.split_conf(self.gff_file, ",")

        zipped_three_pair = list(zip(new_fai_file_list, new_start_with, new_length_file, gff_file_list))
        for fai, start, length, gff in zipped_three_pair:
            self.read_fai_gff(start, fai, length, gff)
