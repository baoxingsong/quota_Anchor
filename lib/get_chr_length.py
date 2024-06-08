import pandas as pd


class Lens:
    def __init__(self, config_pra):
        self.ref_fai = config_pra['length']['ref_fai']
        self.query_fai = config_pra['length']['query_fai']
        self.select_ref_fai_chr_startswith = config_pra['length']['select_ref_fai_chr_startswith']
        self.select_query_fai_chr_startswith = config_pra['length']['select_query_fai_chr_startswith']
        self.ref_length = config_pra['length']['ref_length']
        self.query_length = config_pra['length']['query_length']

    @staticmethod
    def read_fai(selected_prefix, fai_file, output):
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
        lens.columns = ["chr", "length"]
        lens.loc[:, 'chr'] = lens.loc[:, 'chr']
        lens.to_csv(output, sep='\t', index=False, header=True)
        return

    def run(self):
        self.read_fai(self.select_ref_fai_chr_startswith, self.ref_fai, self.ref_length)
        self.read_fai(self.select_query_fai_chr_startswith, self.query_fai, self.query_length)
