import pandas as pd
import matplotlib.pyplot as plt


class ClsVis:
    def __init__(self, config_par):
        self.pair_file = config_par['class_vil']['stats_file']
        self.type = config_par['class_vil']['type']
        self.output = config_par['class_vil']['figure']
        self.title = config_par['class_vil']['title']
        self.ylab = config_par['class_vil']['ylab']

    def get_data(self, pair_file):
        df = pd.read_csv(pair_file, sep="\t", index_col="Type", header=0)
        if self.type == "pair_number" or self.type == "pair" or self.type.startswith("p"):
            df = df.loc[["wgd.pairs", "tandem.pairs", "proximal.pairs", "transposed.pairs", "dispersed.pairs"], :]
        if self.type == "gene_number" or self.type == "gene" or self.type.startswith("g"):
            df = df.loc[["wgd.genes", "tandem.genes", "proximal.genes", "transposed.genes", "dispersed.genes", "singleton.genes"], :]
        return df

    def plot(self, df):
        fig, ax = plt.subplots()
        tp = [ele for ele in df.index]
        length = len(tp)
        counts = [df.loc[tp[i], "Number"] for i in range(length)]
        color_list = ['#FFB6C1', '#C71585', '#FF00FF', '#9932CC', '#00CED1', "#90EE90"]
        # bar_labels = ['red', 'blue', '_red', 'orange']
        # bar_colors = ['tab:red', 'tab:blue', 'tab:red', 'tab:orange']

        ax.bar(tp, counts, label=tp, color=color_list[:length])
        if self.ylab == "" or self.ylab == "None" or self.ylab == "none":
            print(123)
            self.ylab = "Number"
        ax.set_ylabel(str(self.ylab))
        if self.title == "" or self.title == "None" or self.title == "none":
            self.title = "Different gene types of species"
        ax.set_title(str(self.title))
        ax.legend(title='Gene Type')
        plt.savefig(self.output)
        plt.show()

    def run(self):
        df = self.get_data(self.pair_file)
        self.plot(df)
