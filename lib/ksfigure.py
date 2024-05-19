import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


class KsFigure:
    def __init__(self, config_par):
        self.ksfit = config_par['ksfigure']['ksfit']
        self.labelfontsize = config_par['ksfigure']['labelfontsize']
        self.legendfontsize = config_par['ksfigure']['legendfontsize']
        self.xlabel = config_par["ksfigure"]["xlabel"]
        self.ylabel = config_par["ksfigure"]["ylabel"]
        self.title = config_par["ksfigure"]["title"]
        self.area = config_par["ksfigure"]["area"]
        self.figsize = config_par["ksfigure"]["figsize"]
        self.shadow = config_par["ksfigure"]["shadow"]
        self.savefig = config_par["ksfigure"]["savefig"]
        if self.xlabel == 'none' or self.xlabel == '':
            self.xlabel = r'Synonymous nucleotide substitution (${K_{s}}$)'
        if self.ylabel == 'none' or self.ylabel == '':
            self.ylabel = 'kernel density of syntenic blocks'
        if self.title == 'none' or self.title == '':
            self.title = ''
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]

    @staticmethod
    def gaussian_distribution(t, k):
        y = np.zeros(len(t))
        for i in range(0, int((len(k) - 1) / 3)+1):
            if np.isnan(k[3 * i + 2]):
                continue
# 1 / (np.sqrt(2 * np.pi) * sigma )
# ctr u
# np.sqrt(2) * sigma
            k[3 * i + 2] = float(k[3 * i + 2])/np.sqrt(2)              # almost sigma
            k[3 * i + 0] = float(k[3 * i + 0]) * \
                np.sqrt(2*np.pi)*float(k[3 * i + 2])                     # almost equal to 1
            y1 = stats.norm.pdf(
                t, float(k[3 * i + 1]), float(k[3 * i + 2])) * float(k[3 * i + 0])
            y = y+y1
        return y

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['mathtext.default'] = 'regular'
        fig, ax = plt.subplots(figsize=self.figsize)
        ksfit = pd.read_csv(self.ksfit, index_col=0, low_memory=False)
        t = np.arange(self.area[0], self.area[1], 0.0005)
        col = [k for k in ksfit.columns if re.match('Unnamed:', k)]
        for index, row in ksfit.iterrows():
            ax.plot(t, self.gaussian_distribution(
                t, row[col].values), linestyle=row['linestyle'], color=row['color'], alpha=0.8, label=index, linewidth=row['linewidth'])
            if self.shadow.upper() == 'TRUE' or self.shadow.upper() == 'T':
                ax.fill_between(t, 0, self.gaussian_distribution(t, row[col].values),  color=row['color'], alpha=0.15)
        align = dict(verticalalignment="center",
                     horizontalalignment="center")
        ax.set_xlabel(self.xlabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_ylabel(self.ylabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_title(self.title, weight='bold',
                     fontsize=self.labelfontsize, **align)
        plt.tick_params(labelsize=10)
        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles=handles[:int(len(labels))], labels=labels[:int(len(labels))], loc='upper right', prop={
               'family': 'Times New Roman', 'style': 'italic', 'size': self.legendfontsize})
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
