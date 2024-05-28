#Copyright (c) 2018-2018, Pengchuan Sun
#
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification,
#are permitted provided that the following conditions are met:
#
#Redistributions of source code must retain the above copyright notice, this list
#of conditions and the following disclaimer.
#
#Redistributions in binary form must reproduce the above copyright notice, this
#list of conditions and the following disclaimer in the documentation and/or
#other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, linregress


class PeaksFit:
    def __init__(self, config_par):
        self.histogram_only = 'false'

        self.blockinfo = config_par['peaksfit']['blockinfo']
        self.mode = config_par['peaksfit']['mode']
        self.bins_number = int(config_par["peaksfit"]["bins_number"])
        self.ks_area = config_par["peaksfit"]["ks_area"]
        self.fontsize = config_par["peaksfit"]["fontsize"]
        self.area = config_par["peaksfit"]["area"]
        self.figsize = config_par["peaksfit"]["figsize"]
        self.savefig = config_par["peaksfit"]["savefig"]

        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.ks_area = [float(k) for k in self.ks_area.split(',')]
        self.peaks = 1

    @staticmethod
    def ks_values(df):
        ks = df['ks'].str.split('_')
        ks_total = []
        for v in ks.values:
            ks_total.extend([float(k) for k in v])
        ks_average = df['ks_average'].values
        ks_median = df['ks_median'].values
        return [ks_median, ks_average, ks_total]

    @staticmethod
    def gaussian_fuc(x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            amp = float(params[i])
            ctr = float(params[i+1])
            wid = float(params[i+2])
            y = y + amp * np.exp(-((x - ctr)/wid)**2)
        return y

    def kde_fit(self, data, x):
        kde = gaussian_kde(data)
        kde.set_bandwidth(bw_method=kde.factor/3.)
        p = kde(x)
        guess = [1, 1, 1]*self.peaks
        result = curve_fit(self.gaussian_fuc, x, p, guess, maxfev=80000)
        popt = result[0]
        popt = [abs(k) for k in popt]
        data = []
        y = self.gaussian_fuc(x, *popt)
        for i in range(0, len(popt), 3):
            array = [popt[i], popt[i+1], popt[i+2]]
            data.append(self.gaussian_fuc(x, *array))
        slope, intercept, r_value, p_value, std_err = linregress(p, y)
        print("\nR-square: "+str(r_value**2))
        print("The gaussian fitting curve parameters are :")
        print('  |  '.join([str(k) for k in popt]))
        return y, data

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['mathtext.default'] = 'regular'
        fig, ax = plt.subplots(figsize=self.figsize)
        bkinfo = pd.read_csv(self.blockinfo, low_memory=False, header=0)
        ks_median, ks_average, ks_total = self.ks_values(bkinfo)
        data = eval('ks_'+self.mode)
        data = [k for k in data if self.ks_area[0] <= k <= self.ks_area[1]]
        x = np.linspace(self.ks_area[0], self.ks_area[1], self.bins_number)
        ax.hist(data, int(
            self.bins_number), density=1, facecolor='blue', alpha=0.3, label='Histogram')
        if self.histogram_only.upper() == 'TRUE':
            pass
        else:
            y, fit = self.kde_fit(data, x)
            ax.plot(x, y, color='black', linestyle='-', label='Gaussian fitting')
        ax.grid()
        ax.set_xlabel(r'${K_s}$', fontsize=20)
        ax.set_ylabel('Frequency', fontsize=20)
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=20)
        ax.set_xlim(self.ks_area)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.12)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
