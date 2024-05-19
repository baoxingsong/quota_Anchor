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
from scipy.stats.kde import gaussian_kde


class KsPeaks:
    def __init__(self, config_par):
        self.tandem_length = 200

        self.block_length = config_par['kspeaks']['block_length']
        self.blockinfo = config_par['kspeaks']['blockinfo']
        self.tandem = config_par["kspeaks"]["tandem"]
        self.ks_area = config_par["kspeaks"]["ks_area"]
        self.multiple = config_par["kspeaks"]["multiple"]
        self.homo = config_par["kspeaks"]["homo"]
        self.fontsize = config_par["kspeaks"]["fontsize"]
        self.area = config_par["kspeaks"]["area"]
        self.figsize = config_par["kspeaks"]["figsize"]
        self.savefig = config_par["kspeaks"]["savefig"]
        self.savefile = config_par["kspeaks"]["savefile"]

        self.homo = [float(k) for k in self.homo.split(',')]
        self.ks_area = [float(k) for k in self.ks_area.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1']-group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1']-group.loc[:, 'end2']
        index = group[(group['start'].abs() <= int(self.tandem_length)) | (
            group['end'].abs() <= int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    @staticmethod
    def ks_kde(df):
        ks = df['ks'].str.split('_')
        arr = []
        for v in ks.values:
            v = [float(k) for k in v if float(k) >= 0]
            if len(v) == 0:
                continue
            arr.extend(v)
        kdemedian = gaussian_kde(df['ks_median'].values)
        kdemedian.set_bandwidth(bw_method=kdemedian.factor / 3.)
        kdeaverage = gaussian_kde(df['ks_average'])
        kdeaverage.set_bandwidth(bw_method=kdeaverage.factor / 3.)
        kdetotal = gaussian_kde(arr)
        kdetotal.set_bandwidth(bw_method=kdetotal.factor / 3.)
        return [kdemedian, kdeaverage, kdetotal]

    def run(self):
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['mathtext.default'] = 'regular'
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        bkinfo = pd.read_csv(self.blockinfo, low_memory=False, header=0)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[bkinfo['length'] > int(self.block_length)]
        if self.tandem.upper() == 'FALSE':
            bkinfo = self.remove_tandem(bkinfo)
        bkinfo = bkinfo[bkinfo['length'] >= int(self.block_length)]
        bkinfo = bkinfo[bkinfo['homo' + str(self.multiple)] >= float(self.homo[0])]
        bkinfo = bkinfo[bkinfo['homo' + str(self.multiple)] <= float(self.homo[1])]
        bkinfo = bkinfo[bkinfo['ks_median'] >= float(self.ks_area[0])]
        bkinfo = bkinfo[bkinfo['ks_median'] <= float(self.ks_area[1])]
        kdemedian, kdeaverage, kdetotal = self.ks_kde(bkinfo)
        dist_space = np.linspace(self.area[0], self.area[1], 500)
        ax.plot(dist_space, kdemedian(dist_space), color='red', label='block median')
        ax.plot(dist_space, kdeaverage(dist_space), color='black', label='block average')
        ax.plot(dist_space, kdetotal(dist_space), color='blue', label='all pairs')
        ax.grid()
        ax.set_xlabel(r'${K_s}$', fontsize=20)
        ax.set_ylabel('Frequency', fontsize=20)
        ax.tick_params(labelsize=18)
        ax.set_xlim(self.area)
        ax.legend(fontsize=20)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.12)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        bkinfo.to_csv(self.savefile, index=False)
        sys.exit(0)
