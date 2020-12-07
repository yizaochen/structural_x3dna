from os import path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class BasePairAgent:
    hosts = ['a_tract_21mer', 'ctct_21mer', 'gcgc_21mer',
             'g_tract_21mer', 'atat_21mer', 'tgtg_21mer']
    abbr_hosts = {'a_tract_21mer': 'A-tract', 'ctct_21mer': 'CTCT', 'gcgc_21mer': 'GCGC',
                  'g_tract_21mer': 'G-tract', 'atat_21mer': 'ATAT', 'tgtg_21mer': 'TGTG'}
    d_colors = {'a_tract_21mer': 'blue', 'atat_21mer': 'orange', 'ctct_21mer': 'green',
                'g_tract_21mer': 'red', 'gcgc_21mer': 'magenta', 'tgtg_21mer': 'cyan'}
    hosts_group = [['a_tract_21mer', 'g_tract_21mer'],
                   ['atat_21mer', 'gcgc_21mer'],
                   ['ctct_21mer', 'tgtg_21mer']]

    def __init__(self, rootfolder, time_interval):
        self.rootfolder = rootfolder
        self.type_na = 'bdna+bdna'
        self.time_interval = time_interval
        
        self.start_bp = 4
        self.end_bp = 18
        self.n_bp = self.end_bp - self.start_bp + 1

        self.lbfz = 12
        self.lgfz = 12
        self.ticksize = 10

        self.parameters = ['shear', 'buckle', 'stretch', 'propeller', 'stagger', 'opening']

    def histogram_three_groups(self, figsize, parameter, bins, xlines, xlim, ylim):
        nrows = 1
        ncols = 3
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, sharey=True, sharex=True)
        col_id = 0
        xlabel = self.get_xlabel(parameter)
        for hosts in self.hosts_group:
            ax = axes[col_id]
            for host in hosts:
                data = self.get_data(host, parameter)
                ax.hist(data, bins=bins, density=True, color=self.d_colors[host], alpha=0.4, label=self.abbr_hosts[host])
            for xvalue in xlines:
                ax.axvline(xvalue, color='grey', alpha=0.2)
            ax.set_xlabel(xlabel,fontsize=self.lbfz)
            ax.set_ylabel('P', fontsize=self.lbfz)
            ax.legend(frameon=False, fontsize=self.lgfz)
            ax.tick_params(axis='both', labelsize=self.ticksize)
            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
            col_id += 1
        return fig, axes

    def histogram_six_sys_one_para(self, figsize, parameter, bins):
        nrows = 2
        ncols = 3
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, sharey=True, sharex=True)
        d_axes = self.get_daxes(axes, nrows, ncols)
        for host in self.hosts:
            data = self.get_data(host, parameter)
            ax = d_axes[host]
            ax.hist(data, bins=bins, density=True, color=self.d_colors[host], alpha=0.4)
            ax.set_title(self.abbr_hosts[host])
            ax.set_xlabel(parameter)
            ax.set_ylabel('P')
        return fig, axes

    def get_xlabel(self, parameter):
        if parameter in ['shear', 'stretch', 'stagger']:
            return f'{parameter}(Å)'
        else:
            return f'{parameter}(°)'

    def get_daxes(self, axes, nrows, ncols):
        d_axes = dict()
        host_id = 0
        for row_id in range(nrows):
            for col_id in range(ncols):
                host = self.hosts[host_id]
                d_axes[host] = axes[row_id, col_id]
                host_id += 1
        return d_axes

    def get_data(self, host, parameter):
        host_time_folder = path.join(self.rootfolder, host, self.time_interval)
        fname = path.join(host_time_folder, f'{parameter}.csv')
        df = pd.read_csv(fname, index_col='Frame-ID')
        n_frame = df.shape[0]
        temp_array = np.zeros((n_frame, self.n_bp))
        col_id = 0
        for bp_id in range(self.start_bp, self.end_bp+1):
            temp_array[:,col_id] = df[f'bp{bp_id}']
            col_id += 1
        return np.ndarray.flatten(temp_array)

class BaseStepAgent(BasePairAgent):
    def __init__(self, rootfolder, time_interval):
        super().__init__(rootfolder, time_interval)

    def get_xlabel(self, parameter):
        if parameter in ['shift', 'slide', 'rise']:
            return f'{parameter}(Å)'
        else:
            return f'{parameter}(°)'

    def get_data(self, host, parameter):
        host_time_folder = path.join(self.rootfolder, host, self.time_interval)
        fname = path.join(host_time_folder, f'{parameter}.csv')
        df = pd.read_csv(fname, index_col='Frame-ID')
        n_frame = df.shape[0]
        temp_array = np.zeros((n_frame, self.n_bp))
        col_id = 0
        for bp_id in range(self.start_bp, self.end_bp):
            temp_array[:,col_id] = df[f'bp{bp_id}_bp{bp_id+1}']
            col_id += 1
        return np.ndarray.flatten(temp_array)