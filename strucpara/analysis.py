from os import path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class BasePairAgent:
    hosts = ['a_tract_21mer', 'atat_21mer', 'ctct_21mer',
             'g_tract_21mer', 'gcgc_21mer', 'tgtg_21mer']
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
    d_bimodal_label = {'atat_21mer': ('5\'-TA-3\'', '5\'-AT-3\''),
                       'gcgc_21mer': ('5\'-CG-3\'', '5\'-GC-3\''),
                       'ctct_21mer': ('5\'-TC-3\'', '5\'-CT-3\''),
                       'tgtg_21mer': ('5\'-GT-3\'', '5\'-TG-3\'')}

    def __init__(self, rootfolder, time_interval):
        super().__init__(rootfolder, time_interval)

    def get_xlabel(self, parameter):
        if parameter in ['shift', 'slide', 'rise']:
            return f'{parameter}(Å)'
        else:
            return f'{parameter}(°)'

    def histogram_bimodal(self, figsize, parameter, bins, xlines, xlim, ylim):
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=figsize, sharey=True, sharex=True)
        colors = ['blue', 'red']
        xlabel = self.get_xlabel(parameter)
        row_id = 0
        for host in ['atat_21mer', 'gcgc_21mer', 'ctct_21mer']:
            data1, data2 = self.get_bimodal_data(host, parameter)
            for col_id, data in enumerate([data1, data2]):
                ax = axes[row_id, col_id]
                label = f'{self.abbr_hosts[host]}:{self.d_bimodal_label[host][col_id]}'
                ax.hist(data, bins=bins, density=True, color=colors[col_id], alpha=0.4, 
                label=label)
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
            row_id += 1     
        return fig, ax

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
        temp_array_2 = np.ndarray.flatten(temp_array)
        zero_indice = np.isclose(temp_array_2, 0)
        return temp_array_2[~zero_indice]

    def get_bimodal_data(self, host, parameter):
        host_time_folder = path.join(self.rootfolder, host, self.time_interval)
        fname = path.join(host_time_folder, f'{parameter}.csv')
        df = pd.read_csv(fname, index_col='Frame-ID')
        n_frame = df.shape[0]

        temp_array = np.zeros((n_frame, self.n_bp))
        col_id = 0
        for bp_id in range(self.start_bp, self.end_bp):
            temp_array[:,col_id] = df[f'bp{bp_id}_bp{bp_id+1}']
            col_id += 1

        temp_array_1 = temp_array[::2]
        temp_array_2 = temp_array[1::2]

        temp_array_1_flatten = np.ndarray.flatten(temp_array_1)
        temp_array_2_flatten = np.ndarray.flatten(temp_array_2)
        zero_indice_1 = np.isclose(temp_array_1_flatten, 0)
        zero_indice_2 = np.isclose(temp_array_2_flatten, 0)

        return temp_array_1_flatten[~zero_indice_1], temp_array_2_flatten[~zero_indice_2]

class GrooveAgent(BasePairAgent):

    def __init__(self, rootfolder, time_interval):
        super().__init__(rootfolder, time_interval)
        self.start_bp = 4
        self.end_bp = 17
        self.n_bp = self.end_bp - self.start_bp + 1

    def get_xlabel(self, parameter):
        d_para = {'major_gw_pp': 'Major Groove Width', 'minor_gw_pp': 'Minor Groove Width'}
        return f'{d_para[parameter]}(Å)'

    def get_data(self, host, parameter):
        host_time_folder = path.join(self.rootfolder, host, self.time_interval)
        fname = path.join(host_time_folder, f'{parameter}.csv')
        df = pd.read_csv(fname, index_col='Frame-ID')
        n_frame = df.shape[0]
        temp_array = np.zeros((n_frame, self.n_bp))
        col_id = 0
        for bp_id in range(self.start_bp, self.end_bp+1):
            temp_array[:,col_id] = df[f'label{bp_id}']
            col_id += 1
        return np.ndarray.flatten(temp_array)

    def histogram_six(self, figsize, parameter, bins=100):
        n_rows = 2
        n_cols = 3
        fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize, sharex=True, sharey=True)
        d_axes = self.get_d_axes_by_host(axes, n_rows, n_cols)
        for host in self.hosts:
            ax = d_axes[host]
            data = self.get_data(host, parameter)
            mean = np.mean(data)
            std = np.std(data)
            ax.hist(data, color=self.d_colors[host], density=True, bins=bins, label=self.abbr_hosts[host])
            ax.axvline(mean, color='black', alpha=0.5)
            self.plot_assist_x(ax, parameter)
            self.plot_assist_y(ax, parameter)
            ax.set_title(self.get_title(host, mean, std), fontsize=10)
            if host in ['g_tract_21mer', 'gcgc_21mer', 'tgtg_21mer']:
                ax.set_xlabel(self.get_xlabel(parameter))
            if host in ['a_tract_21mer', 'g_tract_21mer']:
                ax.set_ylabel('Probability')
        return fig, d_axes

    def get_title(self, host, mean, std):
        return f'{self.abbr_hosts[host]}\n' + r'$\mu=' + f'{mean:.1f}$Å' + r'$~~\sigma=' + f'{std:.1f}$Å'

    def plot_assist_x(self, ax, parameter):
        d_xvalues = {'major_gw_pp': range(12, 28, 2), 'minor_gw_pp': range(8, 18, 2)}
        xvalues = d_xvalues[parameter]
        for xvalue in xvalues:
            ax.axvline(xvalue, color='grey', alpha=0.1)

    def plot_assist_y(self, ax, parameter):
        d_yvalues = {'major_gw_pp': np.arange(0.05, 0.26, 0.05), 'minor_gw_pp': np.arange(0.1, 0.6, 0.1)}
        yvalues = d_yvalues[parameter]
        for yvalue in yvalues:
            ax.axhline(yvalue, color='grey', alpha=0.1)

    def get_d_axes_by_host(self, axes, n_rows, n_cols):
        d_axes = dict()
        host_id = 0
        for row_id in range(n_rows):
            for col_id in range(n_cols):
                host = self.hosts[host_id]
                d_axes[host] = axes[row_id, col_id]
                host_id += 1
        return d_axes