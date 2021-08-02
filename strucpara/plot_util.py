import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

class FourHostSplitStrand:

    hosts = ['a_tract_21mer', 'tat_21mer', 'g_tract_21mer', 'gcgc_21mer']
    strands = ['STRAND1', 'STRAND2']
    d_outer_coor = {'a_tract_21mer': (0,0), 'tat_21mer': (0,1), 
                    'g_tract_21mer': (1,0), 'gcgc_21mer': (1,1)}
    d_inner_coor = {'STRAND1': 0, 'STRAND2': 1}    

    def __init__(self, figsize, outer_wspace, outer_hspace, inner_wspace, inner_hspace):
        self.figsize = figsize
        self.outer_wspace = outer_wspace
        self.outer_hspace = outer_hspace
        self.inner_wspace = inner_wspace
        self.inner_hspace = inner_hspace

        self.fig = plt.figure(figsize=self.figsize)
        self.outer_grid = GridSpec(2, 2, wspace=self.outer_wspace, hspace=self.outer_hspace, figure=self.fig)
        self.d_inner_grids = self.ini_inner_grids()
        self.d_axes = self.ini_d_axes()

    def ini_inner_grids(self):
        d_inner_grids = dict()
        for host in self.hosts:
            out_coor = self.d_outer_coor[host]
            d_inner_grids[host] = GridSpecFromSubplotSpec(2, 1, subplot_spec=self.outer_grid[out_coor], 
            wspace=self.inner_wspace, hspace=self.inner_hspace)
        return d_inner_grids

    def ini_d_axes(self):
        d_axes = {host: {'STRAND1': None, 'STRAND2': None} for host in self.hosts}
        for host in self.hosts:
            for strand in self.strands:
                inner_coor = self.d_inner_coor[strand]
                d_axes[host][strand] = self.fig.add_subplot(self.d_inner_grids[host][inner_coor])
        return d_axes

    def get_fig(self):
        return self.fig

    def get_d_axes(self):
        return self.d_axes
