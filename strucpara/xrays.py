from os import path, system
from shutil import copyfile
from subprocess import Popen, PIPE
from strucpara.miscell import check_dir_exist_and_make
import numpy as np

class XrayAgent:

    d_n_bp = {'1d98': 12, '1d89': 12, '1bdn': 12, '1d65': 12, '1fzx': 12, '2dnd': 12, '287d': 12,
              '167d': 10, '1ilc': 12, '5mvk': 12, '1cgc': 10, '1dc0': 12, '1dn9': 12}

    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host

        self.host_folder = path.join(self.rootfolder, self.host)
        self.work_folder = path.join(self.host_folder, 'work')
        self.original_folder = path.join(self.rootfolder, 'original_pdbs')
        self.check_folders()

        self.original_pdb = path.join(self.original_folder, f'{host}.pdb')
        self.clean_pdb = path.join(self.host_folder, f'{host}_clean.pdb')

        self.inp_file = path.join(self.host_folder, f'{self.host}.inp')
        self.out_file = path.join(self.work_folder, f'{self.host}_clean.out')

        self.lbp_dat = path.join(self.host_folder, 'local_base_pair.dat')
        self.lbps_dat = path.join(self.host_folder, 'local_base_pair_step.dat')
        self.majorminor_dat = path.join(self.host_folder, 'major_minor.dat')

    def check_folders(self):
        for folder in [self.host_folder, self.work_folder]:
            check_dir_exist_and_make(folder)
    
    def make_clean_pdb(self):
        copyfile(self.original_pdb, self.clean_pdb)
        print(f'cp {self.original_pdb} {self.clean_pdb}')
        print(f'vim {self.clean_pdb}')

    def check_clean_pdb_by_vmd(self):
        print(f'vmd -pdb {self.clean_pdb}')

    def find_pair(self):
        find_pair = '/home/yizaochen/opt/x3dna-v2.3/bin/find_pair'
        cmd = f'{find_pair} {self.clean_pdb} {self.inp_file}'
        system(cmd)
        print(cmd)

    def analyze(self):
        analyze = '/home/yizaochen/opt/x3dna-v2.3/bin/analyze'
        proc = Popen([analyze, self.inp_file], cwd=self.work_folder, stdout=PIPE)
        print(f'{analyze} {self.inp_file}')
        stddata = proc.communicate()
        return stddata

    def get_local_base_pair_dat(self):
        n_after = self.d_n_bp[self.host] + 1
        old_file = self.out_file
        new_file = self.lbp_dat
        cmd = f'grep -A {n_after} \'^Local base-pair parameters\' {old_file} > {new_file}'
        system(cmd)
        print(cmd)

    def get_local_base_pair_step_dat(self):
        n_after = self.d_n_bp[self.host]
        old_file = self.out_file
        new_file = self.lbps_dat
        cmd = f'grep -A {n_after} \'^Local base-pair step parameters\' {old_file} > {new_file}'
        system(cmd)
        print(cmd)

    def get_major_minor_dat(self):
        n_after = self.d_n_bp[self.host]
        old_file = self.out_file
        new_file = self.majorminor_dat
        cmd = f'grep -A {n_after} \'^                  Minor Groove        Major Groove\' {old_file} > {new_file}'
        system(cmd)
        print(cmd)

    def read_local_base_pair_dat(self):
        n_bp = self.d_n_bp[self.host]
        f_in = self.lbp_dat
        parameters = ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = float(lines[bp_id, idx+2])
        return d_result

    def read_local_base_pair_step_dat(self):
        n_bp = self.d_n_bp[self.host] - 1
        f_in = self.lbps_dat
        parameters = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = float(lines[bp_id, idx+2])
        return d_result

    def read_major_minor_dat(self):
        n_bp = self.d_n_bp[self.host] - 1
        f_in = self.lbps_dat
        parameters = ['Minor-P-P', 'Minor-Refined', 'Major-P-P',  'Major-Refined']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = lines[bp_id, idx+2]
        return d_result

    def vim_inp_and_outfile(self):
        print(f'vim {self.inp_file}')
        print(f'vim {self.out_file}')

    def vim_three_dats(self):
        for fname in [self.lbp_dat, self.lbps_dat, self.majorminor_dat]:
            print(f'vim {fname}')