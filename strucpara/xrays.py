from os import path, system
from shutil import copyfile
from subprocess import Popen, PIPE
import numpy as np
from strucpara.miscell import check_dir_exist_and_make

class XrayAgent:

    d_n_bp = {'1d98': 12, '1d89': 12, '1bdn': 12, '1d65': 12, '1fzx': 12, '2dnd': 12, '287d': 12,
              '167d': 10, '1ilc': 12, '5mvk': 12, '1cgc': 10, '1dc0': 12, '1dn9': 12, '1nev': 10}

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
        self.lbph_dat = path.join(self.host_folder, 'local_base_pair_helical.dat')
        self.majorminor_dat = path.join(self.host_folder, 'major_minor.dat')

        self.d_lbp = None
        self.d_lbps = None
        self.d_lbph = None
        self.d_major_minor = None

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
        proc = Popen([find_pair, self.clean_pdb, self.inp_file], cwd=self.work_folder, stdout=PIPE)
        stddata = proc.communicate()
        print(f'{find_pair} {self.clean_pdb} {self.inp_file}')
        return stddata

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

    def get_local_base_pair_h_dat(self):
        n_after = self.d_n_bp[self.host]
        old_file = self.out_file
        new_file = self.lbph_dat
        cmd = f'grep -A {n_after} \'^Local base-pair helical parameters\' {old_file} > {new_file}'
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
        self.d_lbp = d_result
        print(f'{self.host}: Read Local Base Pair.')

    def read_local_base_pair_step_dat(self):
        n_bp = self.d_n_bp[self.host] - 1
        f_in = self.lbps_dat
        parameters = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = float(lines[bp_id, idx+2])
        self.d_lbps = d_result
        print(f'{self.host}: Read Local Base Pair Steps.')

    def read_local_base_pair_h_dat(self):
        n_bp = self.d_n_bp[self.host] - 1
        f_in = self.lbph_dat
        parameters = ['X-disp', 'Y-disp', 'h-Rise', 'Inclination', 'Tip', 'h-Twist']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = float(lines[bp_id, idx+2])
        self.d_lbph = d_result
        print(f'{self.host}: Read Local Base Pair Helical.')

    def read_major_minor_dat(self):
        n_bp = self.d_n_bp[self.host] - 1
        f_in = self.majorminor_dat
        parameters = ['Minor-P-P', 'Minor-Refined', 'Major-P-P',  'Major-Refined']
        d_result = {key: dict() for key in parameters}
        lines = np.genfromtxt(f_in, skip_header=2, skip_footer=0, dtype=object)
        for bp_id in range(n_bp):
            for idx, parameter in enumerate(parameters):
                d_result[parameter][bp_id+1] = lines[bp_id, idx+2]
        self.d_major_minor = d_result
        print(f'{self.host}: Read Major Minor Groove Width.')

    def vim_inp_and_outfile(self):
        print(f'vim {self.inp_file}')
        print(f'vim {self.out_file}')

    def vim_three_dats(self):
        for fname in [self.lbp_dat, self.lbps_dat, self.majorminor_dat]:
            print(f'vim {fname}')

class NMRBigAgnet:
    d_n_bp = {'1nev': 12}
    d_n_model = {'1nev': 10}
    d_n_atom = {'1nev': 632}

    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host
        self.n_model = self.d_n_model[self.host]

        self.host_folder = path.join(self.rootfolder, self.host)
        self.check_folders()

        self.d_model_agent = None

    def check_folders(self):
        for folder in [self.host_folder]:
            check_dir_exist_and_make(folder)

    def make_model_folder(self):
        for model_id in range(1, self.n_model+1):
            model_folder = path.join(self.host_folder, f'{model_id}')
            check_dir_exist_and_make(model_folder)

    def extract_pdb_to_model_folder(self):
        n_atom = self.d_n_atom[self.host]
        ref_pdb = path.join(self.host_folder, f'{self.host}_clean.pdb')
        for model_id in range(1, self.n_model+1):
            new_pdb = path.join(self.host_folder, f'{model_id}', f'{self.host}_clean.pdb')
            model_txt = f'MODEL{model_id:>9}'
            cmd = f'grep -A {n_atom} \'{model_txt}\' {ref_pdb} > {new_pdb}'
            system(cmd)
            print(cmd)

    def get_all_model_agents(self):
        d_x_agent = dict()
        for model_id in range(1, self.n_model+1):
            d_x_agent[model_id] = NMRAgent(self.rootfolder, self.host, model_id)
        self.d_model_agent = d_x_agent

    def find_all_pair(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].find_pair()

    def analyze_all(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].analyze()

    def make_all_dat(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].get_local_base_pair_dat()
            self.d_model_agent[model_id].get_local_base_pair_step_dat()
            self.d_model_agent[model_id].get_local_base_pair_h_dat()
            self.d_model_agent[model_id].get_major_minor_dat()

    def pre_processing(self):
        self.get_all_model_agents()
        self.read_all_local_base_pair_dat()
        self.read_all_local_base_pair_step_dat()
        self.read_all_local_base_pair_h_dat()
        self.read_all_major_minor_dat()

    def read_all_local_base_pair_dat(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].read_local_base_pair_dat()

    def read_all_local_base_pair_step_dat(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].read_local_base_pair_step_dat()

    def read_all_local_base_pair_h_dat(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].read_local_base_pair_h_dat()

    def read_all_major_minor_dat(self):
        for model_id in range(1, self.n_model+1):
            self.d_model_agent[model_id].read_major_minor_dat()

    def get_local_base_pair_list(self, parameter, sele_bps):
        data = list()
        for model_id in range(1, self.n_model+1):
            agent = self.d_model_agent[model_id]
            data += [agent.d_lbp[parameter][bp_id] for bp_id in sele_bps]
        return data

    def get_local_base_pair_step_list(self, parameter, sele_bps):
        data = list()
        for model_id in range(1, self.n_model+1):
            agent = self.d_model_agent[model_id]
            data += [agent.d_lbps[parameter][bp_id] for bp_id in sele_bps]
        return data

    def get_local_base_pair_h_list(self, parameter, sele_bps):
        data = list()
        for model_id in range(1, self.n_model+1):
            agent = self.d_model_agent[model_id]
            data += [agent.d_lbph[parameter][bp_id] for bp_id in sele_bps]
        return data

    def get_major_minor_list(self, parameter, sele_bps):
        data = list()
        for model_id in range(1, self.n_model+1):
            agent = self.d_model_agent[model_id]
            data += [float(agent.d_major_minor[parameter][bp_id]) for bp_id in sele_bps]
        return data


class NMRAgent(XrayAgent):

    def __init__(self, rootfolder, host, model_id):
        self.rootfolder = rootfolder
        self.host = host
        self.model_id = f'{model_id}'

        self.host_folder = path.join(self.rootfolder, self.host, self.model_id) #
        self.work_folder = path.join(self.host_folder, 'work')
        self.original_folder = path.join(self.rootfolder, 'original_pdbs')
        self.check_folders()

        self.original_pdb = path.join(self.original_folder, f'{host}.pdb')
        self.clean_pdb = path.join(self.host_folder, f'{host}_clean.pdb')

        self.inp_file = path.join(self.host_folder, f'{self.host}.inp')
        self.out_file = path.join(self.work_folder, f'{self.host}_clean.out')

        self.lbp_dat = path.join(self.host_folder, 'local_base_pair.dat')
        self.lbps_dat = path.join(self.host_folder, 'local_base_pair_step.dat')
        self.lbph_dat = path.join(self.host_folder, 'local_base_pair_helical.dat')
        self.majorminor_dat = path.join(self.host_folder, 'major_minor.dat')

        self.d_lbp = None
        self.d_lbps = None
        self.d_lbph = None
        self.d_major_minor = None


class XrayCollectAgent:

    hosts = ['1d98', '1d89', '1bdn', '1fzx', '1nev', '1d65', '167d', '1ilc', '2dnd', '287d', '1dn9', '1dc0', '1cgc', '5mvk']
    d_sele_bp_lbp = {'1d98': [5,6,7,8], '1d89': [6,7,8], '1bdn': [5,6,7], '1d65': [5,6,7,8], '1fzx': [5,6,7,8], 
                     '2dnd': [5,6,7,8], '287d': [6,7], '167d': [4,5,6,7], '1ilc': [6,7], '5mvk': [5,6,7,8], 
                     '1cgc': [4,5,6,7], '1dc0': [5,6,7,8], '1dn9': [5,6,7,8], '1nev': [5,6]}

    d_sele_bp_lbps = {'1d98': [5,6,7], '1d89': [6,7,8], '1bdn': [5,6,7], '1d65': [5,6,7], '1fzx': [5,6,7], 
                      '2dnd': [5,6,7], '287d': [6], '167d': [4,5,6], '1ilc': [5,6,7], '5mvk': [4,5,6,7,8], 
                      '1cgc': [4,5,6], '1dc0': [4,5,6,7], '1dn9': [4,5,6,7,8], '1nev': [4,5,6]}

    d_groups = {'Group 1': ['1d98', '1d89', '1bdn'], 'Group 2': ['1fzx', '1nev'], 
                'Group 3': ['1d65', '167d', '1ilc', '2dnd'], 'Group 4': ['287d', '1dn9'],
                'Group 5': ['1dc0', '1cgc', '5mvk']}

    d_colors_group = {'Group 1': 'b', 'Group 2': 'deepskyblue', 'Group 3': 'cyan',  
                      'Group 4': 'blueviolet', 'Group 5': 'red'}

    group_list = ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']

    def __init__(self, rootfolder):
        self.rootfolder = rootfolder
        self.d_x_agent = self.get_all_xagents()

    def get_all_xagents(self):
        d_x_agent = dict()
        for host in self.hosts:
            if host in ['1nev']:
                d_x_agent[host] = NMRBigAgnet(self.rootfolder, host)
                d_x_agent[host].get_all_model_agents()
            else:
                d_x_agent[host] = XrayAgent(self.rootfolder, host)
        return d_x_agent

    def read_all_local_base_pair_dat(self):
        for host in self.hosts:
            if host in ['1nev']:
                self.d_x_agent[host].read_all_local_base_pair_dat()
            else:
                self.d_x_agent[host].read_local_base_pair_dat()

    def read_all_local_base_pair_step_dat(self):
        for host in self.hosts:
            if host in ['1nev']:
                self.d_x_agent[host].read_all_local_base_pair_step_dat()
            else:
                self.d_x_agent[host].read_local_base_pair_step_dat()

    def read_all_local_base_pair_h_dat(self):
        for host in self.hosts:
            if host in ['1nev']:
                self.d_x_agent[host].read_all_local_base_pair_h_dat()
            else:
                self.d_x_agent[host].read_local_base_pair_h_dat()

    def read_all_major_minor_dat(self):
        for host in self.hosts:
            if host in ['1nev']:
                self.d_x_agent[host].read_all_major_minor_dat()
            else:
                self.d_x_agent[host].read_major_minor_dat()

    def get_data_lbp(self, host, parameter):
        agent = self.d_x_agent[host]
        sele_bps = self.d_sele_bp_lbp[host]
        if host in ['1nev']:
            return agent.get_local_base_pair_list(parameter, sele_bps)
        else:
            return [agent.d_lbp[parameter][bp_id] for bp_id in sele_bps]

    def get_data_lbps(self, host, parameter):
        agent = self.d_x_agent[host]
        sele_bps = self.d_sele_bp_lbps[host]
        if host in ['1nev']:
            return agent.get_local_base_pair_step_list(parameter, sele_bps)
        else:
            return [agent.d_lbps[parameter][bp_id] for bp_id in sele_bps]

    def get_data_lbph(self, host, parameter):
        agent = self.d_x_agent[host]
        sele_bps = self.d_sele_bp_lbps[host]
        if host in ['1nev']:
            return agent.get_local_base_pair_h_list(parameter, sele_bps)
        else:
            return [agent.d_lbph[parameter][bp_id] for bp_id in sele_bps]

    def get_data_major_minor(self, host, parameter):
        agent = self.d_x_agent[host]
        sele_bps = self.d_sele_bp_lbps[host]
        if host in ['1nev']:
            return agent.get_major_minor_list(parameter, sele_bps)
        else:
            return [float(agent.d_major_minor[parameter][bp_id]) for bp_id in sele_bps]

    def get_box_positions(self):
        d_position = dict()
        idx = 0
        for group_id in self.group_list:
            n_host = len(self.d_groups[group_id])
            d_position[group_id] = list()
            for i in range(n_host):
                d_position[group_id].append(idx+i)
            idx += n_host
        return d_position
                
    def get_local_base_pair_data_for_boxplot(self, parameter):
        data = dict()
        for group_id in self.group_list:
            hosts = self.d_groups[group_id]
            data[group_id] = [self.get_data_lbp(host, parameter) for host in hosts]
        return data

    def get_local_base_pair_step_data_for_boxplot(self, parameter):
        data = dict()
        for group_id in self.group_list:
            hosts = self.d_groups[group_id]
            data[group_id] = [self.get_data_lbps(host, parameter) for host in hosts]
        return data

    def get_local_base_pair_h_data_for_boxplot(self, parameter):
        data = dict()
        for group_id in self.group_list:
            hosts = self.d_groups[group_id]
            data[group_id] = [self.get_data_lbph(host, parameter) for host in hosts]
        return data

    def get_major_minor_data_for_boxplot(self, parameter):
        data = dict()
        for group_id in self.group_list:
            hosts = self.d_groups[group_id]
            data[group_id] = [self.get_data_major_minor(host, parameter) for host in hosts]
        return data
