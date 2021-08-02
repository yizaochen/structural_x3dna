from os import path, system
from io import StringIO
import pandas as pd
from strucpara.miscell import check_dir_exist_and_make

ensemble_exec = '/home/yizaochen/opt/x3dna-v2.3/bin/x3dna_ensemble'

class BasePairAgent:

    def __init__(self, rootfolder, host, time_interval):
        self.rootfolder = rootfolder
        self.type_na = 'bdna+bdna'
        self.n_bp = 21
        self.host = host
        self.time_interval = time_interval
        self.host_folder = path.join(rootfolder, host)
        self.host_time_folder = path.join(self.host_folder, time_interval)

        self.check_folder()

        self.ensemble_out = path.join(self.host_time_folder, f'{self.type_na}.ensemble.out')

        self.server_root = '/home/yizaochen/x3dna/paper_2021'
        self.host_time_server = path.join(self.server_root, host, time_interval)
        self.ensemble_out_server = path.join(self.host_time_server, f'{self.type_na}.ensemble.out')
        
        self.parameters = ['shear', 'buckle', 'stretch', 'propeller', 'stagger', 'opening']


    def check_folder(self):
        for folder in [self.host_folder, self.host_time_folder]:
            check_dir_exist_and_make(folder)

    def download_ensesmble_out(self, serverip):
        print('Please excute the following in the terminal:')
        cmd = f'scp yizaochen@{serverip}:{self.ensemble_out_server} {self.ensemble_out}'
        print(cmd)
        
    def extract_parameters(self):
        for parameter in self.parameters:
            output_dat = path.join(self.host_time_folder, f'{parameter}.dat')
            cmd = f'{ensemble_exec} extract -f {self.ensemble_out} -p {parameter} -o {output_dat}'
            system(cmd)
            print(cmd)

    def convert_dat_to_csv(self):
        for parameter in self.parameters:
            dat_in = path.join(self.host_time_folder, f'{parameter}.dat')
            f = open(dat_in, 'r')
            lines = f.readlines()
            f.close()
            first_line = self.get_first_line()
            lines = first_line + lines
            buffer = StringIO(''.join(lines))
            df = pd.read_csv(buffer, sep="\t")
            csv_out = path.join(self.host_time_folder, f'{parameter}.csv')
            df.to_csv(csv_out, index=False)
            print(f'Dataframe to csv: {csv_out}')

    def get_first_line(self):
        result = ['Frame-ID']
        for bp_id in range(1, self.n_bp+1):
            if bp_id == self.n_bp:
                result.append(f'bp{bp_id}\n')
            else:
                result.append(f'bp{bp_id}')
        return ['\t'.join(result)]

    def clean_dat_files(self):
        cmd = f'rm {self.host_time_folder}/*.dat'
        print("Please execute the following on Terminal:")
        print(cmd)

class BaseStepAgent(BasePairAgent):
    def __init__(self, rootfolder, host, time_interval):
        super().__init__(rootfolder, host, time_interval)
        self.parameters = ['shift', 'tilt', 'slide', 'roll', 'rise', 'twist']

    def get_first_line(self):
        result = ['Frame-ID']
        for bp_id in range(1, self.n_bp):
            if bp_id == (self.n_bp-1):
                result.append(f'bp{bp_id}_bp{bp_id+1}\n')
            else:
                result.append(f'bp{bp_id}_bp{bp_id+1}')
        return ['\t'.join(result)]

class BaseHelicalAgent(BasePairAgent):
    def __init__(self, rootfolder, host, time_interval):
        super().__init__(rootfolder, host, time_interval)
        self.parameters = ['x_displacement', 'y_displacement', 'inclination', 'tip', 'h_rise', 'h_twist']

    def get_first_line(self):
        result = ['Frame-ID']
        for bp_id in range(1, self.n_bp):
            if bp_id == (self.n_bp-1):
                result.append(f'bp{bp_id}_bp{bp_id+1}\n')
            else:
                result.append(f'bp{bp_id}_bp{bp_id+1}')
        return ['\t'.join(result)]


class GrooveAgent(BasePairAgent):
    def __init__(self, rootfolder, host, time_interval):
        super().__init__(rootfolder, host, time_interval)
        self.parameters = ['major_gw_pp', 'major_gw_refined', 'minor_gw_pp', 'minor_gw_refined']

    def get_first_line(self):
        result = ['Frame-ID']
        for bp_id in range(1, self.n_bp):
            if bp_id == (self.n_bp-1):
                result.append(f'label{bp_id}\n')
            else:
                result.append(f'label{bp_id}')
        return ['\t'.join(result)]

class PhaseChiDeltaAgent(BasePairAgent):
    def __init__(self, rootfolder, host, time_interval):
        super().__init__(rootfolder, host, time_interval)
        self.parameters = ['phase1', 'phase2', 'chi1', 'chi2', 'delta1', 'delta2']

    def get_first_line(self):
        result = ['Frame-ID']
        for bp_id in range(1, self.n_bp+1):
            if bp_id == self.n_bp:
                result.append(f'Base{bp_id}\n')
            else:
                result.append(f'Base{bp_id}')
        return ['\t'.join(result)]