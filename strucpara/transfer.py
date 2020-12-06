from os import path, system
from shutil import copyfile
from strucpara.miscell import check_dir_exist_and_make

class TransferAgent:
    gmx_local = '/usr/bin/gmx'
    type_na = 'bdna+bdna'
    time_series = ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us']
    collect_root = '/home/yizaochen/codes/dna_rna/collect_folder_to_multiscale'

    def __init__(self, allsys_folder, host, time_interval):
        self.allsys_folder = allsys_folder
        self.host = host
        self.time_interval = time_interval # '0_1us', '1_2us', '2_3us', '3_4us', '4_5us'

        self.local_folder = path.join(allsys_folder, host, self.type_na, 'input', 'allatoms')
        self.perfect_gro = path.join(self.local_folder, f'{self.type_na}.perfect.gro')
        self.perfect_pdb = path.join(self.local_folder, f'{self.type_na}.perfect.pdb')
        self.allxtc = path.join(self.local_folder, f'{self.type_na}.all.xtc')
        self.all_fitperfect_xtc = path.join(self.local_folder, f'{self.type_na}.all.fitperfect.xtc')

        self.collect_host_folder = path.join(self.collect_root, host)
        self.collect_folder = path.join(self.collect_host_folder, time_interval)
        self.col_perfect_pdb = path.join(self.collect_folder, f'{self.type_na}.perfect.pdb')
        self.col_all_fitperfect_xtc = path.join(self.collect_folder, f'{self.type_na}.all.fitperfect.xtc')
        self.compress_product = path.join(self.collect_folder, 'x3dna.required.tar.bz2')

        self.check_folders()

        # On server
        self.target_folder = '/home/yizaochen/x3dna/paper_2021'

    def check_folders(self):
        for folder in [self.collect_host_folder, self.collect_folder]:
            check_dir_exist_and_make(folder)

    def perfect_gro_to_pdb(self):
        cmd = f'{self.gmx_local} editconf -f {self.perfect_gro} -o {self.perfect_pdb}'
        system(cmd)
        print(cmd)

    def rmsd_fit_to_perfect(self):
        cmd = f'echo 0 0 | {self.gmx_local} trjconv -fit rot+trans -s {self.perfect_pdb} -f {self.allxtc} -o {self.all_fitperfect_xtc}'
        system(cmd)
        print(cmd)

    def copy_to_collect_folder(self):
        copyfile(self.perfect_pdb, self.col_perfect_pdb)
        print(f'cp {self.perfect_pdb} {self.col_perfect_pdb}')
        copyfile(self.all_fitperfect_xtc, self.col_all_fitperfect_xtc)
        print(f'cp {self.all_fitperfect_xtc} {self.col_all_fitperfect_xtc}')

    def compress_input(self):
        cmd = f'tar -jcv -f {self.compress_product} {self.collect_folder}'
        system(cmd)
        print(cmd)

    def scp_to_server(self, serverip):
        print('Please excute the following in the terminal:')
        cmd = f'scp {self.compress_product} yizaochen@{serverip}:{self.target_folder}'
        print(cmd)

    def decompress_in_server(self):
        target_file = 'x3dna.required.tar.bz2'
        print('Please excute the following in the terminal:')
        cmd = f'cd {self.target_folder}'
        print(cmd)
        cmd = f'tar -jxv -f {target_file} -C ./'
        print(cmd)
        cmd = f'rm {target_file}'
        print(cmd)
