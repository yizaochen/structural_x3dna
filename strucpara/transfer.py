from os import path, system
from shutil import copyfile, make_archive, move
from strucpara.miscell import check_dir_exist_and_make

class TransferAgent:
    gmx_local = '/usr/bin/gmx'
    type_na = 'bdna+bdna'
    time_series = ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us']
    collect_root = '/home/yizaochen/codes/dna_rna/collect_folder_to_multiscale'
    d_mdnum = {'0_1us': (1, 10), '1_2us': (11, 20), '2_3us': (21, 30), 
               '3_4us': (31, 40), '4_5us': (41, 50)}

    def __init__(self, allsys_folder, host, time_interval):
        self.allsys_folder = allsys_folder
        self.host = host
        self.time_interval = time_interval # '0_1us', '1_2us', '2_3us', '3_4us', '4_5us'

        self.local_folder = path.join(allsys_folder, host, self.type_na, 'input', 'allatoms')
        self.perfect_gro = path.join(self.local_folder, f'{self.type_na}.perfect.gro')
        self.perfect_pdb = path.join(self.local_folder, f'{self.type_na}.perfect.pdb')

        self.allxtc = path.join(self.local_folder, f'{self.type_na}.{self.time_interval}.xtc')
        self.all_fitperfect_xtc = path.join(self.local_folder, f'{self.type_na}.{self.time_interval}.fitperfect.xtc')

        self.collect_host_folder = path.join(self.collect_root, host)
        self.collect_folder = path.join(self.collect_host_folder, time_interval)
        self.col_perfect_pdb = path.join(self.collect_folder, f'{self.type_na}.perfect.pdb')
        self.col_all_fitperfect_xtc = path.join(self.collect_folder, f'{self.type_na}.all.fitperfect.xtc')

        self.zipfile = path.join(self.collect_host_folder, f'{self.host}.{self.time_interval}.zip')

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

    def concatenate_trajectory(self, simu_folder):
        start_mdnum = self.d_mdnum[self.time_interval][0]
        end_mdnum = self.d_mdnum[self.time_interval][1]

        na_folder = path.join(simu_folder, self.host, self.type_na)
        roughdir = path.join(na_folder, 'data','roughtrj','1000')
        alltrajfiles = ''
        for mdnum in range(start_mdnum, end_mdnum+1):
            filename = path.join(roughdir, '{0}.nopbc.fit.{1}.1000.xtc'.format(self.type_na, str(mdnum)))
            alltrajfiles = alltrajfiles + filename + ' '
        command = f'{self.gmx_local} trjcat -f {alltrajfiles} -o {self.allxtc} -dt 100'
        print(command)
        system(command)

    def remove_temp_xtc(self):
        command = f'rm {self.allxtc}'
        print(command)
        system(command)

        command = f'rm {self.all_fitperfect_xtc}'
        print(command)
        system(command)

    def rmsd_fit_to_perfect(self):
        cmd = f'echo 0 0 | {self.gmx_local} trjconv -fit rot+trans -s {self.perfect_pdb} -f {self.allxtc} -o {self.all_fitperfect_xtc}'
        system(cmd)
        print(cmd)

    def copy_to_collect_folder(self):
        copyfile(self.perfect_pdb, self.col_perfect_pdb)
        print(f'cp {self.perfect_pdb} {self.col_perfect_pdb}')
        copyfile(self.all_fitperfect_xtc, self.col_all_fitperfect_xtc)
        print(f'cp {self.all_fitperfect_xtc} {self.col_all_fitperfect_xtc}')

    def check_pdb_xtc_by_vmd(self):
        print(f'vmd -pdb {self.col_perfect_pdb} {self.col_all_fitperfect_xtc}')

    def compress_input(self):
        make_archive(self.collect_folder, 'zip', self.collect_host_folder, f'{self.time_interval}')
        old_name = path.join(self.collect_host_folder, f'{self.time_interval}.zip')
        move(old_name, self.zipfile)
        print(f'Archive {self.collect_folder} into {self.zipfile}')

    def scp_to_server(self, serverip):
        print('Please excute the following in the terminal:')
        cmd = f'scp {self.zipfile} yizaochen@{serverip}:{self.target_folder}'
        print(cmd)

    def decompress_in_server(self):
        target_file = f'{self.host}.{self.time_interval}.zip'
        print('Please excute the following in the terminal:')
        cmd = f'cd {self.target_folder}'
        print(cmd)
        cmd = f'unzip {target_file}'
        print(cmd)
        cmd = f'mv ./{self.time_interval} ./{self.host}'
        print(cmd)
        cmd = f'rm {target_file}'
        print(cmd)
