from os import path, system


class InitialLocalAgent:
    gmx = '/usr/bin/gmx'
    rootfolder = '/home/yizaochen/codes/dna_rna/collect_folder_to_multiscale'
    type_na = 'bdna+bdna'

    def __init__(self, host, time_interval):
        self.host = host
        self.time_interval = time_interval
        self.host_folder = path.join(self.rootfolder, host, time_interval)
        self.ref_pdb = path.join(self.host_folder, f'{self.type_na}.perfect.pdb')
        self.in_xtc = path.join(self.host_folder, f'{self.type_na}.all.fitperfect.xtc')
        self.out_pdb = path.join(self.host_folder, f'{self.type_na}.all.fitperfect.pdb')

        self.xtc2pdb_qsub = path.join(self.host_folder, 'xtc2pdb.qsub')

    def xtc2pdb(self):
        cmd = f'echo 0 | {self.gmx} trjconv -s {self.ref_pdb} -f {self.in_xtc} -o {self.out_pdb}'
        print(cmd)
        system(cmd)

    def remove_out_pdb(self):
        cmd = f'rm {self.out_pdb}'
        print(cmd)
        system(cmd)

class InitialAgent:
    gmx_server = '/home/tclick/usr/gromacs/bin/gmx'
    rootfolder = '/home/yizaochen/x3dna/paper_2021'
    type_na = 'bdna+bdna'

    def __init__(self, host, time_interval):
        self.host = host
        self.time_interval = time_interval
        self.host_folder = path.join(self.rootfolder, host, time_interval)
        self.ref_pdb = path.join(self.host_folder, f'{self.type_na}.perfect.pdb')
        self.in_xtc = path.join(self.host_folder, f'{self.type_na}.all.fitperfect.xtc')
        self.out_pdb = path.join(self.host_folder, f'{self.type_na}.all.fitperfect.pdb')

        self.xtc2pdb_qsub = path.join(self.host_folder, 'xtc2pdb.qsub')

    def write_xtc2pdb_qsub(self):
        f = open(self.xtc2pdb_qsub, 'w')
        f.write('#!/bin/bash -l\n')
        f.write(f'#PBS -N {self.host}_xtc2pdb\n')
        f.write('#PBS -l walltime=48:00:00\n')
        f.write('#PBS -q batch\n')
        f.write('#PBS -l nodes=1:ppn=16\n')
        f.write('#PBS -j oe\n')
        f.write(f'#PBS -o /home/yizaochen/log/{self.host}_xtc2pdb.log\n')
        f.write('#PBS -r n\n\n')
        cmd = f'echo 0 | {self.gmx_server} trjconv -s {self.ref_pdb} -f {self.in_xtc} -o {self.out_pdb}'
        f.write(f'{cmd}\n')
        f.close()

    def qsub_xtc2pdb(self):
        cmd = f'qsub {self.xtc2pdb_qsub}'
        system(cmd)
        print(cmd)
        print('Log file is:')
        print(f'/home/yizaochen/log/{self.host}_xtc2pdb.log')