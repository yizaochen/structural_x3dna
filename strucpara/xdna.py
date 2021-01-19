from os import path, system
from miscell import check_dir_exist_and_make
class X3DNAAgent:
    x3dna_root = '/home/yizaochen/x3dna-v2.3/bin'
    ensemble_exec = path.join(x3dna_root, 'x3dna_ensemble')
    python_exec = '/home/yizaochen/miniconda3/envs/x3dna/bin/python'
    rootfolder = '/home/yizaochen/x3dna/paper_2021'
    type_na = 'bdna+bdna'
    yz_scr = '/scratch/yizaochen'

    def __init__(self, host, time_interval):
        self.host = host
        self.time_interval = time_interval
        self.host_folder = path.join(self.rootfolder, host, time_interval)
        self.ref_pdb = path.join(self.host_folder, f'{self.type_na}.perfect.pdb')
        self.out_pdb = path.join(self.host_folder, f'{self.type_na}.all.fitperfect.pdb')
        self.findpair_inp = path.join(self.host_folder, f'{self.type_na}.inp')
        self.ensemble_qsub = path.join(self.host_folder, 'ensemble_analyze.qsub')
        self.ensemble_output = path.join(self.host_folder, f'{self.type_na}.ensemble.out')

        self.check_py = path.join(self.host_folder, 'check_scatch.py')

        self.x3dna_scr = path.join(self.yz_scr, 'x3dna')
        self.host_scr = path.join(self.x3dna_scr, host)
        self.sys_scr = path.join(self.host_scr, time_interval)
        self.ensemble_scr = path.join(self.sys_scr, f'{self.type_na}.ensemble.out')

    def find_pair(self):
        cmd = f'find_pair {self.ref_pdb} {self.findpair_inp}'
        system(cmd)
        print(cmd)
        print(f'Generate input file: {self.findpair_inp}')

    def write_check_scratch_py(self):
        f = open(self.check_py, 'w')
        f.write('from strucpara.xdna import X3DNAAgent\n')
        f.write(f'agent=X3DNAAgent(\'{self.host}\', \'{self.time_interval}\')\n')
        f.write('agent.check_scratch()\n')
        f.close()

    def check_scratch(self):
        for folder in [self.yz_scr, self.x3dna_scr, self.host_scr, self.sys_scr]:
            check_dir_exist_and_make(folder)

    def write_ensemble_analyze_qsub(self):
        f = open(self.ensemble_qsub, 'w')
        f.write('#!/bin/bash -l\n')
        f.write(f'#PBS -N {self.host}_ensemble_analyze\n')
        f.write('#PBS -l walltime=48:00:00\n')
        f.write('#PBS -q batch\n')
        f.write('#PBS -l nodes=1:ppn=16\n')
        f.write('#PBS -j oe\n')
        f.write(f'#PBS -o /home/yizaochen/log/{self.host}_ensemble_analyze.log\n')
        f.write('#PBS -r n\n\n')

        f.write('conda activate x3dna\n')
        f.write(f'{self.python_exec} {self.check_py}\n\n')

        f.write(f'cd {self.sys_scr}\n')

        f.write(f'{self.ensemble_exec} analyze -b {self.findpair_inp} -e {self.out_pdb} -o {self.ensemble_scr}\n')
        f.write(f'cp {self.ensemble_scr} {self.ensemble_output}\n')
        f.close()

    def qsub_ensemble_analyse(self):
        cmd = f'qsub {self.ensemble_qsub}'
        system(cmd)
        print(cmd)
        print('Log file is:')
        print(f'/home/yizaochen/log/{self.host}_ensemble_analyze.log')