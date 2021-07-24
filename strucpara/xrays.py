from os import path
from shutil import copyfile
from strucpara.miscell import check_dir_exist_and_make

class XrayAgent:

    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host

        self.host_folder = path.join(self.rootfolder, self.host)
        self.original_folder = path.join(self.rootfolder, 'original_pdbs')
        self.check_folders()

        self.original_pdb = path.join(self.original_folder, f'{host}.pdb')
        self.clean_pdb = path.join(self.original_folder, f'{host}_clean.pdb')

    def check_folders(self):
        for folder in [self.host_folder]:
            check_dir_exist_and_make(folder)
    
    def make_clean_pdb(self):
        copyfile(self.original_pdb, self.clean_pdb)
        print(f'cp {self.original_pdb} {self.clean_pdb}')
        print(f'vim {self.clean_pdb}')

    def check_clean_pdb_by_vmd(self):
        print(f'vmd -pdb {self.clean_pdb}')