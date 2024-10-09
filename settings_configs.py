'''
To get the list of server paths, just use the 'server_paths' dict.
'''

import os.path, platform

import configparser
parser = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
parser.read('server_paths.ini')
server_paths = {i: {j: parser[i][j] for j in list(parser[i].keys())} for i in parser.sections()}

def server_paths_as_bash():
    '''Prints the server_paths to the console in a manner that can be copied in to a bash script.'''
    sects = list(server_paths.keys())
    for s in sects:
        print(f'### {s}')
        for k,v in server_paths[s].items():
            print(f'{k}={v}')


def full_split_path(mypa):
    '''
    utility to convert a path to a list of folder names so we can deal with the
    WSL/windows platform issues.
    '''
    sppa = []
    pcs = os.path.split(mypa); i=0;
    while pcs[1]!='' and i< 10000:
        sppa.insert(0, pcs[1])
        pcs = os.path.split(pcs[0])
        i += 1
    sppa.insert(0,pcs[0])
    return sppa

def subfolder_list_to_path(sflist):
    return os.path.join(sflist)

# Folders, broadly:
local_wgs_gut_work=['mnt', 'c', 'Users', 'miken', 'Rice', 'Research', 'gut_mb_autoencoder']

shared_paths_doc = '''
bioproject_fastq_locs:  shows where all the RAW fastq.gz files are stored 
                        for each bioproject.
wgs_gut_work:           main working folder for the project
'''


wgs_gut_work = {'laptop': os.path.join(*local_wgs_gut_work),
                'server': '/home/mnute/work/rice/wgs_sra_gut_work'}
metadata = {'server': os.path.join(wgs_gut_work['server'], 'metadata_tracking'),
            'laptop': os.path.join('mnt', 'c', 'Users', 'miken', 'Rice', 'Research', 'gut_mb_autoencoder', 'metadata_tracking')}
dir_paths = [wgs_gut_work, metadata]


# bioproject_fastq_locs= {'filename': 'target_bioproject_path_list.tsv',
#                         'folder': metadata}
bioproject_fastq_locs = {'server': os.path.join(metadata['server'], 'target_bioproject_path_list.tsv'),
                         'laptop': os.path.join(metadata['laptop'], 'target_bioproject_path_list.tsv')}
sra_runslist_default_dumpfile={'filename': 'Master_Project_Dataset_Runlist.txt',
                               'server': os.path.join(metadata['server'], 'Master_Project_Dataset_Runlist.txt'),
                               'laptop': os.path.join(metadata['laptop'], 'Master_Project_Dataset_Runlist.txt')}
fastq_gz_inventory = {'filename': 'download_inventory.tsv',
                      'folder': metadata}
postproc_file_inventory = {'filename': 'postproc_file_inventory.tsv',
                           'folder': metadata}

file_paths = [bioproject_fastq_locs, sra_runslist_default_dumpfile, fastq_gz_inventory, postproc_file_inventory]

megahit_output_folders = {
    'main': '/hdd1/mnute/wgs_sra_gut/megahit_output',
    'transfer': '/hdd1/mnute/wgs_sra_gut/megahit_transfer'
}
megahit_output_inventory = os.path.join(metadata['server'], 'megahit_output_inventory.tsv')


def resolve_dir(var):
    which_machine = 'server' if platform.node()=='nuteserver' else 'laptop'
    if which_machine in var.keys():
        choice = var[which_machine]
        if isinstance(choice,list):
            choice = os.path.join(choice)
    return choice

def resolve_file(var):
    return os.path.join(resolve_dir(var['folder']), var['filename'])



# for i in shared_paths:
#     shared_paths[i]['loc'] = os.path.join(shared_paths[i]['locsubs'])
