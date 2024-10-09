import os, sys, datetime, argparse, logging, glob, re
import shutil

import numpy as np
import phylogeny_utilities.utilities as phy
from settings_configs import *
import main

parser = argparse.ArgumentParser()
rich_format = "[%(filename)s (%(lineno)d) %(asctime)s] %(levelname)s: %(message)s"
logging.basicConfig(format=rich_format, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)

def parse_command_args():
    global parser

    subparsers = parser.add_subparsers(help='specifies the action to take.')

    megahit_stats = subparsers.add_parser('megahit_stats',
                                         help='Goes through a megahit results folder and gets the summary statistics '
                                              'from the `final.contigs.fa` output file.')
    megahit_stats.add_argument('-f', '--megahit_folder', dest='megahit_folder', type=str, required=True,
                              help='Path to the folder containing the megahit output. If this is exactly equal to the '
                                   'string `headers` then the field names are printed instead (e.g. for creating a '
                                   'tabular summary with field names in the first row).')
    megahit_stats.add_argument('-d', '--delimiter', dest='delimiter', type=str, default=',',
                               help='Delimiter to use in separating output fields.')
    megahit_stats.add_argument('-m', '--min_length', dest='min_length', type=int, default=500,
                               help='Minimum contig length to consider when reporting assembly statistics.')
    megahit_stats.set_defaults(func=assembly_summary_stats_to_console)

    megahit_get_output_inventory = subparsers.add_parser('megahit_output_inventory',
                                            help='Goes through the megahit output folders and takes inventory '
                                                 'of all completed assemblies, plus gets statistics, and outputs results '
                                                 'to the megahit_output_inventory file in the metadata folder.')
    megahit_get_output_inventory.set_defaults(func=walk_megahit_output_folders)

    args = parser.parse_args()
    return args

def assembly_summary_stats_to_console(myargs):
    '''Command line function to get the assembly stats and print them to the console in delimited form. For use
    with the subcommand `megahit_stats`.'''
    folder_path = myargs.megahit_folder
    delimiter = myargs.delimiter
    min_length = myargs.min_length

    if folder_path=='headers':
        out_tup = get_assembly_summary_stats('', field_names_only=True)
    else:
        out_tup = get_assembly_summary_stats(folder_path, min_length)

    print(delimiter.join(map(str, out_tup)), flush=True)



def get_assembly_summary_stats(megahit_dir, min_length=500, contigs_filename='final.contigs.fa',
                               field_names_only=False):
    '''Goes into the megahit folder and opens the contigs file and gets several summary statistics about the
    assembly, which are returned as a tuple. Returns the following stats:
        1) # Contigs
        2) Total Aligned Length
        3) Avg Contig Length
        4) Longest Contig Length
        5) N25 (length of contig s/t 25% of total aligned length are in contigs longer than it)
        6) N50
        7) N75
    '''
    field_names = [
        'folder',
        'min_length',
        'num_contigs',
        'total_length',
        'avg_length',
        'longest_contig',
        'N25','N50','N75',
        'filename'
    ]
    if field_names_only:
        return tuple(field_names)
    fa = phy.read_from_fasta(os.path.join(megahit_dir, contigs_filename))
    logger.debug('Read fasta, got %d contigs' % len(fa))
    fa_lens = np.sort(np.array([len(i) for i in fa.values() if len(i)>min_length], dtype=np.int32))[::-1]
    res_nc = fa_lens.shape[0]
    res_totlen = np.sum(fa_lens)
    res_avglen = res_totlen / res_nc
    res_longest = fa_lens[0]

    # Get N25, N50, N75 values:
    fa_lens_cs = np.cumsum(fa_lens)/res_totlen
    ind_n25 = np.where(fa_lens_cs > 0.25)[0][0]
    ind_n50 = np.where(fa_lens_cs > 0.50)[0][0]
    ind_n75 = np.where(fa_lens_cs > 0.75)[0][0]
    res_n25 = fa_lens[ind_n25]
    res_n50 = fa_lens[ind_n50]
    res_n75 = fa_lens[ind_n75]

    return (megahit_dir, min_length, res_nc, res_totlen, res_avglen, res_longest, res_n25, res_n50, res_n75,
            contigs_filename)

def get_megahit_output_contents(rootfold, foldernm):
    '''
    Function for going into a megahit results folder and taking inventory of what's in there.
        1) Checks for the final contigs file present.
            a) If present, gets mtime and size file stats
        2) Checks whether the intermediate_contigs folder has been zipped (i.e. if there is a
            file called 'intermediate_contigs.zip' exists.
        3) Checks whether there is a folder in there called 'intermediate_contigs'
        4) Runs the full get_assembly_summary_stats() routine on the folder.

    Args:
        rootfold (str): path to the output folder where the subfolder was located
        foldernm (str): name of the subfolder we are looking in.

    Returns:
        tsv_values (tuple): results as a tuple, ordered for inserting into a tsv.
        tsv_headers (tuple): headers for the tsv, also ordered.
        new_fields (dict): a full dictionary object containing the results.
    '''
    bioproj, runid = foldernm.split('_')
    conts=os.listdir(os.path.join(rootfold,foldernm))
    is_finished=1 if 'final.contigs.fa' in conts else 0
    int_conts_zipped = 1 if 'intermediate_contigs.zip' in conts else 0
    int_conts_isdir = 1 if os.path.isdir(os.path.join(rootfold,foldernm,'intermediate_contigs')) else 0

    headers = ['bioproj', 'runid', 'is_finished', 'int_conts_zipped', 'int_conts_isdir']
    headers += ['contigs_fa_path', 'contigs_fa_mtime', 'contigs_fa_size']
    mh_ctg_hdrs = get_assembly_summary_stats(os.path.join(rootfold, foldernm), field_names_only=True)
    tsv_headers = tuple(headers) + mh_ctg_hdrs

    new_fields = {'bioproj': bioproj, 'runid': runid, 'conts_folder': conts, 'is_finished': is_finished,
                  'int_conts_zipped': int_conts_zipped, 'int_conts_isdir': int_conts_isdir}

    if is_finished:
        fa_path = os.path.join(rootfold,foldernm,'final.contigs.fa')
        fa_mtime = os.stat(fa_path).st_mtime
        fa_size = os.stat(fa_path).st_size
        mh_ctg_stats = get_assembly_summary_stats(os.path.join(rootfold,foldernm))
        # mh_ctg_stats = dict(zip(mh_ctg_hdrs, mh_ctg_stats))
        # Update dictionary:
        new_fields.update({'contigs_fa_path': fa_path, 'contigs_fa_mtime': fa_mtime, 'contigs_fa_size': fa_size,
                           'megahit_contig_stats': mh_ctg_stats, 'megahit_contig_colheaders': mh_ctg_hdrs})
    else:
        new_fields.update({'contigs_fa_path': None, 'contigs_fa_mtime': None, 'contigs_fa_size': None,
                           'megahit_contig_stats': None})
    #
    tsv_values = [new_fields[hdr] for hdr in headers]
    if is_finished:
        tsv_values = tuple(tsv_values) + mh_ctg_stats
    #
    return tsv_values, tsv_headers, new_fields

def walk_megahit_output_folders(myargs = None):
    '''Goes through the Megahit output folders and makes a TSV with the contents and everything
    needed in general.'''
    mynow = lambda : datetime.datetime.now()
    date_time_8digit=lambda : mynow().strftime('%Y_%m%d_%H%M')
    folds = list(megahit_output_folders.values())

    # backup the old megahit_output_inventory.tsv:
    if os.path.isfile(megahit_output_inventory):
        today = date_time_8digit()
        shutil.copy(megahit_output_inventory, megahit_output_inventory.replace('.tsv', f'_{today}.tsv'))
    mh_out_tsv = open(megahit_output_inventory,'w')
    bioproject_fastq_locs = main.read_bioproject_storage_locs_server()

    def get_generic_tuple_string(num,delim='\t'):
        '''Gets a string that can be used with '''
        return '\t'.join([f'{{{i}}}' for i in range(num)])

    # Walk the folders one at a time:
    megahit_output_vals=[]
    n=1
    for myfold in folds:
        print('Running myfold=%s' % myfold)
        subs=os.listdir(myfold)
        subdirs = [d for d in subs if (os.path.isdir(os.path.join(myfold,d)))]
        d1 = subdirs[0]
        vals, hdrs, _ = get_megahit_output_contents(myfold, d1)
        megahit_output_vals.append(hdrs); megahit_output_vals.append(vals);
        vals_fmtstr = get_generic_tuple_string(len(vals))
        cct = mh_out_tsv.write('\t'.join(hdrs) + '\n')
        cct = mh_out_tsv.write('\t'.join(map(str,vals)) + '\n')

        for d in subdirs[1:]:
            if d.split('_')[0] in bioproject_fastq_locs:
                print(f'{n}: {d}...( {myfold} )')
                try:
                    vals, _, _ = get_megahit_output_contents(myfold, d)
                    megahit_output_vals.append(vals);
                except:
                    logger.info(f'Error occured in subfolder {d}, major_folder: {myfold}')
                cct = mh_out_tsv.write('\t'.join(map(str,vals)) + '\n')
                n+=1

    return megahit_output_vals



if __name__ == '__main__':
    cmd_args = parse_command_args()
    cmd_args.func(cmd_args)
