#!/home/mnute/miniconda3/envs/rice/bin/python

import os, sys, datetime, argparse, logging, glob, shutil
import phylogeny_utilities.utilities as phy
from settings_configs import *
import parsing_utils as parseutils

parser = argparse.ArgumentParser()
rich_format = "[%(filename)s (%(lineno)d) %(asctime)s] %(levelname)s: %(message)s"
logging.basicConfig(format=rich_format, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)

def parse_command_args():
    global parser
    parser.add_argument('-v','--verbose',dest='verbose', action='store_true', help='prints more detailed output to console.')
    parser.add_argument('--debug', dest='debug', action='store_true')

    subparsers = parser.add_subparsers(help='specifies the action to take.')

    fqc_output_p = subparsers.add_parser('fastqc_out_to_tsv',
                                         help='Goes through a root folder containing subfolders with indivdual outputs '
                                              'from fastqc. Each subfolder should be named with the particular bioproject ' 
                                              'it came from. This method will go through each subfolder and find any .zip '
                                              'files, and for each one it will extract the "summary.txt" file, parse it '
                                              'and convert it to a row in a tsv file.')
    fqc_output_p.add_argument('-f','--fastqc_output_folder', dest='fastqc_output_folder', type=str, required=True,
                              help='path to the root folder containing the fastqc outputs.')
    fqc_output_p.add_argument('-o','--output', dest='output_file', type=str, required=True,
                              help='path to the output file to be written to as a tsv.')
    fqc_output_p.add_argument('--no_bioproject_subfolder', action='store_true',
                              help='if given, the method does not assume that the actual outputs are in subfolders named ' 
                                   'by the bioproject. It looks for the .zip files directly in the `fastqc_output_folder` '
                                   'and sets the bioproject column to an empty string.')
    fqc_output_p.add_argument('--no_header', action='store_true', default=False,
                              help='if given, will not print the header to the outputfile (for easier concatenating '
                                   'later).')
    # fqc_output_p.add_argument('-a','--append_mode', dest='append_mode', action='store_true',
    #                           help='if provided, the function will open the output file in append mode and will not '
    #                                'write column headers to start.')
    fqc_output_p.set_defaults(command='fastqc_out_to_tsv', func=fastqc_output_to_tabular)

    # **** FASTQ_INVENTORY_ROUTINE ****
    fastq_count = subparsers.add_parser('inventory_raw_data', help='Routine with go through all the download folders '
                                                                          'and take inventory what has been grabbed and what '
                                                                          'hasnt')
    fastq_count.add_argument('-o', '--output', dest='output_file', type=str, required=False, default=resolve_file(fastq_gz_inventory),
                              help='path to the output file to be written to as a tsv. (default = %s)' % resolve_file(fastq_gz_inventory))
    fastq_count.add_argument('-li', '--last_inventory', dest='last_inventory', type=str, required=False, default='',
                             help='path to the file containing the last inventory produced. Default is the generic inventory '
                                  'path pointed to in the settings_configs.py file.')
    fastq_count.set_defaults(func=take_file_inventory)

    # **** POST-PROCESSING_FILE_INVENTORY ****
    postproc_inventory = subparsers.add_parser('postproc_inventory',
                                               help = 'Routine to go through the fastq file inventory and create a '
                                                      'separate file with the status of all the files for each of the '
                                                      'post-processing steps.')
    postproc_inventory.set_defaults(func=take_postproc_file_inventory)

    # **** BBDUK-DIFFS FILE MAKER ****
    bbduk_diffs = subparsers.add_parser('bbduk_make_diffs_file',
                                        help='Routine to take an original fastq and a post-BBduk fastq and identify all '
                                             'reads with differences, then write those reads from the original to a new '
                                             '`diffs` file.')
    bbduk_diffs.add_argument('-i', '--orig', dest='orig_fastq', type=str, required=True,
                             help='path to the original fastq file. If the argument \'--input_as_file_list\' is '
                                  'provided then this should be the path to a .txt file with three tab-separated '
                                  'columns containing the original-file, bbduk-file, and diffs-file paths '
                                  '(respectively) on each line. If that is given then the \'-b\' and \'-d\' options '
                                  'are ignored.')
    bbduk_diffs.add_argument('-b', '--bbduk', dest='bbduk_fastq', type=str, default='-',
                             help='path to the post-BBduck fastq file.')
    bbduk_diffs.add_argument('-o', '--output', dest='output_fastq', type=str, default='-',
                             help='path to the output diffs fastq file.')
    bbduk_diffs.add_argument('-q', '--quiet', dest='quiet_mode', action='store_true', default=False,
                             help='If given, don\'t print the verbose outputs that are printed by default.')
    bbduk_diffs.add_argument('--input_as_file_list', dest='input_as_file_list', action='store_true', default=False,
                             help='If given, don\'t print the verbose outputs that are printed by default.')
    bbduk_diffs.set_defaults(func=run_make_bbduk_diffs_file)

    bbduk_diffs_test = subparsers.add_parser('bbduk_test_diffs_file',
                                             help='Routine to test that a post-BBduks file plus a bbduk-diffs file '
                                                  'can be put together to make all the same reads as in the original '
                                                  'fastq file.')
    bbduk_diffs_test.add_argument('-i', '--orig', dest='orig_fastq', type=str, required=True)
    bbduk_diffs_test.add_argument('-b', '--bbduk', dest='bbduk_fastq', type=str, required=True)
    bbduk_diffs_test.add_argument('-d', '--diffs', dest='diffs_fastq', type=str, required=True)
    bbduk_diffs_test.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=False)
    bbduk_diffs_test.set_defaults(func=run_test_bbduk_diffs_file)

    bbduk_diffs_delete_batch = subparsers.add_parser('bbduk_diffs_delete_batch',
                                            help = 'Takes a space-separated list of Accessions and a single Bioproject '
                                                   'and runs the diffs PLUS the original fastq deletion. Only to be run '
                                                   'after all the modules are trusted. This subcommand uses the python '
                                                   'argparse modules `nargs` argument to convert a space delimited list '
                                                   'into a list of strings.')
    bbduk_diffs_delete_batch.add_argument('-bp', '--bioproject', dest='bioproject', type=str, required=True)
    bbduk_diffs_delete_batch.add_argument('-ri', '--runids', dest='runids', nargs='+', required=True)
    bbduk_diffs_delete_batch.add_argument('--no_delete', dest='no_delete', action='store_true', default=False)
    bbduk_diffs_delete_batch.add_argument('--fake_run', dest='fake_run', action='store_true', default=False)
    bbduk_diffs_delete_batch.set_defaults(func=batch_bbduk_diffs_delete)

    bbduk_get_logdata = subparsers.add_parser('bbduk_get_logdata',
                                          help = 'Iterates through a large set of bbduk logfiles (typcially all dumped '
                                                 'to a single fodler, but can also be provided as a python list object '
                                                 'containing valid paths.')
    bbduk_get_logdata.add_argument('-i','--input', dest='input', type=str, required=True,
                                   help = 'The path to the folder storing all the bbduk logs. To be parsed file '
                                          'by file. ')
    bbduk_get_logdata.add_argument('-o', '--output', dest='output', type=str, required=True,
                                   help='The path to the file that the output should be written to.')
    bbduk_get_logdata.set_defaults(func=gather_bbduks_results_from_bulk_folder)

    args = parser.parse_args()
    command_args_postprocess(args)
    return args

def command_args_postprocess(cmd_args):
    '''Misc things that should be done as soon as the args are parsed regardless of what's in them.'''
    if cmd_args.verbose:
        logger.setLevel(logging.INFO)
    if cmd_args.debug:
        logger.setLevel(logging.DEBUG)

def fastqc_output_to_tabular(args):
    '''
    Parses the output folders for a batch fastqc run and summarizes them to a tab-delimited file. Target folder should
    contain individual subfolders for each bioproject, though there is an option that doesn't require that.

    The fastqc results are summarized in two particular files within the .zip file that it outputs: summary.txt,
        and fastqc_data.txt. This function goes through all the fastqc result files a) pulls those two files out
        from the .zip, and b) parses them and puts them back with nothing new written to disk. Finally, it takes
        the outputs from both of these and converts it to being tab-delimited which can be explored in Excel.

    REQUIRED COMMAND-LINE ARGUMENTS:
        -f: path to the fastqc_output folder we should summarize.
        -o: output file containing the data in tab-delimited text form.

    '''
    import zipfile
    fastqc_output_folder = args.fastqc_output_folder
    output_tsv_file = args.output_file
    no_bp_subfolder = args.no_bioproject_subfolder
    omit_header = args.no_header
    # append_mode = args.append_mode

    # helper function to parse the summary.txt files:
    def parse_summary_txt(zfobj: zipfile.ZipFile):
        '''
        Parses the summary.txt file within the zipfile.

        Returns a dict-of-dicts in the form:
            { <file_1>: { <field_1>: <pass-warn-fail>, <field-2>: ...}, }

        for as many files as are listed in the summary file.
        '''
        name = [i for i in zfobj.namelist() if i.find('summary.txt')>-1][0]
        summ = [i.decode('utf-8').strip().split('\t') for i in zfobj.open(name).readlines()]
        files_included = list(set([i[2] for i in summ]))
        res = {f: {} for f in files_included}
        for i in summ:
            res[i[2]][i[1]] = i[0]
        return res

    # helper function to parse the fastqc_data.txt files (headers only):
    def parse_fastqc_data_txt(zfobj: zipfile.ZipFile):
        '''parses the file-level header fields in the `fastqc_data.txt` file.'''
        name = [i for i in zfobj.namelist() if i.find('fastqc_data.txt') > -1][0]
        out_args = {}
        with zfobj.open(name) as fdata:
            lns = [i.decode('utf-8').strip() for i in fdata.readlines()]
        stat_lns = [i for i in range(len(lns)) if lns[i].startswith('>>Basic Statistics')]
        for ln_num in stat_lns:
            i = 0; fields=[];
            while lns[ln_num + i].find('END_MODULE')==-1:
                if lns[ln_num+i][0]=='#' or lns[ln_num+i][:2]=='>>':
                    i+=1; continue;
                fields.append(lns[ln_num+i].split('\t'))
                i += 1
            field_dict = dict(fields)
            filename = field_dict.pop('Filename')
            out_args[filename]=field_dict
        return out_args

    # Get an inventory of the subfolders to be covered. (namely bioprojects):
    if no_bp_subfolder:
        subfolders = ['']
    else:
        subfolders = [i for i in os.listdir(fastqc_output_folder) if os.path.isdir(os.path.join(fastqc_output_folder,i))]
    file_list = {}
    for sf in subfolders:   # Collect the zip files in all subfolders into a list of lists
        sfpa = os.path.join(fastqc_output_folder, sf)
        zlist = [i for i in os.listdir(sfpa) if i[-4:]=='.zip']
        file_list[sf] = zlist


    # Workhorse loop: go through each zip file, parse the key outputs and store the results.
    n_zip_files = sum([len(i) for i in file_list.values()])
    logger.info(f'Found {n_zip_files} zip files. Parsing them one by one...')

    all_output_vals = {}; summ_hdrs = set(); data_hdrs = set(); file_num = 1;
    for subfolder in file_list:
        all_output_vals[subfolder] = {}
        for zf_path in file_list[subfolder]:
            zf_fullpath = os.path.join(fastqc_output_folder, subfolder, zf_path)
            with zipfile.ZipFile(zf_fullpath) as zf:
                # storing each of the two files in their own variable for now. Will fixe later.
                summ_vals = parse_summary_txt(zf)
                data_vals = parse_fastqc_data_txt(zf)
                # summ_hdrs.update(summ_vals); data_hdrs.update(data_vals);
            all_output_vals[subfolder][zf_path] = (summ_vals, data_vals)
            print(f'...done with file {file_num} of {n_zip_files}: {subfolder}/{zf_path} ({phy.mynowstr()}).', end='\r')
            file_num += 1

    # Make a single dictionary of summary_headers and data_headers:
    for sf in all_output_vals:
        for zf in all_output_vals[sf]:
            if set(all_output_vals[sf][zf][0].keys())!=set(all_output_vals[sf][zf][1].keys()):
                my_zpath = os.path.join(sf,zf); sfqs=str(set(all_output_vals[sf][zf][0].keys()));
                dfqs=set(all_output_vals[sf][zf][1].keys());
                logger.warning(f'In zip file {my_zpath}, summary/fastqc_data filenames do not overlap: summary.txt={sfqs}, fastqc_data.txt={dfqs}')
            for fq in all_output_vals[sf][zf][0]:
                summ_hdrs.update(all_output_vals[sf][zf][0][fq].keys())
            for fq in all_output_vals[sf][zf][1]:
                data_hdrs.update(all_output_vals[sf][zf][1][fq].keys())

    # Finally we are done with the parsing, now comes the writing portion of practice.
    print(''); logger.info('Done parsing ZipFiles. Writing output file: %s' % output_tsv_file);
    # Set up output file:
    file_loc_hdrs = ['root_folder', 'subfolder', 'zip_file', 'zip_mtime', 'source_file']
    out_f = open(output_tsv_file, 'w')


    # headers = file_loc_hdrs[1:] + list(summ_hdrs) + list(data_hdrs) + file_loc_hdrs[:1]
    headers = ['subfolder', 'zip_file', 'zip_mtime', 'source_file', 'Basic Statistics', 'Per sequence GC content',
     'Overrepresented sequences', 'Per tile sequence quality', 'Per sequence quality scores',
     'Per base sequence quality', 'Per base sequence content', 'Sequence Length Distribution', 'Adapter Content',
     'Per base N content', 'Sequence Duplication Levels', 'Sequences flagged as poor quality', 'File type', 'Encoding',
     'Sequence length', '%GC', 'Total Sequences', 'Total Bases', 'root_folder']
    # This line doesn't look pretty but it works, The main headline for the rest of the way is that we are
    #   trying to combine a bunch of different data types from the data-headers and software.
    dummy_line = '\t'.join([f'{{{i}}}' for i in headers]) + '\n' # triple {{{ and }}} escape the f-string curlys
    logger.debug('headers = ' + str(headers))
    logger.debug('dummy_line = \'' + dummy_line.strip() + '\'')
    if not omit_header:
        cct = out_f.write('\t'.join(headers) + '\n')
    for subfolder in all_output_vals:
        for zf_path in all_output_vals[subfolder]:
            mod_time = datetime.datetime.fromtimestamp(os.stat(os.path.join(fastqc_output_folder, subfolder, zf_path)).st_mtime)
            orig_files_combined = set(all_output_vals[subfolder][zf_path][0].keys())
            orig_files_combined.update(set(all_output_vals[subfolder][zf_path][1].keys()))
            for orig_file in orig_files_combined:
                vals = dict(zip(file_loc_hdrs, [fastqc_output_folder, subfolder, zf_path, phy.my_dt_format(mod_time), orig_file]))
                s_vals = {i: '' for i in summ_hdrs}
                s_vals.update(all_output_vals[subfolder][zf_path][0].get(orig_file, {}))
                d_vals = {i: '' for i in data_hdrs}
                d_vals.update( all_output_vals[subfolder][zf_path][1].get(orig_file, {}) )
                vals.update(s_vals) # summary values
                vals.update(d_vals) # data values
                cct = out_f.write(dummy_line.format(**vals))

    logger.info('Done!')

def read_bioproject_storage_locs_server():
    '''Gets the list of locations for each BioProject.'''
    # bp_locs=shared_paths['bioproject_fastq_locs']['server']
    bp_locs = bioproject_fastq_locs['server']
    return phy.dict_from_tab_delimited_file(bp_locs)

def take_postproc_file_inventory(cmd_args):
    logger.info(f'Running post-processing file inventory routine...')

    out_file = resolve_file(postproc_file_inventory)
    logger.info(f'   output file path: {out_file}')

    args = [
        'bioproj',      'accn',
        'fq_folder_path',
        'fq_file_1',    'fq_file_1_exists', 'fq_file_2',    'fq_file_2_exists',
        'fqc_folder_path',
        'fqc_zip_1_name', 'fqc_zip_1_mtime', 'fqc_zip_2_name', 'fqc_zip_2_mtime',
        'bbduk_folder_path',
        'bbduk_f1_name', 'bbduk_f1_size', 'bbduk_f1_mtime',
        'bbduk_f2_name', 'bbduk_f2_size', 'bbduk_f2_mtime',
        'bbduk_diffs_folder_path',
        'bbduk_diffs_f1_name', 'bbduk_diffs_f1_size', 'bbduk_diffs_f1_mtime',
        'bbduk_diffs_f2_name', 'bbduk_diffs_f2_size', 'bbduk_diffs_f2_mtime',
    ]
    def template_line():
        '''helper function to make a nice string that we can do
        str.format(**args) with to get the line formatted neatly.'''
        return '\t'.join(['{' + i + '}' for i in args]) + '\n'

    fastq_inventory_path = resolve_file(fastq_gz_inventory)
    logger.info(f'    fastq file inventory path: {fastq_inventory_path}')
    proj2path = read_bioproject_storage_locs_server()

    last_fastq_inventory = read_last_fastq_inventory_file(fastq_inventory_path)
    postproc_inventory_records = []
    for bp, accn in last_fastq_inventory:
        bioproj_folder = os.path.join(proj2path[bp], bp)
        fastqc_folder = bioproj_folder.replace('rawdata', 'fastqc_output')
        bbduk_folder = bioproj_folder.replace('rawdata', 'post_bbduk_data')
        bbduk_diffs_folder = bioproj_folder.replace('rawdata', 'bbduk_diffs')

        rec = {i: '' for i in args}
        rec['bioproj'] = bp; rec['accn'] = accn;
        rec['fq_folder_path'] = bioproj_folder
        rec['fq_file_1'] = last_fastq_inventory[(bp, accn)]['file_1']
        rec['fq_file_2'] = last_fastq_inventory[(bp, accn)]['file_2']
        rec['fq_file_1_exists'] = os.path.isfile(os.path.join(bioproj_folder, rec['fq_file_1']))
        rec['fq_file_2_exists'] = os.path.isfile(os.path.join(bioproj_folder, rec['fq_file_2']))
        # FastQC Output Metadata:
        rec['fqc_folder_path'] = fastqc_folder
        fqc_file_paths_1 = glob.glob(os.path.join(fastqc_folder, f'{accn}_1*.zip'))
        fqc_file_paths_2 = glob.glob(os.path.join(fastqc_folder, f'{accn}_2*.zip'))
        if len(fqc_file_paths_1)>0:
            rec['fqc_zip_1_name'] = os.path.split(fqc_file_paths_1[0])[1]
            rec['fqc_zip_1_mtime'] = timestamp_to_str(os.stat(fqc_file_paths_1[0]).st_mtime)
        if len(fqc_file_paths_2)>0:
            rec['fqc_zip_2_name'] = os.path.split(fqc_file_paths_2[0])[1]
            rec['fqc_zip_2_mtime'] = timestamp_to_str(os.stat(fqc_file_paths_2[0]).st_mtime)
        # BBduk Results Metadata:
        rec['bbduk_folder_path'] = bbduk_folder
        bbduk_file_paths_1 = glob.glob(os.path.join(bbduk_folder, f'{accn}_1*'))
        bbduk_file_paths_2 = glob.glob(os.path.join(bbduk_folder, f'{accn}_2*'))
        if len(bbduk_file_paths_1)>0:
            rec['bbduk_f1_name'] = os.path.split(bbduk_file_paths_1[0])[1]
            rec['bbduk_f1_size'] = os.stat(bbduk_file_paths_1[0]).st_size
            rec['bbduk_f1_mtime'] = timestamp_to_str(os.stat(bbduk_file_paths_1[0]).st_mtime)
        if len(bbduk_file_paths_2)>0:
            rec['bbduk_f2_name'] = os.path.split(bbduk_file_paths_2[0])[1]
            rec['bbduk_f2_size'] = os.stat(bbduk_file_paths_2[0]).st_size
            rec['bbduk_f2_mtime'] = timestamp_to_str(os.stat(bbduk_file_paths_2[0]).st_mtime)
        # BBduk Diffs Metadata:
        rec['bbduk_diffs_folder_path'] = bbduk_diffs_folder
        bbduk_diffs_file_paths_1 = glob.glob(os.path.join(bbduk_diffs_folder, f'{accn}_1*'))
        bbduk_diffs_file_paths_2 = glob.glob(os.path.join(bbduk_diffs_folder, f'{accn}_2*'))
        if len(bbduk_diffs_file_paths_1) > 0:
            rec['bbduk_diffs_f1_name'] = os.path.split(bbduk_diffs_file_paths_1[0])[1]
            rec['bbduk_diffs_f1_size'] = os.stat(bbduk_diffs_file_paths_1[0]).st_size
            rec['bbduk_diffs_f1_mtime'] = timestamp_to_str(os.stat(bbduk_diffs_file_paths_1[0]).st_mtime)
        if len(bbduk_diffs_file_paths_2) > 0:
            rec['bbduk_diffs_f2_name'] = os.path.split(bbduk_diffs_file_paths_2[0])[1]
            rec['bbduk_diffs_f2_size'] = os.stat(bbduk_diffs_file_paths_2[0]).st_size
            rec['bbduk_diffs_f2_mtime'] = timestamp_to_str(os.stat(bbduk_diffs_file_paths_2[0]).st_mtime)
        # append it...
        postproc_inventory_records.append(rec)

    # Write all the records to the output file:
    with open(out_file,'w') as of:
        template = template_line()
        cct = of.write(template.format(**{i: i for i in args}))
        ln_ct = 0
        for r in postproc_inventory_records:
            cct = of.write(template.format(**r))
            ln_ct += 1

    logger.info(f'   Done. Write {ln_ct} lines to the output file.')

def read_last_fastq_inventory_file(last_inventory_file=None):
    '''Get a dictionary of all the data in the previous inventory file so that we can account for files that
    were previously downloaded but have since been deleted. Returns a dict keyed by (bioproject, accession)
    where each item is a dictionary with the keys given by the file column headers (so ready to go right into
    the string-format print statment in the new file).'''
    last_inventory_file = resolve_file(fastq_gz_inventory) if last_inventory_file is None else last_inventory_file
    lastinv = {}
    with open(last_inventory_file,'r') as lf:
        last_hdrs = lf.readline().strip().split('\t')
        for inv_line in lf:
            inv_items = inv_line.strip().split('\t')
            # { (bioproject, accession): { <header_1>: <field_1>, ...}, ... }
            lastinv[(inv_items[0], inv_items[1])]=dict(zip(last_hdrs, inv_items))
    return lastinv

def timestamp_to_str(myts):
    '''Helper to  hurry along the timestamp conversion later.

    Args:
        myts (int): linux time stamp as of the form from os.stat.(<file_path>).st_mtime

    Returns:
        myts_string (str): same time-stamp formatted as a nice string as 'YYYY-MM-DD HH:MM:SS'
    '''
    return phy.my_dt_format(datetime.datetime.fromtimestamp(myts))

def take_file_inventory(cmd_args):
    '''
    Go through and get the basic file inventory and metadata, and organize it so that the 1/2 files for each individual
    SRA run are put together where possible.

    No arguments are technically necessary for this. All files should be consisent between runs.

    Args:
        cmd_args:

    Returns:

    '''
    proj2path = read_bioproject_storage_locs_server()
    if cmd_args.last_inventory!='' and os.path.isfile(cmd_args.last_inventory):
        logger.info(f'last inventory pulling from user-supplied path: {cmd_args.last_inventory}')
        last_inventory = read_last_fastq_inventory_file(cmd_args.last_inventory)
    else:
        logger.info(f'last inventory pulling from default path: {resolve_file(fastq_gz_inventory)}')
        last_inventory = read_last_fastq_inventory_file()

    md_args = {'bioproj': '', 'accn': '', 'file_1': '', 'file_2': '',
                'size_1': '', 'atime_1ts': '', 'mtime_1ts': '',
                'ctime_1ts': '', 'atime_1': '', 'mtime_1': '',
                'ctime_1': '','size_2': '', 'atime_2ts': '',
                'mtime_2ts': '', 'ctime_2ts': '', 'atime_2': '',
                'mtime_2': '', 'ctime_2': '',
                'file_1_removed': '', 'file_1_last_seen': '',
                'file_2_removed': '', 'file_2_last_seen': '',}
    hdrs = list(md_args.keys())
    # will use str.format(...) later so this will help:
    def get_md_template_line():
        '''helper function to make a nice string that we can do
        str.format(**args) with to get the line formatted neatly.'''
        return '\t'.join(['{' + i + '}' for i in hdrs]) + '\n'

    # Output file name (should use the default but just in case)
    out_file_name = cmd_args.output_file

    output = []
    # Need a list of all the (bioproject,accession) pairs that we come across in the file sweep so
    #   that we can add anything that has been fully deleted to the inventory from the old one.
    bioproj_accn_fastq_exists = []
    for bioproj in proj2path:
        target_folder = os.path.join(proj2path[bioproj], bioproj)

        if not os.path.isdir(target_folder):
            logger.warning(f'targeted folder does not exist {target_folder}')
            continue
        files = [i for i in os.listdir(target_folder) if i.find('.fastq.gz')>-1]
        accns = list(set([i.split('_')[0] for i in files]))

        for a in accns:
            bioproj_accn_fastq_exists.append((bioproj, a))
            rec = md_args.copy()
            rec['bioproj'] = bioproj
            rec['accn'] = a
            if (bioproj,a) in last_inventory:
                rec.update(last_inventory[(bioproj,a)])
            fwd_filename = f'{a}_1.fastq.gz'
            if fwd_filename in files:
                fpa = os.path.join(target_folder, fwd_filename)
                fstat = os.stat(fpa)
                rec['file_1'] = fwd_filename
                rec['size_1'] = fstat.st_size
                rec['atime_1ts'] = fstat.st_atime
                rec['mtime_1ts'] = fstat.st_mtime
                rec['ctime_1ts'] = fstat.st_ctime
                rec['atime_1'] = timestamp_to_str(fstat.st_atime)
                rec['mtime_1'] = timestamp_to_str(fstat.st_mtime)
                rec['ctime_1'] = timestamp_to_str(fstat.st_ctime)
                rec['file_1_removed'] = False
                rec['file_1_last_seen'] = phy.mynowstr()
            else:
                if (bioproj, a) in last_inventory and last_inventory[(bioproj,a)]['file_1']!='':
                    rec['file_1_removed'] = True
            rev_filename = fwd_filename.replace('_1.fastq.gz','_2.fastq.gz')
            if rev_filename in files:
                fpa = os.path.join(target_folder, rev_filename)
                fstat = os.stat(fpa)
                rec['file_2'] = rev_filename
                rec['size_2'] = fstat.st_size
                rec['atime_2ts'] = fstat.st_atime
                rec['mtime_2ts'] = fstat.st_mtime
                rec['ctime_2ts'] = fstat.st_ctime
                rec['atime_2'] = timestamp_to_str(fstat.st_atime)
                rec['mtime_2'] = timestamp_to_str(fstat.st_mtime)
                rec['ctime_2'] = timestamp_to_str(fstat.st_ctime)
                rec['file_2_removed'] = False
                rec['file_2_last_seen'] = phy.mynowstr()
            else:
                if (bioproj, a) in last_inventory and last_inventory[(bioproj,a)]['file_2']!='':
                    rec['file_2_removed'] = True
            # add it to the list...
            output.append(rec)
    #
    # Write everything to the output file...
    with open(out_file_name, 'w') as my_out_f:
        template = get_md_template_line()
        cct = my_out_f.write(template.format(**dict(zip(hdrs, hdrs))))
        for i in output:
            cct = my_out_f.write(template.format(**i))
    #
    return

def make_bbduk_diffs_file(orig_path, bbduk_path, diffs_out_path, verbose=True, out_as_gzip=True):
    '''
    Reads the original and the post-BBduk fastq files and makes a fastq file containing all and only
    the reads that are **NOT** identical between the two. Also if for any reason some reads are in the
    post-BBduk file that don't match any in the original, adds these at the end with the string `-NEW`
    appended to the readname.
    '''
    start = phy.mynow()
    if verbose:
        print(f'Creating Diffs FastQ file  ({phy.mynowstr()})')
        print(f'    Original:   {orig_path}')
        print(f'    Post-BBduk: {bbduk_path}')

    # 1) Read both input files
    fabb=phy.fastq_read_from_file_simple(bbduk_path)
    faor=phy.fastq_read_from_file_simple(orig_path)

    if verbose:
        print(f'Done reading inputs...     ({phy.mynowstr()})')

    # 2) Compute the set of missing and changed reads:
    n_diff = len(faor) - len(fabb)
    ks_or_only = list(set(faor.keys()) - set(fabb.keys()))
    ks_bb_only = list(set(fabb.keys()) - set(faor.keys()))
    ks_overlap = list(set(fabb.keys()) & set(faor.keys()))
    n_overlap = len(ks_overlap)
    n_or_only = len(ks_or_only)
    n_bb_only = len(ks_bb_only)
    ks_not_equal = [i for i in ks_overlap if fabb[i]!=faor[i]] #print(phy.mynowstr());
    n_not_equal = len(ks_not_equal)

    if verbose:
        print(f'Computed read overlaps...  ({phy.mynowstr()})')
        print(f'    Original # reads: {len(faor):,},  PostBB # reads: {len(fabb):,}  (diff={n_diff:,})')
        print(f'    # overlap: {n_overlap},  # Orig only: {n_or_only},  # BB only: {n_bb_only}')
        print(f'    # overlapping reads not equal: {n_not_equal:,}')

    # 3) Write them to the diffs file.
    lns = ['\n'.join((f'@{k}',)+faor[k]) + '\n' for k in (ks_or_only + ks_not_equal)]
    lns += ['\n'.join((f'@{k}-NEW',)+faor[k]) + '\n' for k in ks_bb_only]
    if not out_as_gzip:
        with open(diffs_out_path,'w') as outf:
            for ln in lns:
                cct = outf.write(ln)
    else:
        diffs_out_path = diffs_out_path if diffs_out_path[-3:]=='.gz' else diffs_out_path + '.gz'
        phy.gz_write_from_text_lines(diffs_out_path, lns)

    end_time = phy.mynow()
    if verbose:
        print(f'Wrote diffs file...        ({phy.mynowstr()})')
        print(f'    total time elapsed:     {end_time - start}')


def test_bbduk_diffs_orig_equivalence(orig_path, bbduk_path, diffs_out_path, verbose=True):
    '''Tests whether the original fastq is equal to the bbduk one plus the diffs (minus
    any with '-NEW' at the end of the readname'''
    assert os.path.isfile(orig_path), f'ERROR: file does not exist:{orig_path}'
    assert os.path.isfile(bbduk_path), f'ERROR: file does not exist:{bbduk_path}'
    assert os.path.isfile(diffs_out_path), f'ERROR: file does not exist:{diffs_out_path}'
    start_time = phy.mynow()
    if verbose:
        print(f'Reading input files...           ({phy.mynowstr()})')

    fabb = phy.fastq_read_from_file_simple(bbduk_path, verbose=False)
    if verbose:
        print(f'    done reading bbduk file. elapsed: {phy.mynow() - start_time}')
    faor = phy.fastq_read_from_file_simple(orig_path, verbose=False)
    if verbose:
        print(f'    done reading orig file. elapsed:  {phy.mynow() - start_time}')
    fadi = phy.fastq_read_from_file_simple(diffs_out_path, verbose=False)
    if verbose:
        print(f'    done reading diffs file. elapsed: {phy.mynow() - start_time}')

    if verbose:
        print(f'    ...done reading files        ({phy.mynowstr()})')
        print(f'    elapsed: {phy.mynow() - start_time}')

    fabb.update(fadi)
    di_ks_new = [k for k in fadi.keys() if k[-4:] == '-NEW']
    if len(di_ks_new) > 0:
        for k in di_ks_new:
            del fabb[k]
            del fabb[k[:-4]]

    assert len(faor) == len(fabb), f'Dict lengths post-merge are not the same: orig={len(faor)}, bb={len(fabb)}'
    assert len(faor) == len(set(fabb.keys()) & set(faor.keys())), f'Dict key sets are not equal'
    ks_not_equal = [k for k in faor.keys() if faor[k] != fabb[k]]
    assert len(ks_not_equal) == 0, f'Not all read/quality-scores are equal between the two (# diffs =  {len(ks_not_equal)})'
    if verbose:
        print(f'Done! Files {bbduk_path} and diffs {diffs_out_path} can recreate {orig_path} successfully')
        print(f'Time Elapsed: {phy.mynow() - start_time}')
    else:
        print('True')

def batch_bbduk_diffs_delete(cmd_args):
    '''Takes the bioproject and the list of accessions. For each one, figures out where the various files should be
    located, and if everything is in order it goes ahead and computes the bbduk-diffs file, and finally it deletes
    the original.'''
    no_del = cmd_args.no_delete
    bioproj = cmd_args.bioproject
    accns = cmd_args.runids
    fake_run = cmd_args.fake_run
    logger.info(f'Running a batch bbduks run with:')
    logger.info(f'    bioproject = {bioproj}')
    logger.info(f'    accessions = {len(accns)} total. First 5: {str(accns[:5])} ...')
    logger.info(f'    Delete originals = {not no_del}')
    logger.info(f'    Fake Run? = {fake_run}')

    proj2path = read_bioproject_storage_locs_server()
    rawdata_fold = os.path.join(proj2path[bioproj], bioproj)
    bbduk_fold = rawdata_fold.replace('rawdata','post_bbduk_data')
    diffs_fold = rawdata_fold.replace('rawdata', 'bbduk_diffs')
    assert os.path.isdir(rawdata_fold), f'raw data folder {rawdata_fold} does not exist...'
    assert os.path.isdir(bbduk_fold), f'raw data folder {bbduk_fold} does not exist...'
    assert os.path.isdir(diffs_fold), f'raw data folder {diffs_fold} does not exist...'

    for accn in accns:
        raw_fwd = glob.glob(os.path.join(rawdata_fold, f'{accn}_1*'))
        raw_rev = glob.glob(os.path.join(rawdata_fold, f'{accn}_2*'))
        try:
            assert len(raw_fwd)==1 and len(raw_rev)==1, f'Number of Raw data files for fwd ({len(raw_fwd)}) and rev ({len(raw_rev)}) ' \
                                                    f'not both equal to 1.'
        except:
            logger.error(f'Accession {accn} encountered an Error because the raw glob files were not both eqaul to 1.')
            logger.error(f'    # Files: fwd ({len(raw_fwd)}), rev ({len(raw_rev)})')
            for i in raw_fwd + raw_rev:
                logger.error(f'      {i}')
            continue
        raw_fwd = raw_fwd[0]; raw_rev = raw_rev[0];
        bbduk_fwd = glob.glob(os.path.join(bbduk_fold, f'{accn}_1*'))
        bbduk_rev = glob.glob(os.path.join(bbduk_fold, f'{accn}_2*'))
        try:
            assert len(bbduk_fwd) == 1 and len(bbduk_rev) == 1, \
            f'Number of BBduk_output data files for fwd ({len(bbduk_fwd)}) and rev ({len(bbduk_rev)}) not both size 1.'
        except:
            logger.error(f'Accession {accn} encountered an Error because the bbduk glob files were not both size 1.')
            logger.error(f'    # Files: fwd ({len(bbduk_fwd)}), rev ({len(bbduk_rev)})')
            for i in bbduk_fwd + raw_rev:
                logger.error(f'      {i}')
            continue
        bbduk_fwd = bbduk_fwd[0]; bbduk_rev = bbduk_rev[0];
        # Get file paths for diffs file that we will create
        diffs_fwd = os.path.join(diffs_fold, os.path.split(raw_fwd)[1].replace('.fastq.gz', '_bbduk_diffs.fastq.gz'))
        diffs_rev = os.path.join(diffs_fold, os.path.split(raw_rev)[1].replace('.fastq.gz', '_bbduk_diffs.fastq.gz'))

        print('\n****************************************************************************************************\n')
        logger.info(f'Making diffs file for Fwd reads for BP={bioproj}, Accn={accn}:')
        logger.info(f'    raw_fwd={raw_fwd}')
        logger.info(f'    bbduk_fwd={bbduk_fwd}')
        logger.info(f'    diffs_fwd={diffs_fwd}\n\n')
        if fake_run:
            logger.info(f'    (FAKE RUN: This is where we would make the fwd diffs file if this were real.)')
        else:
            make_bbduk_diffs_file(raw_fwd, bbduk_fwd, diffs_fwd)
        logger.info(f'   Done.')
        if not no_del:
            if fake_run:
                logger.info(f'    (FAKE RUN: This is where we would delete the fwd rawdata file {raw_fwd}, if this were real.)')
            else:
                logger.info(f'    Deleting {raw_fwd}')
                os.remove(raw_fwd)

        logger.info(f'Making diffs file for Rev reads for BP={bioproj}, Accn={accn}:')
        logger.info(f'    raw_fwd={raw_rev}')
        logger.info(f'    bbduk_fwd={bbduk_rev}')
        logger.info(f'    diffs_fwd={diffs_rev}\n\n')
        if fake_run:
            logger.info(f'    (FAKE RUN: This is where we would make the fwd diffs file if this were real.)')
        else:
            make_bbduk_diffs_file(raw_rev, bbduk_rev, diffs_rev)
        logger.info(f'   Done.')
        if not no_del:
            if fake_run:
                logger.info(
                    f'    (FAKE RUN: This is where we would deleted the reverse rawdata file {raw_rev}, if this were real.)')
            else:
                logger.info(f'    Deleting {raw_rev}')
                os.remove(raw_rev)
        logger.info('\n');
        logger.info(f'Done with accession {accn}!')



def run_make_bbduk_diffs_file(cmd_args):
    '''Dispatcher for the make-file function `make_bbduk_diffs_file(...)`.'''
    verb = not cmd_args.quiet_mode
    if not cmd_args.input_as_file_list:
        assert os.path.isfile(cmd_args.orig_fastq)
        assert os.path.isfile(cmd_args.bbduk_fastq)
        assert os.path.isdir(os.path.split(cmd_args.output_fastq)[0])
        logger.info(f'Running file-combination:')
        logger.info(f'    orig:  {cmd_args.orig_fastq}')
        logger.info(f'    bbduk: {cmd_args.bbduk_fastq}')
        logger.info(f'    diffs: {cmd_args.output_fastq}')
        make_bbduk_diffs_file(cmd_args.orig_fastq, cmd_args.bbduk_fastq, cmd_args.output_fastq, verbose=verb)
    else:
        file_path_list = [tuple(i.split('\t')) for i in phy.get_list_from_file(cmd_args.orig_fastq)]
        logger.info(f'Running with input as a file-path combination list with {len(file_path_list)} entries.')
        logger.info(f'    ...input list file: {cmd_args.orig_fastq}')
        ct = 1
        for orig, bbduk, diffs in file_path_list:
            assert os.path.isfile(orig)
            assert os.path.isfile(bbduk)
            assert os.path.isdir(os.path.split(diffs)[0])
            logger.info(f'Running file-combination {ct}:')
            logger.info(f'    orig:  {orig}')
            logger.info(f'    bbduk: {bbduk}')
            logger.info(f'    diffs: {diffs}')
            make_bbduk_diffs_file(cmd_args.orig_fastq, cmd_args.bbduk_fastq, cmd_args.output_fastq, verbose=verb)
            ct += 1

def run_test_bbduk_diffs_file(cmd_args):
    '''dispatcher for the function `test_bbduk_diffs_orig_equivalence`'''
    verb = not cmd_args.quiet
    isgood = test_bbduk_diffs_orig_equivalence(cmd_args.orig_fastq, cmd_args.bbduk_fastq, cmd_args.diffs_fastq, verbose=verb )

def gather_bbduks_results_from_bulk_folder(cmd_args):
    i_path = cmd_args.input
    o_path = cmd_args.output
    logger.info(f'Running bbduk logfile gathering on target folder: {i_path}')
    assert os.path.isdir(i_path), f'The provided input with `-i` was not a valid folder: {i_path}'
    o_path_dir, o_path_name = os.path.split(o_path)
    assert os.path.isdir(o_path_dir), f'The provided input does not split into a valid directory: {o_path} became ' \
                                      f'{o_path_dir} using os.path.split().'
    str_out = parseutils.parse_bbduk_logfile_folder_to_tsv(i_path)
    #
    with open(o_path, 'w') as out_f:
        cct = out_f.write(str_out)
    #
    logger.info(f'...done reading all log files. Wrote final output to specified output path: {o_path}')

if __name__ == '__main__':
    cmd_args = parse_command_args()
    cmd_args.func(cmd_args)
