
import os, sys, datetime, argparse, logging, glob, re, subprocess
import phylogeny_utilities.utilities as phy
from settings_configs import *


def parse_bbduk_logfile(logfile_path=None, get_descriptions_only=False, get_groupkeys_only=False):
    '''Takes a bbduk logfile path and parses it out, where possible, into a dict of values
    for all of the fields given in the logfile. If `get_descriptions_only` is true, a dict object
    of only a string with a description of the field (e.g. as a column header) is returned.'''
    #
    fpa, fnm = os.path.split(logfile_path) if logfile_path is not None else ('','')
    re1_ptime  = 'Processing time:\s+(?P<ptime>\d+\.\d+) seconds.'
    re2_input  = 'Input:\s+(?P<input_rct>\d+) reads\s+(?P<input_bct>\d+) bases.'
    re3_ktrim  = 'KTrimmed:\s+(?P<ktrim_rct>\d+) reads [^\s]+\s+(?P<ktrim_bct>\d+) bases'
    re4_otrim  = 'Trimmed by overlap:\s+(?P<otrim_rct>\d+) reads [^\s]+\s+(?P<otrim_bct>\d+) bases'
    re5_totrem = 'Total Removed:\s+(?P<totrem_rct>\d+) reads [^\s]+\s+(?P<totrem_bct>\d+) bases'
    re6_result = 'Result:\s+(?P<res_rct>\d+) reads [^\s]+\s+(?P<res_bct>\d+) bases'
    re7_time   = 'Time:\s+(?P<time_sec>\d+\.\d+) seconds.'
    re8_rdproc = 'Reads Processed:\s+(?P<readsproc>\d+[kmgt])\s+'
    re9_baproc = 'Bases Processed:\s+(?P<basesproc>\d+[kmgt])\s+'
    re_err     = '(?P<err_msg>Exception in thread .*$)\n'
    gkeys = ['filename', 'folder', 'ptime', 'input_rct', 'input_bct', 'ktrim_rct', 'ktrim_bct', 'otrim_rct', 'otrim_bct',
             'totrem_rct', 'totrem_bct', 'res_rct', 'res_bct', 'time_sec', 'readsproc', 'basesproc',
             'err_msg']
    gkey_descr = {'filename': 'file_name', 'folder': 'folder_location_containig_input_file',
                    'ptime': 'processing_time', 'input_rct': 'input_readct', 'input_bct': 'input_basect',
                  'ktrim_rct': 'ktrimmed_readct', 'ktrim_bct': 'ktrimmed_basect',
                  'otrim_rct': 'trimmed_by_overlap_readct', 'otrim_bct': 'trimmed_by_overlap_basect',
                  'totrem_rct': 'total_removed_readct', 'totrem_bct': 'total_removed_basect',
                  'res_rct': 'result_readct', 'res_bct': 'result_basect',
                  'time_sec': 'time_elapsed_sec', 'readsproc': 'reads_processed', 'basesproc': 'bases_processed',
                  'err_msg': 'error_message_full_line'}
    if get_descriptions_only:
        return gkey_descr
    if get_groupkeys_only:
        return gkeys
    assert os.path.isfile(logfile_path), f'File provided is not a valid file path: {logfile_path}'
    with open(logfile_path,'r') as f:
        ftext = f.read()
    #
    vals_dict = {i: '' for i in gkeys}
    vals_dict['filename'] = fnm
    vals_dict['folder'] = fpa
    re_list = [re1_ptime, re2_input, re3_ktrim, re4_otrim, re5_totrem, re6_result, re7_time, re8_rdproc, re9_baproc]
    for my_re in re_list:
        recomp = re.compile(my_re)
        vals = recomp.search(ftext)
        if vals is not None:
            vals_dict.update(vals.groupdict())
    #
    return vals_dict

def parse_bbduk_logfile_folder_to_tsv(folder_path_or_file_list, delim='\t', include_headers=True, verbose=True, quote_char='"'):
    '''
    Goes through a folder containing a big dump of bbduk logfiles and parses each one into the main fields
    using the function above. Compiles all these records into a single string representing a tab-separated-values
    file. The job of actually writing it to a file is left to the calling function.

    Args:
        folder_path_or_file_list (str or list): Should either be a string representing a valid path to a folder or a list of valid
                        paths for a set of files to be parsed in bulk.
        delim (str):                Data field delimiter (default = '\t')
        include_headers (bool):     If false, does not include a header line
        verbose (bool):             If true, reports on progress as it iterates through the files.
        quote_char (bool):          This character is used to enclose the 'error_message' field in particular if
                                    it is provided.

    Returns:
        tsv_string (str):       Python string object that can be directly written to a file as a valid TSV.
    '''
    attach_quote_chars = lambda x: quote_char + str(x) + quote_char if len(str(x))>0 else str(x)
    if isinstance(folder_path_or_file_list,list):
        file_path_list = [i for i in folder_path_or_file_list if os.path.isfile(i)]
        file_input_display_str = str(file_path_list[:5]) + (' ...' if len(file_path_list)>5 else '')
        not_file_path_list = [i for i in folder_path_or_file_list if not os.path.isfile(i)]
        if len(not_file_path_list)>0 and verbose:
            print('At least one of the file paths provided is not valid:')
            print('\n'.join([f'\t{f}' for f in not_file_path_list]))
    else:
        file_input_display_str = folder_path_or_file_list
        assert os.path.isdir(folder_path_or_file_list), f"Folder path provide is not valid: {folder_path_or_file_list}"
        if verbose:
            print('Folder/file-list input is not a \'list\' type, so we are treating it like a folder path. (%s)' % folder_path_or_file_list)
        file_path_list = [os.path.join(folder_path_or_file_list,i) for i in os.listdir(folder_path_or_file_list)]\
    #
    n_files = len(file_path_list)
    header_key_list = parse_bbduk_logfile(get_groupkeys_only=True)
    fmt_string = delim.join(['{%s}' % i for i in header_key_list]) + '\n'
    output_str = ''
    if include_headers:
        header_dict = parse_bbduk_logfile(get_descriptions_only=True)
        output_str += fmt_string.format(**header_dict)
    #
    if verbose:
        print(f'Found {n_files} files in the provided folder/file_list ({file_input_display_str}). Iterating through them...')
    ct=0
    for f in file_path_list:
        bbduk_log_fields = parse_bbduk_logfile(f)
        bbduk_log_fields[header_key_list[-1]] = attach_quote_chars(bbduk_log_fields[header_key_list[-1]])
        output_str += fmt_string.format(**bbduk_log_fields)
        del bbduk_log_fields
        ct += 1
        if verbose:
            print(f'...done with {ct} of {n_files}...', end='\r')
    #
    if verbose:
        print('')
    #
    return output_str

def count_fastq_gz_reads(folder, filename, timing=True):
    '''
    Uses the subprocess module to execute the equivalend of:
        gunzip -c <path> | sed -n '1~4p' | wc -l
    '''
    start_time = phy.mynow()
    fqpath = os.path.join(folder,filename)
    if not os.path.isfile(fqpath):
        print(f'Path {fqpath} does not contain a valid file.')
        return ''
    gunzip_c = subprocess.Popen(['gunzip', '-c', fqpath],
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                text=True)
    # head_10 = subprocess.Popen(['head', '-n', '12'],
    #                            stdin=gunzip_c.stdout,
    #                            stdout=subprocess.PIPE,
    #                            stderr=subprocess.DEVNULL)
    sed_every_fourth = subprocess.Popen(['sed','-n','1~4p'],
                                        stdin=gunzip_c.stdout, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True)
    word_count_lines = subprocess.Popen(['wc','-l'],
                                        stdin=sed_every_fourth.stdout, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True)
    word_count_lines.wait()
    wc_comms = word_count_lines.communicate()
    wc_lines = int(wc_comms[0].strip())
    end_time = phy.mynow()
    if timing:
        print(f'{fqpath} finished in {end_time-start_time}')
    return int(wc_lines), end_time-start_time

def make_post_bbduk_read_count_tsv():
    '''Gets the list of post_bbduk files from the json with all the paths. Goes out to each one and
    counts the reads with the function above. Outputs everything as a tsv in the metadata folder.'''
    post_bbduk_file_locs_json = os.path.join(metadata['server'], 'post_bbduk_file_locations.json')
    post_bbduk_read_cts_tsv = os.path.join(metadata['server'], 'post_bbduk_fastq_read_cts.tsv')
    headers = ['key','folder','file_1','f1_readct','f1_telapse','file_2','f2_readct','f2_telapse']
    post_bbduk_file_locs = phy.json_read(post_bbduk_file_locs_json)
    keylist = list(post_bbduk_file_locs.keys())

    def get_readct_fields(k,as_tuple=False):
        fo, f1nm = post_bbduk_file_locs[k]
        f2nm = f1nm.replace('_1','_2')
        if os.path.isfile(os.path.join(fo,f1nm)):
            f1_readct, f1_te = count_fastq_gz_reads(fo, f1nm, False)
        else:
            f1_readct = -1
            f1_te = 'file_DNE'
        if os.path.isfile(os.path.join(fo,f2nm)):
            f2_readct, f2_te = count_fastq_gz_reads(fo, f2nm, False)
        else:
            f2_readct = -1
            f2_te = 'file_DNE'
        if as_tuple:
            return (k, fo, f1nm, f1_readct, f1_te, f2nm, f2_readct, f2_te)
        else:
            return f'{k}\t{fo}\t{f1nm}\t{f1_readct}\t{f1_te}\t{f2nm}\t{f2_readct}\t{f2_te}'

    with open(post_bbduk_read_cts_tsv,'w') as out_tsv:
        cct = out_tsv.write('\t'.join(headers) + '\n')
        donect = 0; st_time = phy.mynow();
        for key in keylist:
            print(f'running {key} at {phy.mynow()} ({phy.mynow()-st_time}, {donect} of {len(keylist)})', end='\r')
            try:
                ln = get_readct_fields(key)
            except:
                ln = f'{key}\tERROR\t\t\t\t\t\t'
            cct = out_tsv.write(ln + '\n')
            donect += 1

if __name__=='__main__':
    # Meant to run only once, on 3/23/24
    make_post_bbduk_read_count_tsv()




