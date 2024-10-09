import os.path

import argparse, gzip
import numpy as np

import settings_configs
import phylogeny_utilities.utilities as phy


parser = argparse.ArgumentParser()

# Helper functions up here:
#
# converts a (string) row of the depth file to a tuple of type 
#   ('str','int','int')
frow_to_tup = lambda x: (x[0],int(x[1]),int(x[2]))
nparray_is_consecutive = lambda x: np.all(x==np.arange(x[0],x[0]+x.shape[0]))

def build_contig_coverage_dict(depth_npstruct: np.ndarray, f2_min_row_index: np.ndarray,
                               f2_max_row_index: np.ndarray, coverage_only: bool = True) -> tuple:
    '''
    Intermediate function between the samtools depth file and the numpy dict output.

    Takes a structured numpy array representing the samtools depth file, plus a set of
    indices representing the min and max rows for each contig in that array (which
    requires that it has been sorted already). Converts that to a 'dict' object
    where keys are the contig names and values are np.int32 arrays with a coverage
    value at each position.

    The way this function operates requires making some assumptions about the depth file, and
    by extension the np array that is input. For one, we assume that it has been sorted so
    that anything between the min_row_index and max_row_index for a particular contig *ALL*
    belongs to the same contig. For another, we assume that *ALL* positions in a contig are
    listed in the depth file because we have forced that option at the command line. Still,
    we do a check for both of these conditions along the way and the truth values are stored
    in the Assertion Array, the third return value.

    Args:
        depth_npstruct (np.ndarray): samtools depth file as a 3-field numpy.struct array (see
                                below).
        f2_min_row_index (np.ndarray): numpy array with one entry for every contig. represents
                                the 0-indexed starting row for that contig in the depth file.
        f2_max_row_index (np.ndarray): same as previous, but for the *last* row.
        coverage_only (bool): if True, only returns the coverage dict with None for the other two.

    Returns:
        f_coverage (dict):  dict of the form {<contig_name>: <coverage_array>} where <coverage_array>
                            is a numpy int32 arrayhas the same shape as the length of the original contig.
        f_positions (dict): dict of the form {<contig_name>: <position_array>} where <position_array>
                            contains the 'position' field from the depth file. If we've run samtools
                            correctly this should just be a consective array of integers, and so far it
                            has been.
        assertion_array (dict): For each contig in the assembly, a 1 if more than 1 'key' value were
                            between the min and max contig positions (i.e. the array was not sorted on
                            'key'). Also for each contig, a 1 if the 'position' values do not form a
                            consecutive integer array (also could be caused by non-sortedness, but also
                            by incomplete depth reporting).
    '''
    # Initialize
    update_user_every=100
    f_coverage = {}; f_positions={};
    n_seqs = f2_min_row_index.shape[0]
    assertion_array = np.zeros((n_seqs,2),dtype=np.uint8)

    for i in range(n_seqs):
        # for each contig, get the subarray between its bounds:
        subarray=depth_npstruct[f2_min_row_index[i]:(f2_max_row_index[i] + 1)]
        #Make sure sorting was done correctly:
        if np.unique(subarray['key']).shape[0]!=1:
            assertion_array[i,0]=1
        # Make sure no columns were left out:
        if nparray_is_consecutive(subarray['col']):
            assertion_array[i,1]=1
        my_key = str(subarray['key'][0])
        # if the assertions are all good, we just pop that struct array between
        #   the boudns into a dict value and report that out.
        f_coverage[my_key] = subarray['depth']
        f_positions[my_key] = subarray['col']
        if i % update_user_every ==0:
            print(f' i={i} out of {n_seqs}', end='\r')

    if coverage_only:
        return f_coverage, None, None
    else:
        return f_coverage, f_positions, assertion_array


def samtools_depth_to_numpy_dict(depth_file_path: str, out_numpy_path=None, return_dicts=False):
    '''
    Reads the samtools depth file and converts it to a numpy '.npz' file which is essentially a 'dict'
    where the values are np.arrays. Here the 'keys' in the output are the names of the contigs in
    the contig file and the values are np.int32 arrays containing the coverage at each position
    along the original contig.

    Args:
        depth_file_path (str):  path to the samtools '.depth' file.
        out_numpy_path (str):   path to the output file, fed straight to np.savez(...)
        return_dicts (bool):    if True, returns the dictionaries that were created in the process. If
                                False, returns (None, None, None).

    Returns:
        depth_nps_sorted (np.ndarray):  structured numpy array of the depth file, sorted by 'key' and 'col'
        ctg_min_row (np.ndarray):       0-indexed row numbers where each contig starts in the sorted array
        ctg_max_row (np.ndarray):       0-indexed row numbers where each contig ends in the sorted array

    '''
    # **(1)** Read the depth file as a list of -tuples where the first column is a string
    #           and the second two columns are converted to int.
    if os.path.splitext(depth_file_path)[1]=='.gz':
        depth_file_lns = phy.gz_read_to_text_lines(depth_file_path)
        depth_tuples = [frow_to_tup(i.split('\t')) for i in depth_file_lns]
    else:
        with open(depth_file_path,'r') as depth_file:
            depth_tuples = [frow_to_tup(i.strip().split('\t')) for i in depth_file.readlines()]

    # Get the width of the key-column to make an np.struct:
    kw = max([len(f[0]) for f in depth_tuples])
    # **(2)** Turn the list of tuples into a np.structured array.
    depth_npstruct=np.array(depth_tuples,dtype=[('key',f'<U{kw}'),('col','i4'),('depth','i4')])
    #     **  Sort the array first by 'key' (i.e. seq name) and then column number:
    depth_nps_sorted=np.sort(depth_npstruct, order = ['key','col'])

    # **(3)** Get the relevant start/stop positions of every contig in the tall depth file:
    #     ...max rows for each contig = points where the key column changes:
    ctg_max_row=np.where(depth_nps_sorted['key'][:-1]!=depth_nps_sorted['key'][1:])[0]
    #     ...min rows for each contig = row 0, plus the rows right after the max-rows:
    ctg_min_row=np.hstack(([0,],ctg_max_row[0:-1]+1))

    # **(4)** Call the specialized function for constructing the output once we have the key inputs:
    depth_by_seq, depth_positions, asserts = build_contig_coverage_dict(depth_nps_sorted, ctg_min_row, ctg_max_row)
    if out_numpy_path is not None:
        np.savez(out_numpy_path,**depth_by_seq)
    # np.savez(out_extradata_path, assertion_matrix=asserts)

    if return_dicts:
        return depth_nps_sorted, ctg_min_row, ctg_max_row

def get_assembly_and_depth_filepaths(bioproject: str, runid: str) -> tuple:
    '''
    Given a bioproject and runid, returns the paths to the depth file and the assembly file.

    Args:
        bioproject (str):   bioproject name
        runid (str):        runid name

    Returns:
        depth_file_path (str):  path to the depth file
        assembly_file_path (str):   path to the assembly file
    '''
    rm_out = settings_configs.server_paths['output_folders_homeserver']['read_mapping_out']
    rm_trans = os.path.abspath(os.path.join(rm_out, '..', 'read_mapping_transfer'))
    megahit_out = settings_configs.server_paths['output_folders_homeserver']['megahit_out']
    depth_file_path1 = os.path.join(rm_out, f'{bioproject}_{runid}', f'{bioproject}_{runid}.depth.gz')
    depth_file_path2 = os.path.join(rm_trans, f'{bioproject}_{runid}', f'{bioproject}_{runid}.depth.gz')
    depth_file_path3 = os.path.join(rm_trans, f'{bioproject}_{runid}', f'{bioproject}_{runid}.depth')
    if os.path.isfile(depth_file_path1):
        depth_file_path = depth_file_path1
    elif os.path.isfile(depth_file_path2):
        depth_file_path = depth_file_path2
    elif os.path.isfile(depth_file_path3):
        depth_file_path = depth_file_path3
    else:
        raise FileNotFoundError(f'Could not find the depth file for {bioproject}_{runid}.')
    assembly_file_path = os.path.join(megahit_out, f'{bioproject}_{runid}', 'final.contigs.fa')
    return assembly_file_path, depth_file_path


    
def parse_args():
    global parser
    parser.add_argument('-i','--input', dest='input', type=str, help='Path to the .depth file that should be converted to a '
                                                                            'coverage vector along the length of each contig.')
    parser.add_argument('-o','--output', dest='output', type=str, help='Path to the .npz file that the result should be written '
                                                            'to. If this file doesn ot end with \'.npz\' that will be added.')
    parser.add_argument('-a','--auxoutput', dest='auxoutput', type=str, 
                        help='Path to another .npz file that can have the auxiliary data written to it, namely the asserts '
                            ' matrix if it\'s needed for debugging.')
    
    args = parser.parse_args()
    return args

if __name__=='__main__':
    cmd_args = parse_args()
    in_depth_f = cmd_args.input
    out_depth_npz = cmd_args.output
    aux_out_npz = cmd_args.auxoutput
    samtools_depth_to_numpy_dict(in_depth_f, out_depth_npz)

