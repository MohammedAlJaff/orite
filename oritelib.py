import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, misc, ndimage
from copy import deepcopy
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re



'''
input: file path to a fasta file contaning a single sequence for the complete genome of a bacteria

output: A string with only the base sequence
'''
def seq_from_fasta(file_path):
    f=open(file_path, "r")

    all_lines = f.readlines()

    sequence = ''

    for i in range(1, len(all_lines)):
        sequence =  sequence + all_lines[i][0:len(all_lines[i])-1]

    return sequence

#-------------------------

'''
input: a string representing the genome sequence
output: a 4 by N matrix. N = genome size/length.
    each row reprensts the cumilative count of a base nucleotide at a given position.
    first row: A, second row: C, third row: G, fourth row: T
'''

def cumilative_base_count(seq):
    cum_base_array = np.zeros([4, len(seq)])

    cum_a_i = 0
    cum_c_i = 0
    cum_g_i = 0
    cum_t_i = 0

    sec_length = len(seq)

    for indx in range(len(seq)):
        base = seq[indx]
        if base == 'A':
            cum_a_i = cum_a_i + 1

        if base == 'C':
            cum_c_i = cum_c_i + 1

        if base == 'G':
            cum_g_i = cum_g_i + 1

        if base == 'T':
            cum_t_i = cum_t_i + 1

        cum_base_array[0,indx] = cum_a_i
        cum_base_array[1,indx] = cum_c_i
        cum_base_array[2,indx] = cum_g_i
        cum_base_array[3,indx] = cum_t_i

    return cum_base_array

#-------------------------

'''

'''

def gc_skew_sliding_window(seq, window_rad=500):

    seq_length = len(seq)

    if window_rad > seq_length:
        raise Exception('ERROR - Window size is larger than sequence length')


    extetend_seq =  seq[seq_length-window_rad:] + seq + seq[0:window_rad]
    extetend_seq_length = len(extetend_seq)

    c_b_c = cumilative_base_count(extetend_seq)

    print('seq len: ', seq_length)
    print('window_rad: ', window_rad)
    print('extended seq length: ', extetend_seq_length)


    gc_vals = np.zeros([seq_length])

    for i in range(window_rad, window_rad+seq_length):

        #print(i)
        left_indx = (i-window_rad)
        right_indx = (i+window_rad)



        n_c_in_window = (c_b_c[1,right_indx] - c_b_c[1, left_indx])
        n_g_in_window = (c_b_c[2,right_indx] - c_b_c[2, left_indx])

        if (n_g_in_window+n_c_in_window == 0):
            print('opps at: ', i)

        gc_skew = (n_g_in_window-n_c_in_window)/(n_g_in_window+n_c_in_window)
        gc_vals[i-window_rad] = gc_skew

        '''
        if ((i-window_rad) % 100000 )== 0:
            print('left_indx: ', left_indx)
            print('current base indx: ', i)
            print('right_indx: ', right_indx)
        '''

    c_gc_vals = cumilative_skew(gc_vals)

    return gc_vals, c_gc_vals

#-------------------------

'''
'''

def at_skew_sliding_window(seq, window_rad=500):

    seq_length = len(seq)

    if window_rad > seq_length:
        raise Exception('ERROR - Window size is larger than sequence length')


    extetend_seq =  seq[seq_length-window_rad:] + seq + seq[0:window_rad]
    extetend_seq_length = len(extetend_seq)

    c_b_c = cumilative_base_count(extetend_seq)

    print('seq len: ', seq_length)
    print('window_rad: ', window_rad)
    print('extended seq length: ', extetend_seq_length)


    gc_vals = np.zeros([seq_length])

    for i in range(window_rad, window_rad+seq_length):

        #print(i)
        left_indx = (i-window_rad)
        right_indx = (i+window_rad)

        n_c_in_window = (c_b_c[0,right_indx] - c_b_c[0, left_indx])
        n_g_in_window = (c_b_c[3,right_indx] - c_b_c[3, left_indx])

        if (n_g_in_window+n_c_in_window == 0):
            print('opps at: ', i)

        gc_skew = (n_g_in_window-n_c_in_window)/(n_g_in_window+n_c_in_window)
        gc_vals[i-window_rad] = gc_skew

        '''
        if ((i-window_rad) % 100000 )== 0:
            print('left_indx: ', left_indx)
            print('current base indx: ', i)
            print('right_indx: ', right_indx)
        '''

    c_gc_vals = cumilative_skew(gc_vals)

    return gc_vals, c_gc_vals




#-------------------------

'''

return a cumilative skew array
'''
def cumilative_skew(base_skew_array):
    return np.cumsum(base_skew_array)

# ----------------


'''
return: a smoothed version of input array

Note:
smoothing param is dependent on resolution.
ndi.gaussian_filter1d(cGCsk, sigma=110)'''

def smooth_curve(raw_curve, smoothing_param=60):

    smoothed_curve = ndimage.gaussian_filter1d(raw_curve, sigma=smoothing_param)

    return smoothed_curve

'''
Function: Returns the sum of the squared distances for each element in the arrays. Arrays must be of same length
x: numpy array
y: numpy array
'''
def sum_of_squared_distances(x, y):
    return np.sum((x-y)*(x-y))


'''
Returns the value of the linear function with paramters k (slope), m (intersect) over the specified position.
k: double, float.
m: double, float.
pos: Numpy array.
'''
def line_array(k,m,pos):
    y = k*pos + m
    return y

'''
sequence: string of genome sequence
kmer_length: int on how long kmers to look for
circular: Boolean, True by default. Use False if genome is not circular or if a subsequence is used

Output: A list of tuples. First element in tuple represents the kmer (String), the other elelment is a list. The first element in the list
represents how many occurences there are of the kmer in the investigated sequence. The following indices store the
positions on which the first base in the kmer is located on the sequence. The tuple are based in descending order
by how many occurences there are of the kmer in the sequence.
'''
def get_kmers(sequence, kmer_length, circular = True):
    kmer = ""
    kmer_dict = {}
    sequence_length = len(sequence)

    if circular:
        sequence = sequence + sequence[0:kmer_length-1]

        for i in range(sequence_length):
            kmer = sequence[i:(i+kmer_length)]

            if kmer in kmer_dict:
                value_list = kmer_dict[kmer]
                value_list[0] = value_list[0]+1
                value_list.append(i)
                kmer_dict[kmer] = value_list

            else:
                kmer_dict.update({kmer: [1,i]})
    else:
        for i in range(sequence_length-(kmer_length-1)):
            kmer = sequence[i:(i+kmer_length)]

            if kmer in kmer_dict:
                value_list = kmer_dict[kmer]
                value_list[0] = value_list[0]+1
                value_list.append(i)
                kmer_dict[kmer] = value_list

            else:
                kmer_dict.update({kmer: [1,i]})

    kmer_dict = sorted(kmer_dict.items(), key =
             lambda kv:(kv[1], kv[0]), reverse = True)

    return kmer_dict

'''
kmer_list: A list of tuples. First element in tuple represents the kmer (String), the other elelment is a list. The first element in the list
represents how many occurences there are of the kmer in the investigated sequence. The following indices store the
positions on which the first base in the kmer is located on the sequence. The tuple are based in descending order
by how many occurences there are of the kmer in the sequence.

n: Number of tuples to return

Output: A list of the n tuples with kmers that occured the most in a sequence
'''
def get_top_n_kmers(kmer_list, n):
    return deepcopy(kmer_list[0:n])


'''
kmer_list: A tuple. First element represents the kmer (String), the other elelment is a list. The first element in the list
represents how many occurences there are of the kmer in the investigated sequence. The following indices store the
positions on which the first base in the kmer is located on the sequence. The tuple are based in descending order
by how many occurences there are of the kmer in the sequence.

n: A threshold for how many times a kmer must occur in a sequence for it to be considered relevant

Output: A list of the tuples with kmers that occured more or equal to n number of times in the sequence
'''

def get_kmer_by_occurence(kmer_list, n):
    new_kmer_list = []
    for i in range(len(kmer_list)):
        if (kmer_list[i][1][0] >= n):
            new_tuple = deepcopy(kmer_list[i])
            new_kmer_list.append(new_tuple)
    return new_kmer_list




#-------------- Genbank functions

'''

'''

def genbank_to_non_coding_intervals(file_path):

    raw_p = raw_positions_strings_from_genbank(file_path)

    raw_plus, raw_neg = split_strands(raw_p)


    valid_plus = make_into_valid_pos(raw_plus)
    valid_neg = make_into_valid_pos(raw_neg)


    non_coding_plus = get_non_coding_intervals(valid_plus)
    non_coding_neg = get_non_coding_intervals(valid_neg)

    return non_coding_plus, non_coding_neg


'''
PRIVATE / HELPER
'''

def raw_positions_strings_from_genbank(file_path):

    recs = [rec for rec in SeqIO.parse(file_path, "genbank")]


    raw_positions_str = []

    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            x = str((feat.location))
            raw_positions_str.append(x)


    return raw_positions_str

'''
PRIVATE / HELPER
'''
def split_strands(raw_pos_list):
    plus_pos = []
    neg_pos  = []

    for row in raw_pos_list:


        if row.find('+') == -1:
            neg_pos.append(row)

        if row.find('-') == -1:
            plus_pos.append(row)

    return plus_pos, neg_pos



'''
PRIVATE / HELPER
'''

def make_into_valid_pos(strand_pos_raw):

    interval_list = []

    for row in strand_pos_raw:
        m = re.findall(r'\d+', row)

        #print(len(m))

        start_p = int(m[0])
        stop_p = int(m[len(m)-1])
        interval_list.append((start_p, stop_p))

    return interval_list


'''
PRIVATE / HELPER
'''
def get_non_coding_intervals(coding_intervals):
    non_coding_intervals = []



    if coding_intervals[0][0] != 0:

        ith_non_coding_start = 0
        ith_non_coding_stop = coding_intervals[0][0]-1

        non_coding_intervals.append((ith_non_coding_start, ith_non_coding_stop))


    for i in range(1, len(coding_intervals)):
        ith_non_coding_start = coding_intervals[i-1][1]+1
        ith_non_coding_stop = coding_intervals[i][0]-1

        if ith_non_coding_stop-ith_non_coding_start>1:
            non_coding_intervals.append((ith_non_coding_start, ith_non_coding_stop))

    return non_coding_intervals





'''
PUBLIC
'''
def extract_seq_from_non_coding_intervals(non_coding_intervals, seq):

    non_coding_seqs = []

    for interval in non_coding_intervals:

        non_coding_seqs.append(seq[interval[0]:interval[1]])

    return non_coding_seqs
