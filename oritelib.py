import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, misc, ndimage


'''
input: file path to a fasta file contaning a single sequence for the complete genome of a bacteria

output: A string with only the base sequence
'''
def seq_from_fasta(file_path):
    f=open(file_path, "r")

    all_lines = f.readlines()

    sequence = ''

    for i in range(1, len(all_lines)):
        sequence =  sequence + all_lines[i][0:len(all_lines[i])-2]

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
    gc_vals = np.zeros([seq_length])

    c_b_c = cumilative_base_count(seq)

    for i in range(window_rad, seq_length-window_rad):

        left_indx = (i-window_rad) % seq_length
        right_indx = (i+window_rad) % seq_length

        n_c_in_window = (c_b_c[1,right_indx] - c_b_c[1, left_indx])
        n_g_in_window = (c_b_c[2,right_indx] - c_b_c[2, left_indx])

        gc_skew = (n_g_in_window-n_c_in_window)/(n_g_in_window+n_c_in_window)
        gc_vals[i] = gc_skew


    return gc_vals

#-------------------------

'''
'''

def at_skew_sliding_window(seq, window_rad=500):

    seq_length = len(seq)
    at_vals = np.zeros([seq_length])

    c_b_c = cumilative_base_count(seq)

    for i in range(window_rad, seq_length-window_rad):

        left_indx = (i-window_rad) % seq_length
        right_indx = (i+window_rad) % seq_length

        n_a_in_window = (c_b_c[0,right_indx] - c_b_c[0, left_indx])
        n_t_in_window = (c_b_c[3,right_indx] - c_b_c[3, left_indx])

        at_skew = (n_a_in_window-n_t_in_window)/(n_a_in_window+n_t_in_window)
        at_vals[i] = at_skew


    return at_vals


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
