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
from sklearn.preprocessing import minmax_scale


'''
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
FUNCTIONS THAT EXTRACT FASTA FILE DATA: START
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''


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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ----- FUNCTIONS THAT EXTRACT FASTA FILE DATA: END -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::





'''
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SLIDING WINDOW SKEW AND CUMILATIVE SKEW CALCULAATION FUNCTIONS.
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''


'''
input:
    Fasta sequence (string)
    window_rad (int). Default = 50 000

Output:

    final_gc: GC Skew values of max rotaed sequence

    final_cgc: cumilative GC-skew value of max rotaded

    max_rotated_fasta: Rotaed fasta sequence based on

    max_cCG_indx_original_offset: Index of maximum cGC skeq value in the original fasta sqeuence.


'''
def max_rotate_seq_and_skew_calc(f, window_radius = 50000):

    initial_pos_index = list(range(len(f)))


    # Inital gc skew calc on original unrotaed sequence and find
    # position at max CG skew value.
    gc, cgc = gc_skew_sliding_window(f, window_rad=window_radius)
    max_indx = np.argmax(gc)
    print('inital max gc skew indx', max_indx)


    # Rotade sequence to begin at max cg-skew value and index.
    new = f[max_indx:] + f[0:max_indx]
    new_pos_index = initial_pos_index[max_indx:] + initial_pos_index[0:max_indx]


    # Calc cg-skew for rotated sequence and find position at max cumilative CG skew value.
    new_gc, new_cgc = gc_skew_sliding_window(new, window_rad=window_radius)
    cGC_max_indx = np.argmax(new_cgc)
    print('max cgc skew indx', cGC_max_indx)


    # Final rotation of sequence and final gc calc results
    final = new[cGC_max_indx:] + new[0:cGC_max_indx]
    final_pos_index = new_pos_index[cGC_max_indx:]+new_pos_index[0:cGC_max_indx]


    final_gc, final_cgc = gc_skew_sliding_window(final, window_rad=window_radius)


    # Return sequence rotated to max
    max_rotated_fasta = final

    max_cCG_indx_original_offset = (max_indx + cGC_max_indx)% len(f)

    return final_gc, final_cgc, max_rotated_fasta, max_cCG_indx_original_offset

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

    #print('seq len: ', seq_length)
    #print('window_rad: ', window_rad)
    #print('extended seq length: ', extetend_seq_length)


    gc_vals = np.zeros([seq_length])

    for i in range(window_rad, window_rad+seq_length):

        left_indx = (i-window_rad)
        right_indx = (i+window_rad)



        n_c_in_window = (c_b_c[1,right_indx] - c_b_c[1, left_indx])
        n_g_in_window = (c_b_c[2,right_indx] - c_b_c[2, left_indx])

        if (n_g_in_window+n_c_in_window == 0):
            print('opps at: ', i)

        gc_skew = (n_g_in_window-n_c_in_window)/(n_g_in_window+n_c_in_window)
        gc_vals[i-window_rad] = gc_skew


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
UNUSED AS OF 5TH DEC

Function: Returns the sum of the squared distances for each element in the arrays. Arrays must be of same length
x: numpy array
y: numpy array
'''
def sum_of_squared_distances(x, y):
    return np.sum((x-y)*(x-y))


'''
UNUSED AS OF 5TH DEC

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
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
NON CODING REGIONS CLASS: START
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''


'''
Class: Non coding region object.
'''
class NC_region:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        self.length = stop - start +1
        self.kmer_info = {}
        self.max_relative_start = -1
        self.max_relative_stop = -1

    #def __str__(self):
        #return str(self.start)

    def add_sequence(self,string):
        self.sequence = string

    def add_cgc_val(self,cgc_val):
        self.cgc_val = cgc_val

    def add_gc_val(self,gc_val):
        self.gc_val = gc_val

    def add_kmer_counts(self, k):
        kmer_dict = get_kmers(self.sequence, k, circular=False)
        self.kmer_info.update({k: kmer_dict})

    def filter_kmer_by_occurence(self, n):
        for key, value in self.kmer_info.items():
            self.kmer_info[key] = get_kmer_by_occurence(value, n)

    def filter_top_n_kmers(self, n):
        for key, value in self.kmer_info.items():
            self.kmer_info[key] = get_top_n_kmers(value, n)
    def calc_kmer_density(self):
        for key, value in self.kmer_info.items():
            self.kmer_info[key] = calc_kmer_density(value)

    def add_max_relative_pos(self, max_offset, whole_genome_length):
        self.max_relative_start = (self.start - max_offset)%whole_genome_length
        self.max_relative_stop = (self.stop - max_offset)%whole_genome_length







'''
Function that assigns a cgc VALUE BASED SCORE TO A REGION
'''

def calc_score_for_NC_region(cgc_curve, nc_obj):
    inter = list(range(nc_obj.start, nc_obj.stop+1))
    region_cgc_vals = cgc_curve[inter]

    score= -np.average(region_cgc_vals)

    nc_obj.add_cgc_val(score)



'''
EXTRACTS A NON CODING SEQUENCE AND ASSIGNS IT TO THE NON CODING REGION BASED ON
ITS START AND STOP PLACES.
'''
def add_sequence_to_region(seq, nc_obj):
    nc_obj.add_sequence(seq[nc_obj.start:nc_obj.stop +1])








'''
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
FUNCTIONS THAT OPERATE ON KMER_LISTS
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''


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
kmer_list: A list of tuples
[('ATTTAT', [#OCCURANCE, POS1, POS2, POS,3]),
('GGGTTTT', [#OCCURANCE, POS1, POS2, POS,3]),
('GGGGCCG', [#OCCURANCE, POS1, POS2, POS,3]),
...
...
    ]

For each tuple A tuple. First element represents the kmer (String), the 2nd elelment is a list. The first element in the list
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





'''
Input: A kmer_list
Output:
'''

def calc_kmer_density(kmer_list):
    kmer_list_plus_d = []
    for touple in kmer_list:
        this_kmer = touple[0]
        this_list = touple[1]
        this_occurance = this_list[0]
        start = this_list[1]
        stop = this_list[len(this_list)-1]
        this_density = this_list[0] /((stop-start))
        kmer_list_plus_d.append((this_kmer, deepcopy(this_list), this_density))
    return kmer_list_plus_d





'''
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
FUNCTIONS THAT OPERATE ON A REGION_LISTS
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''

'''
Filters a regions_list by length
'''

def filter_regions_by_length(regions_list, length):
    new_region_list = []

    for region in regions_list:
        if region.length > length:
            new_region_list.append(region)
    return new_region_list



'''
Sorts a regions_list by each region objects score attribute.
NOTE:
THIS FUNCTION OPERATES ON THE INPUTED region_list itself.
It doesnt return anything.
'''
def sort_regions_by_score(regions_list):
    regions_list.sort(key=lambda x: x.cgc_val, reverse = True)
    #regions_list.sort(reversed = True)










'''
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
FUNCTIONS THAT ARE USED TO EXTRACT AND PARSE GENBANK DATA: START
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
'''

'''
Input: file path to a genbank file
Output: a tuple where the first element contains a list of tuples of the non-coding regions on both strands,
the second being a list of the same regions but with each base position as an element in the list.

The function utilizes various other functions contained in oritelib.
'''

def genbank_to_non_coding_intervals(file_path):
    length = get_length_genbank(file_path)
    raw_p = raw_positions_strings_from_genbank(file_path)

    raw_plus, raw_neg = split_strands(raw_p)


    valid_plus = make_into_valid_pos(raw_plus)
    valid_neg = make_into_valid_pos(raw_neg)


    non_coding_plus = get_non_coding_intervals(valid_plus, length)
    non_coding_neg = get_non_coding_intervals(valid_neg, length)



    true_nc_positions = get_true_nc_positions(non_coding_plus, non_coding_neg)
    true_nc_intervals = position_list_to_intervals(true_nc_positions)





    return true_nc_intervals, true_nc_positions, non_coding_plus, non_coding_neg


'''
Input two lists of non-coding regions.
Output set of all true non-coding positions

HELPER FUNCTION 1: Input a list of intervals, output a set of all positions in intervals

'''
def interval_list_to_position_set(interval_list):
    pos_set = set()

    for interval in interval_list:
        start = interval[0]
        stop = interval[1]

        interval_range = list(range(start,stop+1))
        pos_set.update(interval_range)

    return pos_set

#-------------------------

'''
Input: non-coding intervals for each strand
Output: common non-coding regions for both strands
'''
def get_true_nc_positions(nc_plus_intervals, nc_neg_intervals):
    nc_plus_set = interval_list_to_position_set(nc_plus_intervals)
    nc_neg_set = interval_list_to_position_set(nc_neg_intervals)

    intersection_set = nc_plus_set.intersection(nc_neg_set)
    intersect_set_list = list(intersection_set)
    intersect_set_list.sort()
    return intersect_set_list

#-------------------------



'''

Input: list of all non-coding positions.
Output: List of tuples with non-coding intervals'''

def position_list_to_intervals(pos_list):

    crap_bag = []

    current_start = pos_list[0]
    current_stop = -1

    for i in  range(1, len(pos_list)-1):

        if pos_list[i-1] +1 == pos_list[i] and  pos_list[i]+1!=pos_list[i+1]:
            current_stop = pos_list[i]

            crap_bag.append((current_start, current_stop))

            #current_start = pos_list[i+1]

        elif pos_list[i-1] +1 != pos_list[i] and  pos_list[i]+1==pos_list[i+1]:
            current_start = pos_list[i]

        #elif pos_list[i-1] +1 != pos_list[i] and  pos_list[i]+1!=pos_list[i+1]:



    if pos_list[i+1] == current_start:
        current_stop = pos_list[i+1]
        crap_bag.append((current_start, current_stop))

    if pos_list[i]+1==pos_list[i+1]:
        current_stop = pos_list[i+1]
        crap_bag.append((current_start,current_stop))

    return crap_bag

#-------------------------


#-------------------------
'''
Input: List of tuples with non-coding intervals
Output: list of range of each non-coding interval
'''
def interval_list_to_range_list(interval_list):

    crap_list = []
    for touple in interval_list:
        x = list(range(touple[0], touple[1]+1))
        crap_list.append(x)

    return crap_list

#-------------------------


'''
INTERNAL / PRIVATE / HELPER
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

#-------------------------


'''
INTERNAL / PRIVATE / HELPER
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

#-------------------------


'''
INTERNAL / PRIVATE / HELPER
'''

def make_into_valid_pos(strand_pos_raw):

    interval_list = []

    for row in strand_pos_raw:
        m = re.findall(r'\d+', row)


        start_p = int(m[0])
        stop_p = int(m[len(m)-1])
        interval_list.append((start_p, stop_p))

    return interval_list

#-------------------------

'''
INTERNAL / PRIVATE / HELPER
'''
def get_non_coding_intervals(coding_intervals, length):
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

    if coding_intervals[i][1] != length:
        non_coding_start = coding_intervals[i][1]+1
        non_coding_stop = length-1
        non_coding_intervals.append((non_coding_start, non_coding_stop))
    return non_coding_intervals

#-------------------------


'''
input: genbank path (string)
output: genome length/Size (int)

Returns the length of the organism genome 'contained' in the genbank file
'''
def get_length_genbank(file_path):
    recs = [rec for rec in SeqIO.parse(file_path, "genbank")]
    for rec in recs:
        length = len(rec.seq)
    return length

#-------------------------



'''
PUBLIC
'''
def extract_seq_from_non_coding_intervals(non_coding_intervals, seq):

    non_coding_seqs = []

    for interval in non_coding_intervals:

        non_coding_seqs.append(seq[interval[0]:interval[1]])

    return non_coding_seqs

'''

'''

'''
Input: a list of non-coding intervals
Output: a list of the objects form NC_region class, containing sequence and cgc values
'''

def nc_intervals_to_nc_objects(nc_intervals, seq, cgc):
    nc_objcts = []
    for i in range(len(nc_intervals)):
        nc_obj = orite.NC_region(nc_intervals[i][0], nc_intervals[i][1])
        orite.add_sequence_to_region(seq, nc_obj)
        orite.calc_score_for_NC_region(cgc, nc_obj)
        nc_objcts.append(nc_obj)
    return nc_objcts


def get_kmers_from_region_list(region_list, k_array):
    new_list = []
    for region in region_list:
        for k in k_array:
            region.add_kmer_counts(k)
        new_list.append(region)
    return new_list


def filter_region_list_by_kmer_occurence(region_list, n):
    new_list = []
    for region in region_list:
        region.filter_kmer_by_occurence(n)
        new_list.append(region)
    return new_list


def has_empty_kmer_info(nc_region):

    for key, value in nc_region.kmer_info.items():
        if len(nc_region.kmer_info[key]) != 0:
            return False

    return True


def filter_empty_kmer_regions(region_list):
    new_list = []

    for region in region_list:
        if not has_empty_kmer_info(region):
            new_list.append(region)
    return new_list



def plot_region_list(region_list, curve):

    region_intervals = []

    for region in region_list:
        this_interval = (region.start, region.stop)
        region_intervals.append(this_interval)

    regions_pos_set = orite.interval_list_to_position_set(region_intervals)
    region_pos_list = list(regions_pos_set)
    region_pos_list.sort()

    relevant_pos_list = np.array(region_pos_list)

    relevant_curve_point = curve[relevant_pos_list]


    plt.figure(figsize=[40,10])
    plt.plot(curve)
    plt.plot(relevant_pos_list, relevant_curve_point, 'x')
    plt.title(str(len(region_list)))

def print_region_list_kmer_info(region_list):
    i = 0
    for region in region_list:
        print('region:', i, '---',' score: ', region.cgc_val, '---- pos: ', region.start, '---- max_relative_start_pos', region.max_relative_start)

        for key, value in region.kmer_info.items():
            print('\tk=', key,)
            for thing in value:
                print('\t', thing[0], ' - ', thing[1])

        print('\t-------')

        i = i+1

def add_max_relative_position(region_list, genome_length, max_offset):
    new_list = []
    for region in region_list:
        add_max_relative_position(region, max_offset, genome_length)
        new_list.append(region)

    return new_list




'''
Z-curve calculation of a genome sequence sequence

input: A sequence string of a genome

output: The xn, yn, and zn components of the Z curve of the genome.
'''
def calc_z_curve(seq):

    if type(seq[0]) != str :
        raise Exception('Must input a valid strng sequence')

    else:
        # Acumilative base counts across genome
        cum_base_count = cumilative_base_count(seq)

        an = cum_base_count[0,:]
        cn = cum_base_count[1,:]
        gn = cum_base_count[2,:]
        tn = cum_base_count[3,:]

        # z-curve component calculations
        xn = (an+gn) - (cn+tn) #RY
        yn = (an+cn) - (gn+tn) # MK
        zn = (an+tn) - (cn+gn)


    return xn, yn,





'''
Scales a cruve so that lowest point is at -1 and largest value is +1 by default

input: curve (numpy array)
output:

'''

def scale_skew(curve, default_range=(-1,1)):
    return minmax_scale(X=curve, feature_range=default_range)
