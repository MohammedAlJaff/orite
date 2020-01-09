

'''
Input: a string representing the genome sequence
Output: a 4 by N matrix. N = genome size/length.
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



'''
Return a cumilative skew array from the output of cumilative_base_count(seq)
See numpy's 'cumsum' function.
'''
def cumilative_skew(base_skew_array):
    return np.cumsum(base_skew_array)


'''
Computes the CG-skew and cumilative CG-skew curves for a given sequences.
Input: Seq string representing a fasta file like the one returned from from seq_from_fasta(file_path)

Output: 
    cg_vals: numpy array with CG-skew values
    cgc_vals: numpy array with cumilative CG-skew values for the given window length 
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


'''
Scales a cruve so that lowest point is at -1 and largest value is +1 by default

Input: curve (numpy array)
Output: Scaled version of curve. Numpy array. 
'''
def scale_skew(curve, default_range=(-1,1)):
    return minmax_scale(X=curve, feature_range=default_range)



'''
Input:
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

