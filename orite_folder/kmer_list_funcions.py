# 4. FUNCTIONS THAT OPERATE ON KMER_LISTS

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
Output: Computes the density for each kmer row in the kmer_info attribute of a region list
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


