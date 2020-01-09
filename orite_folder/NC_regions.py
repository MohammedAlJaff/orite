# 3. NON CODING REGIONS OBJECTS AND THINGS THAT OPERATE ON THESE
'''
Class: Non coding region object.
'''

class NC_region:

    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        self.length = stop - start +1   # Sequence length.
        self.kmer_info = {} # Kmer content info
        self.max_relative_start = -1 # Position in max-rotated fasta file
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
        if self.start >= max_offset:
            self.max_relative_start = self.start - max_offset
            self.max_relative_stop = self.stop - max_offset
        else:
            self.max_relative_start = whole_genome_length - (max_offset - self.start)
            self.max_relative_stop = whole_genome_length - (max_offset - self.stop)


    def filter_out_empty_kmer_lists_in_kmer_dict(self):

        new_dict = {}

        for key, value in self.kmer_info.items():
            if len(value) != 0:
                new_dict.update({key:value})

        self.kmer_info = new_dict


    def sort_kmer_info_by_density(self):

        new_dict = {}
        for key, value in self.kmer_info.items():

            value.sort(key= lambda x:x[2], reverse=True)

            new_dict.update({key:value})

        self.kmer_info = new_dict

    def remove_kmer_overlap(self):
        new_dict = {}
        for key, value in self.kmer_info.items():
            new_value = []
            for listy in value:
                overlap = False
                for i in range(1,len(listy[1])-1):
                    if (listy[1][i+1] - listy[1][i]) < key:
                        overlap = True
                        break

                if not overlap:
                    new_tuple = (listy[0],listy[1],listy[2])
                    new_value.append(new_tuple)
            new_dict.update({key:new_value})
        self.kmer_info = new_dict




'''
Function that assigns a cgc VALUE BASED SCORE TO A REGION

No return value. 
'''
def calc_score_for_NC_region(cgc_curve, nc_obj, rotated = False):

    if rotated:
        inter = list(range(nc_obj.max_relative_start, nc_obj.max_relative_stop+1))
        region_cgc_vals = cgc_curve[inter]

        score= -np.average(region_cgc_vals)

        nc_obj.add_cgc_val(score)
    else:
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
This operates on a single nc-region object. 
returns true of the kmer 
'''
def has_empty_kmer_info(nc_region):

    for key, value in nc_region.kmer_info.items():
        if len(nc_region.kmer_info[key]) != 0:
            return False

    return True