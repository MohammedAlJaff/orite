
# FUNCTIONS TO SORT AND FILTER LISTS OF NC_REGION OBJECTS

'''
Input: a list of non-coding intervals from genbank_to_non_coding_intervals()
Output: a list of the objects form NC_region class. 
'''
def nc_intervals_to_nc_objects(nc_intervals, og_seq):
    nc_objcts = []
    for i in range(len(nc_intervals)):
        nc_obj = NC_region(nc_intervals[i][0], nc_intervals[i][1])
        add_sequence_to_region(og_seq, nc_obj)
        #orite.calc_score_for_NC_region(cgc, nc_obj)
        nc_objcts.append(nc_obj)
    return nc_objcts


'''
Sortednc kmers in a region list.
NOTE: The objects themselves do not get resorted, 
Only the rows of each kmers
'''
def sort_region_list_on_density(region_list):

    new_list = []

    for region in region_list:
        region.sort_kmer_info_by_density()
        new_list.append(region)
    return new_list



'''
Input:
Outout: 
Inputs: A list of nc-objects.
Outpout: Removes kmer rows for each region in a region_L and returns the 
'surviving' regions.
'''

def remove_overlapping_kmers_from_region_list(region_list):

    new_list = []

    for region in region_list:
        region.remove_kmer_overlap()
        new_list.append(region)
    return new_list


'''
'''
def calc_score_over_region_list(region_list, curve, rotated):
    new_list = []
    for region in region_list:
        calc_score_for_NC_region(curve, region, rotated)
        new_list.append(region)
    return new_list




'''
Filters a regions_list by region sequence length
'''

def filter_regions_by_length(regions_list, length):
    new_region_list = []

    for region in regions_list:
        if region.length > length:
            new_region_list.append(region)
    return new_region_list



'''
Sorts a regions_list by each region objects score attribute.
Note:
THIS FUNCTION OPERATES ON THE INPUTED region_list itself.
It doesnt return anything.
'''
def sort_regions_by_score(regions_list):
    regions_list.sort(key=lambda x: x.cgc_val, reverse = True)
    #regions_list.sort(reversed = True)



'''
Function that makes each NC_region compute its kmer stats and properties.
'''
def calc_kmers_from_region_list(region_list, k_array):
    new_list = []
    for region in region_list:
        for k in k_array:
            region.add_kmer_counts(k)
        new_list.append(region)
    return new_list






'''
'''
def filter_region_list_by_kmer_occurence(region_list, n):
    new_list = []
    for region in region_list:
        region.filter_kmer_by_occurence(n)
        new_list.append(region)

    # take our all regions with "empty" kmer_info
    new_new_list = filter_empty_kmer_regions(new_list)

    final_list = filter_out_empty_kmer_key_in_region_list(new_new_list)
    return final_list






'''
Filters out / removes regions that have a completly empty kmer info content. 
Bascaially, if all keys have a zero size value  
then this regions gets removed. 
'''
def filter_empty_kmer_regions(region_list):
    new_list = []

    for region in region_list:
        if not has_empty_kmer_info(region):
            new_list.append(region)
    return new_list




'''
'''
def plot_region_list(region_list, curve, rotated = False):

    region_intervals = []

    if rotated:
        for region in region_list:
            this_interval = (region.max_relative_start, region.max_relative_stop)
            region_intervals.append(this_interval)
    else:
        for region in region_list:
            this_interval = (region.start, region.stop)
            region_intervals.append(this_interval)

    regions_pos_set = interval_list_to_position_set(region_intervals)
    region_pos_list = list(regions_pos_set)
    region_pos_list.sort()

    relevant_pos_list = np.array(region_pos_list)

    relevant_curve_point = curve[relevant_pos_list]


    plt.figure(figsize=[40,10])
    plt.plot(curve)
    plt.plot(relevant_pos_list, relevant_curve_point, 'x')
    plt.title(str(len(region_list)))





'''
'''
def print_region_list_kmer_info(region_list):
    i = 0
    for region in region_list:
        print('region:', i, '---',' score: ', region.cgc_val, '---- pos: ', region.start, '---- max_relative_start_pos', region.max_relative_start)

        for key, value in region.kmer_info.items():
            print('\tk=', key,)
            for thing in value:
                print('\t', thing[0], ' - ', thing[1], ' - ', 'density:', thing[2]) 

        print('\t-------')

        i = i+1




def add_max_relative_position(region_list, genome_length, max_offset):
    new_list = []
    for region in region_list:
        region.add_max_relative_pos(max_offset, genome_length)
        new_list.append(region)

    return new_list






# TRICKY
'''
Incorporates the relative positions of each nc region
from max rotated fasta

input:
    1 - NC-INTERVALS LIST. list with interval touples specifing

'''
def get_phased_nc_region_list(nc_intervals, og_fasta, max_offset, score_curve):

    # Create and get the region_list
    nc_objects = nc_intervals_to_nc_objects(nc_intervals, og_fasta)

    # Add the max_rotated positions to the regoin leists
    phased_nc_objects = add_max_relative_position(nc_objects, len(og_fasta), max_offset)


    max_scored_nc_objects = calc_score_over_region_list(phased_nc_objects, score_curve, rotated = True)


    return max_scored_nc_objects




'''
Addsa density calculation for each ncregion in a region list
'''
def calc_density_for_region_list(region_list):
    new_list = []

    for region in region_list:
        region.calc_kmer_density()
        new_list.append(region)
    return new_list



'''
self explanitory
'''
def filter_out_empty_kmer_key_in_region_list(region_list):

    new_list = []

    for region in region_list:
        region.filter_out_empty_kmer_lists_in_kmer_dict()
        new_list.append(region)
    return new_list



'''
Sortednc kmers in a region list.
NOTE: The objects are notsorted in the region list
'''

def sort_region_list_on_density(region_list):

    new_list = []

    for region in region_list:
        region.sort_kmer_info_by_density()
        new_list.append(region)
    return new_list



'''
Removes kmer rows for each region in a region_List
'''

def remove_overlapping_kmers_from_region_list(region_list):

    new_list = []

    for region in region_list:
        region.remove_kmer_overlap()
        new_list.append(region)
    return new_list
