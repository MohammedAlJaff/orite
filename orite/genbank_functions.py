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
    
    # Extract and store genome length
    length = get_length_genbank(file_path)
    
    # Raw coding region start and stop position in sting format.
    raw_p = raw_positions_strings_from_genbank(file_path)

    # Split into 
    raw_plus, raw_neg = split_strands(raw_p)

    # 
    valid_plus = make_into_valid_pos(raw_plus)
    valid_neg = make_into_valid_pos(raw_neg)


    # Compute the inbetween non coding interval start and stop positions 
    # given the coding interval positions. 
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




'''
Extracts intervals that are non-coding on both strands

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




'''
COMMENT! 
Input: list of all non-coding positions.
Output: List of tuples with non-coding intervals
'''
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




'''
Returns the length of the organism genome 'contained' in the genbank file

input: 
    Genbank path (string)
output:
    Genome length/Size (int)
'''
def get_length_genbank(file_path):
    recs = [rec for rec in SeqIO.parse(file_path, "genbank")]
    for rec in recs:
        length = len(rec.seq)
    return length




'''
COMMENT
'''
def extract_seq_from_non_coding_intervals(non_coding_intervals, seq):

    non_coding_seqs = []

    for interval in non_coding_intervals:

        non_coding_seqs.append(seq[interval[0]:interval[1]])

    return non_coding_seqs

