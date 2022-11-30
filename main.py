import re
import string
import csv
import os
import argparse
import shutil
import matplotlib.pyplot as plt
import sys
'''
FUNCTIONS FOR GETTING INPUTS
'''
def parse_args():
    parser = argparse.ArgumentParser(description="Analysis program.")
    parser.add_argument("input", help="Directory containing all input files")
    parser.add_argument("--out", help="Directory to place all output files", default="output")
    args = parser.parse_args()

    # Check if output directory is safe to delete
    if not output_dir(args.out):
        print("Halting...")
        exit()

    return args

'''
FUNCTIONS FOR PARSING FILE
'''

def get_meta(line):
    '''Gets protein name from first line of text.
    Params:
        line (str): A single line (newline separated) from input
    Return:
        (str): id
    '''

    # Basic checks to see if file is formatted correctly
    assert(line[0] == '>')
    

    # Look for the string between > and the first whitespace or till end of string
    search_object = re.search(r'(?<=\>).*?(?=\s+|$)', line)
    assert(search_object != None)
    fid = str(search_object.group())

    return fid

def line_range(line):
    '''Gets line range string (e.g. 123-234) from a line
    Params:
        line (str): A single line (newline separated) from input
    Return:
        str: line range string
    '''
    return re.search(r'\d+-\d+', line)

def get_letters(line):
    '''Extracts all lower case or upper case letters from a line
    Params:
        line (str): A single line (newline separated) from input
    Return:
        str: all alphabetic characters in the single line
    '''
    search_object = re.search(r'[a-zA-Z]+', line)

    if search_object:
        return search_object.group()

    return ""

def add_current_line(logical_lines, current_line):
    '''Add line to list
    Params:
        logical_lines (list[str]): List of all "logical" lines
        current_line (str): Logical line to add into the list
    Return:
        (list[str], str): List with new line appended and an empty current_line
    '''
    if current_line != "":
        if len(logical_lines) > 0 and is_upper(logical_lines[-1]) == is_upper(current_line):
            # If last case same as current case, just append it (same logical line)
            logical_lines[-1] += current_line
        else:
            logical_lines.append(current_line)

    current_line = ""
    return logical_lines, current_line

def parse_file(filepath, verbose=False):
    '''Get Protein Name, list of lines, and list of ranges from a file
    Params:
        filepath (str): Path to input file
        verbose (bool): Flag for outputting log messages
    Return:
        (str, list[str], list[str]): Protein Name, list of lines, list of ranges
    '''
    
    if verbose:
        print('Reading file "{}"...'.format(filepath))


    # Read all raw files
    with open(filepath, 'r') as f:
        f_raw = f.read()
        lines = f_raw.splitlines()

    # Get Protein Name
    fid = get_meta(lines[0])

    # Get logical lines and ranges
    logical_lines, ranges = [], []
    current_line = ""

    for index, line in enumerate(lines[1:]):
        # New logical line and range
        if line_range(line) != None:
            logical_lines, current_line = add_current_line(logical_lines, current_line)
            ranges.append(line_range(line)[0])

        # Add line to current line
        current_line += get_letters(line)

    # Add the final line and range
    logical_lines, current_line = add_current_line(logical_lines, current_line)

    return fid, logical_lines, ranges

'''
FUNCTIONS TO EXTRACT INFO
'''
def is_upper(line):
    '''Gets whether line is upper case or lower case
    Params:
        line (str): A single logical line
    Return:
        bool: True if line is upper case, False if line is lower case
    '''
    upper = line[0].isupper()

    # Double check that the entire line is the same case
    for char in line:
        assert(char.isupper() == upper)

    return upper

def edge_exception(line_upper, line_lower_1, line_lower_2):
    """Criteria 1 Exception: terminus line is both a high-complexity segment and longer than 10 amino acids
    Condition: 
        HCR/(HCR+LCR) <= 50% 
        Top 4 amino acids within the HCR+LCR stretch >= 50% of the entire amino acids composition of the combined segments.
    Return:
        top_aa (float): aggregate percentage of top 4 amino acids
        int: 0 if failed condotion
    """
    if (len(line_upper)/(len(line_lower_1) + len(line_lower_2) + len(line_upper))) <= 0.5:
        combo_dict = dict.fromkeys((line_upper + line_lower_1 + line_lower_2).lower(),0)
        for i in (line_upper + line_lower_1 + line_lower_2).lower():
            combo_dict[i] += 1
        num = sum(sorted(combo_dict.values(), reverse=True)[:4])
        top_aa = num/(len(line_lower_1) + len(line_lower_2) + len(line_upper))
        return top_aa
    else:
        return 0 

def criteria1_match(motif, line_1, line_2, line_3):
    '''Check termini cases
    Params:
        upper_lower_upper (str): trio segments motif
        line_1 (str): A string for segment 1
        line_2 (str): A string for segment 2
        line_3 (str): A string for segment 3
    Condition:
        HCR < 10 amino acids (upper case letters)
        LCR next to it >= 12 amino acids (lower case letters)
    Return:
        bool: True if satisfy condition, False if does not match criteria 1
        str: append line 1 with either line 2 or line 3
        "" (str): If does not match criteria 1
        int: 1 if HCR is N-terminus, 3 if HCR is C-terminus, None if does not match criteria 1
    '''

    # HCR - LCR - HCR
    if motif == "upper_lower_upper":

        if (is_upper(line_2) or not len(line_2) >= 12):
            return False, "", None, None

        # HCR-LCR-Null
        if line_3 == None:
            return True,line_2, 2, 1
        else:
            if len(line_1) <= 10:
                return True, line_1 + line_2, 1, 1

            if len(line_1) > 10:
                if edge_exception(line_1, line_2, '') >= 0.5:
                    return True, line_1 + line_2, 1, 1
                else:
                    return False, "", None, None

        # Null-LCR-HCR
        if line_3 != None:
            if len(line_3) <= 10:
                return True, line_2 + line_3, 3, 1
            elif len(line_3) > 10:
                if edge_exception(line_3, line_2, "") >= 0.5:
                    return True, line_2 + line_3, 3, 1
                else:
                    return False, "", None, None

    # LCR - HCR - LCR
    elif motif == "lower_upper_lower":

        #if is_upper(line_2) and (not len(line_2) <= 10):
        #    return False, "", None

        if line_3 == None:
            if len(line_1) >= 12 and len(line_2) <= 10:
                return True, line_1 + line_2, 1, 1

            elif len(line_1) >= 12 and len(line_2) > 10:
                if edge_exception(line_2, line_1, "") >= 0.5:
                    return True, line_1 + line_2, 1, 1
                else:
                    return False, "", None, None

        # This is a special case that return 0 for the index_1, since we only want to get the index of the first segment and not anything else.
        if len(line_1) >= 12 and len(line_2) > 10:
            return True, line_1, 1, 0

        if line_3 != None and len(line_3) >= 12:
            return True, line_2 + line_3, 3, 1

    return False, "", None, None
def criteria2_match(upper_lower_upper, line_1, line_2, line_3):
    ''' Check Island Case
    Params:
        upper_lower_upper (str): trio segments motif
        line_1 (str): A string for segment 1
        line_2 (str): A string for segment 2
        line_3 (str): A string for segment 3
    Condition:
        LCR >= 12 amino acids (lower case letters)
    Return: 
        bool: True if satisfy condition, False if: 1) line 3 is an empty string or 2) not upper_lower_upper motif or 3) line 2 has less than 12 letters
        str: line 2 if match criteria 2, empty string if failed to match
    '''
    if line_3 == None:
        return False, ""

    if not upper_lower_upper:
        return False, ""

    if len(line_2) < 12:
        return False, ""

    return True, line_2

def criteria3_match(upper_lower_upper, line_1, line_2, line_3):
    ''' Check Main Body Case
    Params:
        upper_lower_upper (str): trio segments motif
        line_1 (str): A string for segment 1
        line_2 (str): A string for segment 2
        line_3 (str): A string for segment 3
    Condition:
        LCR - HCR - LCR
        (LCR >= 12 - HCR <= 10 - LCR >= 12) or (HCR > 10 and < 100, is HCR/(LCR+HCR+LCR) <= 50%)
    Res:
        bool: True if satisfy conditions, False if: 1) line 3 is NULL or 2) not upper_lower_upper or 3) line 2 > 100 letters
        Str: append line 1, line 2, and line 3
        "" (str): failed to match criteria 3

    '''
    if line_3 == None:
        return False, ""

    if upper_lower_upper:
        return False, ""

    if len(line_2) <= 10:
        if len(line_1) >= 12 or len(line_3) >= 12:
            return True, line_1 + line_2 + line_3
    # LCR-HCR-LCR
    else:
        total_len = len(line_1) + len(line_2) + len(line_3)
        if len(line_2) > 70:
            return False, ""
        else:
            if((len(line_2) / total_len <= 0.5) and (edge_exception(line_2, line_1, line_3) >= 0.5)):
                return True, line_1 + line_2 + line_3
            else:
                return False, ""
    
    return False, ""

def select_nonzero_string(strings):
    '''Only look at non-NULL strings
    '''
    for string in strings:
        if string != "":
            return string

    return ""

def select_element_ifexists(logical_lines, index):
    if index + 2 <= len(logical_lines) - 1:
        return logical_lines[index + 2]

    return None



def criteria_match(logical_lines, ranges):
    """ Matching Criterias
    Params:
        logical_lines (list[str]): List of all "logical" lines
        ranges (list[str]): List of ranges
    technique:
        Uses a variation of the window sliding technique 
    Logic:
        Criteria precedence: 3 -> 1 -> 2
        Look at three adjacent segments (line_1 - line_2 - line_3: trio block) and check if block satisfy any criteria:
            If satisfy any criteria:
                trio block is now considered the new line_1 -> look at another set of trio block (new line_1 - line_2 - line_3)
            If does not satisfy any criteria:
                move down a segment: old line_2 -> line_1, old line_3 -> line_2, new line -> line_3
            Iterates through the segments until no more new segments available
    Return:
        out_lines (list[str]): List of segments that satisfied the criterions
        out_range (list[list[str]]): Nested list of strings that denote the positions of the segments
        out_criteria (list[str]): List of the criterions the segments satisfied

    Note: This is the main body of the SEG-Filtered Algorithm
    """
    out_lines, out_ranges, out_criterias = [], [], []

    string_1, string_2, string_3 = "", "", ""
    old_string_1, old_string_2, old_string_3 = "", "", ""
    string_ranges = []
    string_criterias = []

    iterator = iter(range(len(logical_lines) - 1))


    for index in iterator:
        # Store old values (if string_1 needs to be added)
        old_string_1, old_string_2, old_string_3 = string_1, string_2, string_3

        combo_dict = []

        matched_1 =False

        # Get current lines
        line_1 = select_nonzero_string([string_3, string_1, string_2, logical_lines[index]])
        line_2 = logical_lines[index + 1]
        line_3 = select_element_ifexists(logical_lines, index)

        # Get state of lines9
        upper_lower_upper = not is_upper(line_2) and (line_3 == None or is_upper(line_3))
        lower_upper_lower = is_upper(line_2) and (line_3 == None or not is_upper(line_3))
        assert(upper_lower_upper or lower_upper_lower)

        # Check criterias


        # If line_3 == None or (index == 0 and upper_lower_upper):
        if line_3 == None:
            # LCR-HCR-NULL
            if is_upper(line_2):
                    matched_1, string_1, line_matched, index_1 = criteria1_match("lower_upper_lower", line_1, line_2, line_3)
            else:
                # HCR-LCR--NULL
                    matched_1, string_1, line_matched, index_1 = criteria1_match("upper_lower_upper", line_1, line_2, line_3)

        elif index == 0:
            #HCR-LCR--HCR
            if is_upper(line_1): 
                    matched_1, string_1, line_matched, index_1 = criteria1_match("upper_lower_upper", line_1, line_2, line_3)
            else:
                #LCR-HCR-LCR
                    matched_1, string_1, line_matched, index_1 = criteria1_match("lower_upper_lower", line_1, line_2, line_3)

        else:
            matched_1 = False

        matched_2, string_2 = criteria2_match(upper_lower_upper, line_1, line_2, line_3) 
        matched_3, string_3 = criteria3_match(upper_lower_upper, line_1, line_2, line_3) 

        # Extracting information about the matched trio segments: Segments' positions and criteria matched  
        # Note: If criteria match, the segment(s) is not stored in the outputted list right away 
        # The segment(s) is use in the next iteration until either: 1) no new match is found or 2) no more segments available
        if matched_3:
            string_1, string_2 = "", ""

            string_ranges.append(ranges[index + 0])
            string_ranges.append(ranges[index + 1])
            string_ranges.append(ranges[index + 2])
            string_criterias.append(str(3))

            # Skip next 2 iterations (and just use string_3 as line_1)
            next(iterator, 0)
            continue

        if matched_1:
            string_2, string_3 = "", ""

            if line_matched == 1:
                string_ranges.append(ranges[index + 0])
            elif line_matched == 2:
                string_ranges.append(ranges[index + 1])
            elif line_matched == 3:
                string_ranges.append(ranges[index + 2])
                next(iterator, 0)
            string_ranges.append(ranges[index + index_1])
            string_criterias.append(str(1))
            continue

        if matched_2:
            # Since matched_2 only applies to the middle line, previous matches have to be added
            if old_string_1 != "":
                out_lines.append(old_string_1)
                out_ranges.append(list(set(string_ranges)))
                out_criterias.append("+".join(string_criterias))
                
            if old_string_3 != "":
                out_lines.append(old_string_3)
                out_ranges.append(list(set(string_ranges)))
                out_criterias.append("+".join(string_criterias))

            string_ranges = []
            string_criterias = []
            string_1, string_3 = "", ""

            string_ranges.append(ranges[index + 1])
            string_criterias.append(str(2))
            continue
        # Add matches
        for string in [old_string_3, old_string_1, old_string_2]:
            if string != "":
                out_lines.append(string)
                out_ranges.append(list(set(string_ranges)))
                out_criterias.append("+".join(string_criterias))

        # Reset since no criteria matched
        string_1, string_2, string_3 = "", "", ""
        string_ranges = []
        string_criterias = []

    # Add matches
    for string in [string_3, string_1, string_2]:
        if string != "":
            out_lines.append(string)
            out_ranges.append(list(set(string_ranges)))
            out_criterias.append("+".join(string_criterias))
            break

    return out_lines, out_ranges, out_criterias

def get_common_letters(lines):
    '''Count frequency of characters
    Params:
        lines (list[str]): List of logical lines
    Return:
        dict{char, int}: Mapping of character (upper case) to frequency
    '''
    letters_freq = dict.fromkeys(string.ascii_uppercase, 0)

    for line in lines:
        for char in line:
            letters_freq[char.upper()] += 1

    return letters_freq

def sigfig(num, accuracy=0.01):
    '''Round number to some accuracy.
    Params:
        num (int): Number to round
        accuracy (float): The accuracy of the rounded number.
    Return:
        int: Rounded number.
    '''
    multiplier = 1 / accuracy
    return round(num * multiplier) / multiplier

def format_common_letters(letters_freq, num=4):
    '''Sort letters for printing
    Params:
        letters_freq (dict{char, int}): Mapping of character to frequency
        num (int): The number of characters to return
    Return:
        List[Tuple(char, int, int)]: List of the top "num" amount of characters, containing info on the count and relative percentage.
    '''
    # Get total number of chars
    total = 0
    for letter in letters_freq:
        total += letters_freq[letter]

    if total == 0:
        return None

    sorted_freq = sorted(letters_freq.items(), key=lambda item: item[1], reverse=True)

    # Get percentages of chars
    out = []
    for freq in sorted_freq:
        out.append((freq[0], freq[1], sigfig(freq[1]/total)))

    return out[:num]

def get_cluster_filenames(args, common_letters_list):
    '''Get cluster name for putting into directories
    Params:
        common_letters_list (List[dict{char, int}]): List of mappings of character to frequency
    Return:
        List[str]: Gives top 2 letters in a string for each mapping
    '''
    fmt_letters_list = list(map(lambda x : format_common_letters(x, num=2), common_letters_list))
    clusters = []
    for ll in fmt_letters_list:
        top_letters = list(map(lambda x: x[0], ll))
        clusters.append("".join(sorted("".join(top_letters))))

    f_clusters = []
    for cluster_name in clusters:
        f_cluster = os.path.join(args.out, "clusters", cluster_name + ".csv")
        f_clusters.append(f_cluster)

    return f_clusters

'''
FUNCTIONS TO OUTPUT
'''
def y_n_prompt(directory):
    '''Safety mechanism before deleting output folder
    Params:
        directory (str): Directory to be deleted.
    Return:
        bool: True if directory can be deleted, False if it can't
    '''
    while True:
        prompt = 'WARNING: This action will replace all contents that are currently in the folder "{}". Continue? y/n: '.format(directory)
        reply = str(input(prompt)).lower().strip()
        if len(reply) != 1:
            continue
        if reply[0] == 'y':
            print('Clearing directory "{}"...'.format(directory))
            shutil.rmtree(directory)
            return True
        if reply[0] == 'n':
            return False

def create_dir(directory):
    '''Create directory if it doesn't exist
    Params:
        directory (str): Directory to be created.
    '''
    if os.path.exists(directory):
        return
    os.makedirs(directory)

def output_dir(directory):
    '''Creates all the necessary basic directory structures.
    Params:
        directory (str): Base output directory specified by user
    Return:
        bool: True if directories were created successfully, False if there was a fatal error
    '''
    if os.path.exists(directory):
        if not y_n_prompt(directory):
            return False
    os.makedirs(directory)
    os.makedirs(os.path.join(directory, "sequences", "by_id"))
    os.makedirs(os.path.join(directory, "sequences", "by_overall"))

    os.makedirs(os.path.join(directory, "clusters"))

    os.makedirs(os.path.join(directory, "graphs", "by_id"))
    os.makedirs(os.path.join(directory, "graphs", "by_overall"))
    return True

def write_line_csv(filepaths, fid, fmt_letters, strings, ranges, criteria, verbose=False):
    '''Append a new line to the csv file
    Params:
        filepaths (List[str]): List of filepaths to write to
        fid (str): Protein Name
        fmt_letters (List[List[Tuple(str, int, float)]]): List of 4 most common letters and metadata per each string
        strings (List[str]): List of culled strings
        ranges (List[str]): List of ranges for each string
        criteria (List[int]): List of criterias for each string
        verbose (Bool): Whether to list error messages.
    '''
    if len(strings) < 1:
        return

    if verbose:
        print('Writing to file "{}"...'.format(filepath))


    # Add header if file doesn't exist yet
    for filepath in filepaths:
        if os.path.isfile(filepath):
            continue

        with open(filepath, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow([
                "ID", "Sequence", "Ranges", "Criteria",
                "Common 1", "Count 1", "Percent 1",
                "Common 2", "Count 2", "Percent 2",
                "Common 3", "Count 3", "Percent 3",
                "Common 4", "Count 4", "Percent 4",])

    # Append data
    for filepath in filepaths:
        with open(filepath, 'a+', newline='') as outfile:
            for index, letter_dat in enumerate(fmt_letters):
                writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                ranges_str = ",".join(ranges[index])
                fid_indx = "{}.{}".format(fid, index + 1)
                row = [item for sublist in letter_dat for item in sublist]
                writer.writerow([fid_indx, strings[index], ranges_str, criteria[index]] + row)

def write_cluster_line_csv(filepaths, fid, fmt_letters, strings, ranges, criteria, verbose=False):
    '''Append a new line to the csv file
    Params:
        filepaths (List[str]): List of filepaths to write to
        fid (str): Protein Name
        fmt_letters (List[List[Tuple(str, int, float)]]): List of 4 most common letters and metadata per each string
        strings (List[str]): List of culled strings
        ranges (List[str]): List of ranges for each string
        criteria (List[int]): List of criterias for each string
        cluster (Bool): whether filenames are being passed in a cluster
        verbose (Bool): Whether to list error messages.
    '''
    if len(strings) < 1:
        return

    if verbose:
        print('Writing to file "{}"...'.format(filepath))


    # Add header if file doesn't exist yet
    for filepath in filepaths:
        if os.path.isfile(filepath):
            continue

        with open(filepath, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow([
                "ID", "Sequence", "Ranges", "Criteria",
                "Common 1", "Count 1", "Percent 1",
                "Common 2", "Count 2", "Percent 2",
                "Common 3", "Count 3", "Percent 3",
                "Common 4", "Count 4", "Percent 4",])

    # Append data
    for index, filepath in enumerate(filepaths):
        with open(filepath, 'a+', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            ranges_str = ",".join(ranges[index])
            fid_indx = "{}.{}".format(fid, index + 1)
            row = [item for sublist in fmt_letters[index] for item in sublist]
            writer.writerow([fid_indx, strings[index], ranges_str, criteria[index]] + row)



def write_summary_line_csv(filepaths, fid, data, verbose=False):
    '''Append a new line to the csv file for a fid/family count
    Params:
        filepaths (List[str]): List of filepaths to write to
        fid (str): Protein Name
        data (List[List[Tuple(str, int, float)]]): List of 4 most common letters and metadata (wrapped by a list)
        verbose (Bool): Whether to list error messages.
    '''
    if verbose:
        print('Writing to file "{}"...'.format(filepaths))

    # Add header if file doesn't exist yet
    for filepath in filepaths:
        if os.path.isfile(filepath):
            continue

        with open(filepath, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            header = [
                "ID", 
                "Common 1", "Count 1", "Percent 1",
                "Common 2", "Count 2", "Percent 2",
                "Common 3", "Count 3", "Percent 3",
                "Common 4", "Count 4", "Percent 4", ]
            writer.writerow(header)

    for filepath in filepaths:
        # Append data
        with open(filepath, 'a+', newline='') as outfile:
            for index, letter_dat in enumerate(data):
                if letter_dat == None:
                    print("Did not detect any regions of Low-Complexity within the dataset.")
                    sys.exit()
                else:
                    writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    row = [item for sublist in letter_dat for item in sublist]

                    fid_indx = "{}.{}".format(fid, index + 1)
                    row = [fid_indx] + row 
                    writer.writerow(row)


def generate_seq_filenames(args, fid):
    '''Generates all file names for the given sequence.
    Params:
        args (Dict{str, }): Args from parseargs.
        fid (str): Protein Name
    Returns:
        List[str, str, str,...]: List of all filenames for sequences
    '''
    f_overall = os.path.join(args.out, "sequences", "by_overall", "overall.csv")
    f_fid = os.path.join(args.out, "sequences", "by_id",  str(fid) + ".csv")
    return [f_overall, f_fid]

def generate_counts_filenames(args, fid):
    '''Generates all summaries filenames.
    Params:
        args (Dict{str, }): Args from parseargs.
        fid (str): Protein Name
    Returns:
        List[str, str, str,...]: List of all filenames for summaries
    '''
    f_overall = os.path.join(args.out, "sequences", "by_overall", "overall-summary.csv")
    f_fid = os.path.join(args.out, "sequences", "by_id",  str(fid) + "-summary.csv")
    return [f_overall, f_fid]

'''
FUNCTIONS FOR GRAPHING
'''


def bar_graph(rootpath, data, fid=None, family=None, overall=None):
    '''Createts bar graph
    Params:
        rootpath (str): Rootpath of graph directories.
        data (List[Dict{str, int}]): List of mappinigs from letter to count
    '''
    filepath = os.path.join(rootpath)

    # Combine data for bar graph
    if fid != None:
        temp_data = dict.fromkeys(string.ascii_uppercase, 0)
        for index, dat in enumerate(data):
            temp_data = {x: dat.get(x) + temp_data.get(x)
                for x in dat }
        data = temp_data
    else:
        data = data[0]


    if fid != None:
        title = fid
    if overall != None:
        title = overall

    impath = os.path.join(filepath, str(title) + ".jpg")

    fig = plt.figure(figsize = (9, 7))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    sorted_dat = sorted(data.items(), key=lambda x: x[1], reverse=True)
    sorted_dat = list(filter(lambda x: x[1] > 0, sorted_dat))

    sorted_dat_labels = list(map(lambda x: x[0], sorted_dat))
    sorted_dat_data = list(map(lambda x: x[1], sorted_dat))



    ax.set_title(title)
    ax.set_xlabel("Amino acids")
    ax.set_ylabel("Count")
    #ax.set_xticks(list(range(len(sorted_dat))))
    #ax.set_xticklabels(sorted_dat_labels)
    ax.barh(sorted_dat_labels , sorted_dat_data, color = "lightskyblue")
    for index, value in enumerate(sorted_dat_data):
        plt.text(value/2, index, str(value), ha = "center", fontweight = 'bold')
        plt.text(value, index, "{0:.2f}%".format(100 * value/sum(sorted_dat_data)), fontweight = 'bold')
    plt.savefig(impath)
    plt.close()


'''
MAIN
'''
if __name__ == "__main__":
    args = parse_args()

    # Create dict to keep track of most common letters
    overall_common_letters = dict.fromkeys(string.ascii_uppercase, 0)

    for filename in os.listdir(args.input):
        input_fname = os.path.join(args.input, filename)
        try:
            fid, logical_lines, ranges = parse_file(input_fname)
        except:
            print("Invalid file format found: {}".format(filename))
            continue

        # Get culled string and range of letters
        crit_str, crit_range, criteria = criteria_match(logical_lines, ranges)

        if len(crit_str) == 0:
            continue

        # Count for each culled string
        common_letters_list = []
        for culled_str in crit_str:
            common_letters = get_common_letters(culled_str)
            common_letters_list.append(common_letters)

            # Add common letters to family and overall
            overall_common_letters = {x: common_letters.get(x) + overall_common_letters.get(x)
                for x in common_letters }

        # Create culled strings csv
        print("Creating csv for fid", fid)
        seq_filenames = generate_seq_filenames(args, fid)
        fmt_letters_list = list(map(format_common_letters, common_letters_list))
        write_line_csv(seq_filenames, fid, fmt_letters_list, crit_str, crit_range, criteria)

        cluster_names = get_cluster_filenames(args, common_letters_list)
        write_cluster_line_csv(cluster_names, fid, fmt_letters_list, crit_str, crit_range, criteria)

        # Creating graph for id
        print("Creating graph for id", fid)
        f_image = os.path.join(args.out, "graphs", "by_id")
        bar_graph(f_image, common_letters_list, fid=fid)


    print("Creating graph for overall")
    counts_filenames = generate_counts_filenames(args, 0)
    fmt_letters = format_common_letters(overall_common_letters)
    write_summary_line_csv([counts_filenames[0]], -1, [fmt_letters])
    f_image = os.path.join(args.out, "graphs",  "by_overall")
    bar_graph(f_image, [overall_common_letters], overall="overall")

