from __future__ import division
from decimal import *
from Bio import SeqIO
import subprocess
import scipy
import sys
from scipy import stats
import numpy as np
import math

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import time
from matplotlib import rcParams


#    This file is part of 3D TMH Complexity.
#
#    3D TMH Complexity is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    3D TMH Complexity is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with 3D TMH Complexity.  If not, see <http://www.gnu.org/licenses/>.

print("This file is part of 3D TMH Complexity.")

print("3D TMH Complexity is free software: you can redistribute it and/or modify")
print("it under the terms of the GNU General Public License as published by")
print("the Free Software Foundation, either version 3 of the License, or")
print("(at your option) any later version.")

print("3D TMH Complexity is distributed in the hope that it will be useful,")
print("but WITHOUT ANY WARRANTY; without even the implied warranty of")
print("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the")
print("GNU General Public License for more details.")

print("You should have received a copy of the GNU General Public License")
print("along with 3D TMH Complexity.  If not, see <http://www.gnu.org/licenses/>.")

print(
    "\n\n\nUSAGE:python3 [Script] [file] [TMH number for graph] [OPTIONAL TMH filter number]\nEXAMPLE:python3 Uniprot_to_graphs_and_stats.py file.txt 7 7")
# To do:
# Separate hydrophobicity calculation - more efficient and flexible
# Length restrictions and statistics for final results

# Input files should be obtained in text format downloaded from Uniprot
# and moved to the same directory as this script.

input_filenames = [
    str(sys.argv[1])
]


try:
    tmh_filter_number = sys.argv[3]
    try:
        test = int(sys.argv[3]) + 1
        tmh_filter_number = int(sys.argv[3])
        tmh_filter = True
        print("TMH filter set. Only considering proteins with",
              tmh_filter_number, " tmhs")
    except TypeError:
        print("length filter must me an integer.")
        tmh_filter = False
except IndexError:
    print("no tmh filter set.")
    tmh_filter = False


# Parameters for tmh allowances
minimum_tmd_length = 16
maximum_tmd_length = 38
feature_type = "TRANSMEM"
alternative_feature = "INTRAMEM"

# Several dictionaries of hydrophobicity values are needed to verify findings.

white_and_wimley_hydrophobicity = {
    'type': ' White and Wimley',
    'A': 0.33,
    'C': 0.22,
    'D': 2.41,
    'E': 1.61,
    'F': -0.58,
    'G': 0.12,
    'H': 1.37,
    'I': -0.81,
    'K': 1.81,
    'L': -0.69,
    'M': -0.44,
    'N': 0.43,
    'P': -0.31,
    'Q': 0.19,
    'R': 1.00,
    'S': 0.33,
    'T': 0.11,
    'V': -0.53,
    'W': -0.24,
    'Y': 0.23,
}

eisenberg_hydrophobicity = {
    'type': "Eisenberg",
    'A': 0.620,
    'C': 0.290,
    'D': -0.900,
    'E': -0.740,
    'F': 1.190,
    'G': 0.480,
    'H': -0.400,
    'I': 1.380,
    'K': -1.500,
    'L': 1.060,
    'M': 0.640,
    'N': -0.900,
    'P': 0.120,
    'Q': -0.850,
    'R': -2.530,
    'S': -0.180,
    'T': -0.050,
    'V': 1.080,
    'W': 0.810,
    'Y': 0.260,
}

kyte_doolittle_hydrophobicity = {
    'type': "Kyte & Doolittle",
    'A'	: 1.800,
    'C'	: 2.500,
    'D'	: -3.500,
    'E'	: -3.500,
    'F'	: 2.800,
    'G'	: -0.400,
    'H'	: -3.200,
    'I'	: 4.500,
    'K'	: -3.900,
    'L'	: 3.800,
    'M'	: 1.900,
    'N'	: -3.500,
    'P'	: -1.600,
    'Q'	: -3.500,
    'R'	: -4.500,
    'S'	: -0.800,
    'T'	: -0.700,
    'V'	: 4.200,
    'W'	: -0.900,
    'Y'	: -1.300,
}

### Hydrophobicity calculation code###


def hydrophobicity_calculation(sequence, tmh_locations):
    '''
    Returns a list of lists, and a list of the hydrophocity scales used in the list of lists.
    Each child lists contains a further list of hydrophobicity values
    ordered by tmh number. The parent list is used to split the child lists up according to
    the different hydrophobicity scales.
    '''

    list_of_hydrophobicity_dictionaries = [
        kyte_doolittle_hydrophobicity, eisenberg_hydrophobicity, white_and_wimley_hydrophobicity]
    list_of_tmhs_from_different_hyrodobicity_dictionaries = []
    list_of_dictionary_types = []

    # tmh locations exists as a string in the format: start,end start,end
    tmh_locations = tmh_locations.split(" ")
    for hydrophobicity_type in list_of_hydrophobicity_dictionaries:
        list_of_dictionary_types.append(hydrophobicity_type.get("type"))
        list_of_hydrophobicity_by_tmh_number = []
        for n, i in enumerate(tmh_locations):
            location = i.split(",")
            if len(location) == 2:
                # Some sequence positions given are reported as <0 or >N. These cases
                # will be treated as length restrictions.
                if '<' in str(location):
                    list_of_hydrophobicity_by_tmh_number.append("null")
                elif '>' in str(location):
                    list_of_hydrophobicity_by_tmh_number.append("null")
                else:
                    start = int(location[0])
                    end = int(location[1])

                    # location start and end are human readable locations,
                    # whereas the slices start counting at 0.
                    tmh_sequence = sequence[start - 1:end - 1]

                    # Length restriction and 'X' residues are counted together.
                    if len(tmh_sequence) < maximum_tmd_length and len(tmh_sequence) > minimum_tmd_length and 'X' not in tmh_sequence:
                        hydrophobicity_values_in_tmh = []
                        for residue in tmh_sequence:
                            hydrophobicity_values_in_tmh.append(
                                hydrophobicity_type.get(residue))
                        print(hydrophobicity_values_in_tmh)
                        try:
                            for i in hydrophobicity_values_in_tmh:
                                test_if_number=i+1


                            tmh_hydrophobicity = np.mean(
                                hydrophobicity_values_in_tmh)
                            list_of_hydrophobicity_by_tmh_number.append(
                                tmh_hydrophobicity)
                        except(TypeError):
                            # There was no number!
                            list_of_hydrophobicity_by_tmh_number.append("null")
                    else:
                        list_of_hydrophobicity_by_tmh_number.append("null")
        list_of_tmhs_from_different_hyrodobicity_dictionaries.append(
            list_of_hydrophobicity_by_tmh_number)
    return(list_of_tmhs_from_different_hyrodobicity_dictionaries, list_of_dictionary_types)

# The bahadur statistic allows statistical comparison of P-values.


def bahadur_statistic(p_value, sample_size):
    bahadur_value = (abs(np.log(p_value))) / sample_size
    return(bahadur_value)

# Shannon entropy term of a string

def entropy(string):
    '''
    Use the following code for a custom command. via "Shannon's entropy equation is the standard method of calculation. Here is a simple implementation in Python, shamelessly copied from the Revelation codebase, and thus GPL licensed:"
    '''
    # get probability of chars in string
    prob = [float(string.count(c)) / len(string)
            for c in dict.fromkeys(list(string))]
    # calculate the entropy
    entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])
    return entropy


# Violin plots for datasets.
rcParams.update({'figure.autolayout': True})


def violin_plot(dataset, name, input_file):
    '''
    Gausian Kernel Density Estimation
    Violin plots are similar to histograms and box plots in that they show
    an abstract representation of the probability distribution of the
    sample. Rather than showing counts of data points that fall into bins
    or order statistics, violin plots use kernel density estimation (KDE) to
    compute an empirical distribution of the sample. That computation
    is controlled by several parameters. This example demonstrates how to
    modify the number of points at which the KDE is evaluated (``points``)
    and how to modify the band-width of the KDE (``bw_method``).

    For more information on violin plots and KDE, the scikit-learn docs
    have a great section: http://scikit-learn.org/stable/modules/density.html
    '''

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # plot violin plot
    parts = axes.violinplot(dataset, showmeans=True,
                            showmedians=True, showextrema=False)
    # axes[0].set_title('violin plot')
    # Colour palette is colour blind friendly according to Wong B 2011 Points
    # of view: Color blindness Nature Methods 8:441.
    for pc in parts['bodies']:
        pc.set_facecolor('#56b4e9')
        pc.set_alpha(0.5)
        pc.set_edgecolor('black')
        pc.set_linewidth(0.2)

    parts['cmedians'].set_color('black')
    parts['cmeans'].set_color('#e69f00')

    axes.yaxis.grid(False)
    axes.set_xticks([y + 1 for y in range(len(dataset))])
    axes.set_xlabel('TMH number')
    axes.set_ylabel(name)

    # add x-tick labels
    plt.setp(axes, xticks=[y + 1 for y in range(len(dataset))],
             xticklabels=[x + 1 for x in range(len(dataset))])
    # plt.show()

    # Save bitmap and vector images of figure, then flush figure.
    date_time = time.strftime("%H_%M_%S__%d_%m_%Y")
    file_name = str(date_time + "_" + input_file + "_" + name + ".pdf")
    plt.savefig(file_name, dpi=700)
    file_name = str(date_time + "_" + input_file + "_" + name + ".png")
    plt.savefig(file_name, dpi=700)
    plt.close(fig)

#### Length sorting for TMHs ####


def length_sorting(sequence, tmh_locations):
    '''
    Returns a list with values of complexities of TMHs.
    '''
    list_of_tmh_features = []

    tmh_locations = tmh_locations.split(" ")

    for n, i in enumerate(tmh_locations):
        location = i.split(",")
        if len(location) == 2:
            # print(str(location))
            if '<' in str(location):
                list_of_tmh_features.append("null")
            elif '>' in str(location):
                list_of_tmh_features.append("null")
            elif 'Unknown' in str(location):
                list_of_tmh_features.append("null")
            else:
                start = int(location[0])
                end = int(location[1])

                # Lenght restrictions-these should obey parameters given at
                # begining of script.
                if abs(start - end) < minimum_tmd_length or abs(start - end) > maximum_tmd_length:
                    list_of_tmh_features.append(str("null"))
                else:
                    list_of_tmh_features.append(abs(start - end))

    return(list_of_tmh_features)


def tmsoc_calculation(sequence, tmh_locations):
    '''
    Returns a list with values of complexities of TMHs.
    '''

    list_of_tmh_features = []

    # Generates the TMSOC input files from the source file.
    with open("TMsegments.txt", 'w') as tm_segments_file:
        tm_segments_file.write(tmh_locations)

    with open("sequence.fasta", 'w') as fasta_file:
        fasta_file.write(">placeholderheader\n")
        fasta_file.write(str(sequence))
    # Runs TMSOC.
    perl_script_output = subprocess.check_output(
        ["perl", "TMSOC.pl", "sequence.fasta", "TMsegments.txt"])
    # Processes the TMSOC input files. Normally the TMSOC output would be a
    # multi-line file, not a string.
    list_of_output_lines = str(perl_script_output.decode('UTF-8'))
    list_of_output_lines = list_of_output_lines.replace("\n", ":")
    list_of_output_lines = list_of_output_lines.split(":")
    for line in list_of_output_lines:
        line = str(line.replace(',', ";"))
        line = line.split(";")
        # First we can ignore the "junk" lines. This 7 is not referring to the
        # number of TMHs, but the number of items in the TMSOC output.
        if len(line) == 7:
            tmh_sequence = line[0]
            start_position = int(line[1])
            end_position = int(line[2])
            complexity_score = float(line[3])
            hydrophobicity_score = float(line[4])
            z_score = line[5]
            complexity = line[6]
            # Lenght restrictions-these should obey parameters given at
            # begining of script.
            if abs(start_position - end_position) < minimum_tmd_length or abs(start_position - end_position) > maximum_tmd_length:
                complexity_score = str("null")
            else:
                list_of_tmh_features.append(complexity_score)

    return(list_of_tmh_features)


# Before we conduct the analysis we must determine how many empty TMH
# lists to make. For example if the protein with the most TMHs in the
# dataset has 13 TMHs, then there will be 13 empty lists made in a list
# representing the 1-13 possible TMH locations.
for input_file in input_filenames:
    max_tmd_to_print = int(sys.argv[2])
    # These are the parameters used by the biopython Seq.IO module
    filename = input_file
    input_format = "swiss"

    # We need to check against nearby features to prevent overlapping
    # flanking regions. Note here we want to avoid clashing with INTRAMEM
    # regions, since their flanking regions may be similar, however INTRAMEM regions
    # should not be included in the logged transmembrane features since
    # they may have extremely short transmembrane sequences which would
    # disrupt very much so the alignment of the flanking regions.

    # We iterate through each record, parsed by biopython to find the most
    # number of TMDs. This avoids list index exceeded errors later on.
    tmd_count = 0
    for record in SeqIO.parse(filename, input_format):
        this_record_tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                this_record_tmd_count = this_record_tmd_count + 1
        if this_record_tmd_count > tmd_count:
            tmd_count = this_record_tmd_count
            if this_record_tmd_count > max_tmd_to_print:
                print(record.id, " contained more than ",
                      max_tmd_to_print, "TMHs.")

    # Generate a list of empty lists, one for each tmh set.
    print("Maximum tmh count in", input_file, "is", tmd_count)

    #Generates the empty lists of scores
    list_of_complexity_scores_in_tmh = []
    list_of_lengths_in_tmh = []
    list_of_hydrophobicity_scores_in_tmh = [[]]
    list_of_entropy_scores_in_tmh=[]

    # Fills the empty lists of scores with empty lists, one for each tmh.
    for n in range(tmd_count):
        list_of_complexity_scores_in_tmh.append([])

    for n in range(tmd_count):
        list_of_lengths_in_tmh.append([])

    for n in range(tmd_count):
        list_of_entropy_scores_in_tmh.append([])


    # The hydrophobicity needs to be done separately since it's also dependent
    # on the number of scales being used. ### ADDRESS THIS ISSUE. IT CAN BE
    # OPTIMISED.###

    # Now we can iterate through the records inserting complexity scores into
    # the empty lists.
    for record in SeqIO.parse(filename, input_format):
        transmembrane_record = False
        fuzzy_records = False
        this_record_tmd_count = 0

        # Sequence fasta file
        sequence = record.seq

        # TMH positions file
        tmh_positions = str("")
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                this_record_tmd_count = this_record_tmd_count + 1
                tmh_positions = tmh_positions + \
                    (str(f.location.start) + "," + str(f.location.end) + " ")
                # The -1 is due to a mismatch between the biopython sequence in the record numbering and the pythonic 0 base counting.

                try:
                    # We need to test to see if there are numbers in the position values. This explicit check will help save time removing bad entries.
                    f.location.start+1
                    f.location.end+1
                    tmh_sequence = sequence[f.location.start -
                                            1:f.location.end - 1]
                    if "<" in str(tmh_sequence) or ">" in str(tmh_sequence) or "X" in str(tmh_sequence):
                        fuzzy_records = True
                except(TypeError):
                    print("None integer/float detected in positions in", record.id)

        if tmh_filter == True and int(this_record_tmd_count) == int(tmh_filter_number) and fuzzy_records == False:
            transmembrane_record = True
        elif tmh_filter == False and fuzzy_records == False:
            transmembrane_record = True
        else:
            pass

        # Avoids adding uncleared scores to lists if no TM regions were in
        # protein record.

        if transmembrane_record == True:
            # Adds the complexity score of a helix at (for example the 4th helix)
            # to the complexity list of lists (for example in the 4th position)
            record_complexity_scores = tmsoc_calculation(
                sequence, tmh_positions)
            for n, i in enumerate(record_complexity_scores):
                # Null entries are added for TMHs that are not within length
                # restrictions.
                if i == "null":
                    pass
                else:
                    list_of_complexity_scores_in_tmh[n].append(i)

            #This is the sequence entropy by itself. The complexity is specifically the TMSOC z-score.
            #Any non-TMH features that might get through here won't because the filter has already removed records with non-conforming tmhs.
            record_entropy_scores=[]
            for i, f in enumerate(record.features):
                if f.type == feature_type:
                    # location start and end are human readable locations,
                    # whereas the slices start counting at 0.
                    feature_entropy_score=entropy(record.seq[f.location.start-1:f.location.end-1])
                    record_entropy_scores.append(feature_entropy_score)
            #adding entropy scores to list of lists.
            for tmh_number, i in enumerate(record_entropy_scores):
                if i == "null":
                    pass
                else:
                    list_of_entropy_scores_in_tmh[tmh_number].append(i)

            # Adds the lengths of a helix at (for example the 4th helix)
            # to the complexity list of lists (for example in the 4th position)
            record_lengths = length_sorting(sequence, tmh_positions)
            for n, i in enumerate(record_lengths):
                # Null entries are added for TMHs that are not within length
                # restrictions.
                if i == "null":
                    print("incorrect length found after filtering in ", str(record.id))
                    pass
                else:
                    list_of_lengths_in_tmh[n].append(i)

            # Adds to the hydrophobicity list of lists. The structure is different
            # to complexity since there are 3 different hydrophobicity scales.
            hydrophobicity_for_record = hydrophobicity_calculation(
                sequence, tmh_positions)
            for scale_number, hydrophobicity_scale in enumerate(range(len(hydrophobicity_for_record[1]))):
                list_of_hydrophobicity_scores_in_tmh.append([])
                for n in range(tmd_count):
                    list_of_hydrophobicity_scores_in_tmh[scale_number].append([
                    ])
                for tmh_number, i in enumerate(hydrophobicity_for_record[0][scale_number]):
                    if i == "null":
                        pass
                    else:
                        list_of_hydrophobicity_scores_in_tmh[scale_number][tmh_number].append(
                            i)

    # Graphs
    violin_plot(list_of_complexity_scores_in_tmh[
                0:max_tmd_to_print], str("TMSOC z-score"), input_file)
    violin_plot(list_of_lengths_in_tmh[0:max_tmd_to_print], str(
        "Length"), input_file)
    violin_plot(list_of_entropy_scores_in_tmh[0:max_tmd_to_print], str(
        "Sequence Entropy"), input_file)
    for scale_number, scales in enumerate(list_of_hydrophobicity_scores_in_tmh[0:len(hydrophobicity_for_record[1])]):
        violin_plot(list_of_hydrophobicity_scores_in_tmh[scale_number][0:max_tmd_to_print], str(
            hydrophobicity_for_record[1][scale_number] + " Hydrophobiciity Scale"), input_file)






    stat_tests_list = [scipy.stats.kruskal, scipy.stats.ks_2samp]

    print("\n\n", "###", input_file, "###\n")

    for stat_tests in stat_tests_list:

        # stats for complexity
        for n, i in enumerate(list_of_complexity_scores_in_tmh):

            # n is the index, so for human readable numbers we need to add 1. i.e
            # the first helix is n=0, so we report it as n+1.
            print("TMH ", n + 1)
            # print(i)
            print("Mean complexity:", np.mean(i), ", N:", len(i))
            # list_of_complexity_scores_in_tmh = [
            #    x for x in list_of_hydrophobicity_scores_in_tmh if x != []]

            # Only calculating stats if the sample size is large enough to be worth it.
            sample_size_threshold = 10
            if len(i) > sample_size_threshold:
                if n + 1 < len(list_of_complexity_scores_in_tmh):
                    for x in range(len(list_of_complexity_scores_in_tmh)):
                        if x > n and len(list_of_complexity_scores_in_tmh[x]) > sample_size_threshold:
                            sample_size_for_bahadur = len(
                                list_of_complexity_scores_in_tmh[x]) + len(i)
                            statistic_test_result = stat_tests(
                                list_of_complexity_scores_in_tmh[n], list_of_complexity_scores_in_tmh[x])
                            print(sample_size_for_bahadur,
                                  statistic_test_result[1])
                            print("TMH ", n + 1, " to ", x + 1, ":",  statistic_test_result, ", Bahadur slope=",
                                  bahadur_statistic(statistic_test_result[1], sample_size_for_bahadur))

        print("\n")

        # stats for entropy
        for n, i in enumerate(list_of_entropy_scores_in_tmh):

            # n is the index, so for human readable numbers we need to add 1. i.e
            # the first helix is n=0, so we report it as n+1.
            print("TMH ", n + 1)
            # print(i)
            print("Mean entropy:", np.mean(i), ", N:", len(i))
            # list_of_entropy_scores_in_tmh = [
            #    x for x in list_of_hydrophobicity_scores_in_tmh if x != []]

            # Only calculating stats if the sample size is large enough to be worth it.
            sample_size_threshold = 10
            if len(i) > sample_size_threshold:
                if n + 1 < len(list_of_entropy_scores_in_tmh):
                    for x in range(len(list_of_entropy_scores_in_tmh)):
                        if x > n and len(list_of_entropy_scores_in_tmh[x]) > sample_size_threshold:
                            sample_size_for_bahadur = len(
                                list_of_entropy_scores_in_tmh[x]) + len(i)
                            statistic_test_result = stat_tests(
                                list_of_entropy_scores_in_tmh[n], list_of_entropy_scores_in_tmh[x])
                            print(sample_size_for_bahadur,
                                  statistic_test_result[1])
                            print("TMH ", n + 1, " to ", x + 1, ":",  statistic_test_result, ", Bahadur slope=",
                                  bahadur_statistic(statistic_test_result[1], sample_size_for_bahadur))

        print("\n")

        # Stats for lengths
        for n, i in enumerate(list_of_lengths_in_tmh):
            print("TMH ", n + 1)
            print("Mean length:", np.mean(i), ", N:", len(i))

            if len(i) > 10:
                if n + 1 < len(list_of_lengths_in_tmh):
                    print("TMH ", n + 1, " to ", n + 2, ":", stat_tests(
                        list_of_lengths_in_tmh[n], list_of_lengths_in_tmh[n + 1]))

        print("\n")

        # stats for hydrophobicity.
        # Remove empty lists created previously as a one liner.
        list_of_hydrophobicity_scores_in_tmh = [
            x for x in list_of_hydrophobicity_scores_in_tmh if x != []]
        for scale_number, scales in enumerate(list_of_hydrophobicity_scores_in_tmh):
            print("Hydrophobicity scale:",
                  hydrophobicity_for_record[1][scale_number])

            no_empty_tmh_list_of_hydrophobicity_scores_in_tmh = [
                x for x in list_of_hydrophobicity_scores_in_tmh[scale_number] if x != []]
            for n, i in enumerate(no_empty_tmh_list_of_hydrophobicity_scores_in_tmh):
                # n is the index, so for human readable numbers we need to add 1. i.e
                # the first helix is n=0, so we report it as n+1.
                print("TMH ", n + 1)
                print("Mean Hydrophobicity:", np.mean(i), ", N:", len(i))
                #print(i)
                if len(i) > 10:
                    if n + 1 < len(no_empty_tmh_list_of_hydrophobicity_scores_in_tmh):
                        print("TMH ", n + 1, " to ", n + 2, ":", stat_tests(
                            no_empty_tmh_list_of_hydrophobicity_scores_in_tmh[n], no_empty_tmh_list_of_hydrophobicity_scores_in_tmh[n + 1]))
                print("\n")
