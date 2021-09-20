# **********************************************************************
#  Copyright notice
#
#  (c) 2015-2021 Dominic Simm <dominic.simm@cs.uni-goettingen.de>
#  All rights reserved.
#
#  This file is part of Odysseus.
#
#  Odysseus is free software: you can redistribute it and/or modify
#  it under the terms of the Creative Commons BY-NC-SA 4.0 License as
#  published by the Creative Commons Corporation, either version 2 of the
#  License, or (at your option) any later version.
#
#  Odysseus is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  Creative Commons BY-NC-SA 4.0 License for more details.
#
#  You should have received a copy of the Creative Commons BY-NC-SA 4.0 License
#  along with Odysseus.  If not, see <https://creativecommons.org/licenses>.
# **********************************************************************

# Changes 2015-2021 by Dominic Simm <dominic.simm@cs.uni-goettingen.de>
# See the ChangeLog or git repository for details.

#
# Created by dsimm on 04/08/15.
# Updated by dsimm on 18/12/17
#

# Inclusions
# Load BioPython packages
from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import CodonUsage, GC

# Load PyParsing package
from pyparsing import *
import pyparsing as pp

# Load general Python packages
from io import StringIO
import csv
import operator
import subprocess

# Load Package for Database support
import psycopg2  # only system on that runs the code

# Load Package for pretty print
import pprint

# Local Inclusions
from .globals import *
from .figure_utilities import *
from .help_utilities import *
from .html_utilities import *
from .markov_chain import *

# Numpy: No scientific number representation
np.set_printoptions(suppress=True)
# Numpy: Set div.zero handling
np.seterr(divide='ignore')


def calculate_abs_codon_usages(dna_fasta_file):
    """Calculate codon-usage of a set of DNA sequences (nucleotides) stored in a passed file
    => Absolute values / total codon abundances

    Args:
        dna_fasta_file (str): Absolute path to FASTA file of DNA sequences

    Returns:
        dict, dict: Absolute codon usages (mono, di-codon usage)
    """
    # Initialize new local codon-lists
    codons, dicodons = initialize_codons()
    first_codon = None
    # Check file existence
    if os.path.isfile(dna_fasta_file):
        for seq_record in SeqIO.parse(dna_fasta_file, "fasta"):
            step = 3  # 1 (all 3 frames) | 3 (specified "coding" frame)
            # print(seq_record.seq)
            for n in range(0, len(seq_record), step):
                first_codon = seq_record.seq[n:n+3]
                sec_codon = seq_record.seq[n+3:n+6]
                # Ignore internal occuring stop_codons
                if first_codon in all_codons and sec_codon in all_codons:
                    codons[str(first_codon)] += 1
                    dicodons[str(first_codon)][str(sec_codon)] += 1
            # Add last codon (stop-codon)
            if first_codon in all_codons:
                codons[str(first_codon)] += 1
    else:
        print("Can't find file at the specified location.")

    return codons, dicodons


def calculate_abs_codon_usages_weighted(dna_fasta_file):
    """Calculate weighted codon-usage of a set of DNA sequences (nucleotides) stored in a passed file
    => Absolute values / total codon abundances

    Args:
        dna_fasta_file (str): Absolute path to FASTA file of DNA sequences

    Returns:
        dict, dict: Weighted absolute codon usages (mono, di-codon usage)
    """
    # Initialize new local codon-lists
    codons, dicodons = initialize_codons()
    first_codon = None
    # Check file existence
    if os.path.isfile(dna_fasta_file):
        for seq_record in SeqIO.parse(dna_fasta_file, "fasta"):
            ppm_weight = int(float(seq_record.description.split(" ")[1]))
            step = 3  # 1 (all 3 frames) | 3 (specified "coding" frame)
            for n in range(0, len(seq_record)+1, step):
                first_codon = seq_record.seq[n:n+3]
                sec_codon = seq_record.seq[n+3:n+6]
                # Ignore internal occuring stop_codons
                if first_codon in all_codons and sec_codon in all_codons:
                    # Weighting by ppm
                    codons[str(first_codon)] += ppm_weight
                    dicodons[str(first_codon)][str(sec_codon)] += ppm_weight
            # Add last codon (stop-codon)
            if first_codon in all_codons:
                codons[str(first_codon)] += ppm_weight
    else:
        print("Can't find file at the specified location.")
    return codons, dicodons


def calculate_abs_codon_usages_string(dna_fasta_arr):
    """Calculate codon-usage of a set of DNA sequences (nucleotides) stored in passed string
    => Absolute values / total codon abundances

    Args:
        dna_fasta_arr (str): String with DNA sequences (FASTA formatted)

    Returns:
        dict, dict: Absolute codon usages (mono, di-codon usage)
    """
    # Initialize variables
    start_time = datetime.datetime.now()
    print('> Calc. codon usages ...')

    # Initialize new local codon-lists
    codons, dicodons = initialize_codons()
    first_codon = None
    # Count codon frequences
    if isinstance(dna_fasta_arr, list):
        for seq_record in dna_fasta_arr:
            step = 3  # 1 (all 3 frames) | 3 (specified "coding" frame)
            for n in range(0, len(seq_record)+1, step):
                first_codon = seq_record.seq[n:n+3]; sec_codon = seq_record.seq[n+3:n+6]
                # Ignore internal occuring stop_codons
                if first_codon in all_codons and sec_codon in all_codons:
                    codons[str(first_codon)] += 1
                    dicodons[str(first_codon)][str(sec_codon)] += 1
            # Add last codon (stop-codon)
            if first_codon in all_codons:
                codons[str(first_codon)] += 1
    else:
        print("The data has been passed in a wrong format.")

    # Measure time
    print('> Calc. time: ' + str(millis(start_time)/1000) + ' sec.')
    return codons, dicodons


def calculate_abs_codon_usages_weighted_string(dna_fasta_arr):
    """Calculate weighted codon-usage of a set of DNA sequences (nucleotides) stored in passed string
    => Absolute values / total codon abundances

    Args:
        dna_fasta_arr (str): String with DNA sequences (FASTA formatted)

    Returns:
        dict, dict: Weighted absolute codon usages (mono, di-codon usage)
    """
    # Initialize variables
    start_time = datetime.datetime.now()
    print('> Calc. codon usages ...')

    # Initialize new local codon-lists
    codons, dicodons = initialize_codons()
    first_codon = None
    # Count codon frequences
    if isinstance(dna_fasta_arr, list):
        for seq_record in dna_fasta_arr:
            ppm_weight = int(float(seq_record.description.split(" ")[1]))
            step = 3  # 1 (all 3 frames) | 3 (specified "coding" frame)
            for n in range(0, len(seq_record)+1, step):
                first_codon = seq_record.seq[n:n+3]
                sec_codon = seq_record.seq[n+3:n+6]
                # Ignore internal occuring stop_codons
                if first_codon in all_codons and sec_codon in all_codons:
                    # Weighting by ppm
                    codons[str(first_codon)] += ppm_weight
                    dicodons[str(first_codon)][str(sec_codon)] += ppm_weight
            # Add last codon (stop-codon)
            if first_codon in all_codons:
                codons[str(first_codon)] += ppm_weight
    else:
        print("The data has been passed in a wrong format.")

    # Measure time
    print('> Calc. time: ' + str(millis(start_time)/1000) + ' sec.')
    return codons, dicodons


def check_codons_rcu(codons):
    """Check mono codon usage (ACU or RCU) for sufficient data (Null)

    Args:
        codons (dict): Codon usage

    Returns:
        bool
    """
    # Check single codons
    for codon, amount in codons.iteritems():
        if amount == 0:
            print(codon + ' does not occur ...')
            return False


def check_codons(codons, dicodons):
    """Check mono- and di-codon usage (ACU or RCU) for sufficient data (Null)

    Args:
        codons (dict): Mono-Codon usage (ACU or RCU)
        dicodons (dict): Di-Codon usage (ACU or RCU)

    Returns:
        bool
    """
    # Check single codons
    for codon, amount in codons.iteritems():
        if amount == 0 and codon in valid_codons:
            print(codon + ' does not occur ...')
            return False

        # Select valid dicodon-combinations
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        subset_dicodons = dict((select, dicodons[codon][select]) for select in alt_codons)
        subset_len = len(subset_dicodons); cnt = 0
        # Check di-codons
        for dicodon, diamount in subset_dicodons.iteritems():
            if diamount == 0 and dicodon in valid_codons:
                # print(codon + dicodon + ' doesnt occur ...')
                cnt += 1
            # At least one dicodon opportunity needs to be greater Zero
            if subset_len == cnt:
                print('No valid dicodon-combinations found for: ' + codon)
                return False

    # Reach end, return True
    return True


def codons_rcu_fill_up(codons_rcu, ref_codons_rcu):
    """Fill up empty fields in mono codon usage (RCU) with data from reference MCU

    Args:
        codons_rcu (dict): Mono codon usage (RCU)
        ref_codons_rcu (dict): Reference mono codon usage (RCU)

    Returns:
        dict: Filled up codons (mono RCU)
    """
    for c, val in codons_rcu.iteritems():
        if float(val) == 0.0:
            codons_rcu[c] = ref_codons_rcu[c]
    # Return refilled RCU
    return codons_rcu


def dicodons_rcu_fill_up(dicodons_rcu, ref_dicodons_rcu):
    """Fill up empty fields in di-codon usage (RCU) with data from reference MCU

    Args:
        dicodons_rcu (dict): Di-codon usage (RCU)
        ref_dicodons_rcu (dict): Reference di-codon usage (RCU)

    Returns:
        dict: Filled up dicodons (DiRCU)
    """
    for start_c, codons in dicodons_rcu.iteritems():
        for end_c, val in codons.iteritems():
            if float(val) == 0.0:
                dicodons_rcu[start_c][end_c] = ref_dicodons_rcu[start_c][end_c]
    # Return refilled RCU
    return dicodons_rcu


def create_ordered_array_from_cu(codons):
    """Create codon-ordered array (Numpy) for further calculations

    Args:
        codons (dict): Mono codon usage (RCU)

    Returns:
        nparray: Ordered codons dict
    """
    codons_keys = sorted(codons.keys())
    codons_values = []
    for c in codons_keys:
        codons_values.append(codons[c])
    # Numpy array for further calculations
    return np.asarray(codons_values)


def calculate_rel_codon_usage(codons_acu):
    """Calculate relative global codon usage (RCU) of the pre-calculated absolute codon usage (total codon abundances)
    (normalized in total sum to 1.0)

    Args:
        codons_acu (dict): Mono codon usage (RCU)

    Returns:
        codons_rel: Relative  codon usage
    """
    # IMPORTANT: Assure that pairings are the same as before
    # Calculate Codon usage - relative values ()
    codons_values = []; codons_keys = sorted(codons_acu.keys())
    for c in codons_keys:
        codons_values.append(codons_acu[c])
    sums = sum(list(codons_values))
    freqs = np.divide(list(codons_values), float(sums))
    codons_rel = dict(zip(codons_keys, freqs))
    return codons_rel


def calculate_rel_prop_codon_usage(codons, codons_proportions):
    """Calculate relative global codon usage (RCU) of the pre-calculated absolute codon usage (total codon abundances)
    proportional to previous group distribution (normalized in total sum to 1.0)

    Args:
        codons (dict): Mono codon usage (RCU)
        codons_proportions:

    Returns:
        codonsPropRel: Relative proportional mono codon usage
    """
    # Calculate Codon usage - relative values ()
    codonsPropRel = {}
    # Calc. relative synonymous codon usage
    for codon, value in codons.iteritems():
        # Define all alternative codons for amino acid
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        # Build sum over all alternative codons of one group
        sum_alt_prop_codons = sum([codons_proportions[c] for c in alt_codons])
        # Calc. weighted/proportional values
        codonsPropRel[codon] = np.multiply(codons[codon], float(sum_alt_prop_codons))
    return codonsPropRel


def calculate_rel_syn_codon_usage(codons):
    """Calculate relative synonymous codon usage (RSCU) of the pre-calculated absolute codon usage
    [or relative codon usage (RCU)]
    Synonymous = Group-based; Normalization to equivalent/all alternative codons for an aminoacid

    Args:
        codons (dict): Mono codon usage (ACU or RCU)

    Returns:
        codonsPropRel: Relative synonymous mono codon usage
    """
    # Calculate Codon usage - relative values ()
    codonsRel = {}
    # Calc. relative synonymous codon usage
    for codon, value in codons.iteritems():
        # Define all altenative codons for amino acid
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        # Build sum over all alternative codons of one group
        sum_alt_codons = sum([codons[c] for c in alt_codons])
        # Calc. RSCU values
        # print("Codon: %s - Group sum: %0.2f" % (codon, float(sum_alt_codons)))
        codonsRel[codon] = np.divide(codons[codon], float(sum_alt_codons if sum_alt_codons > 0.0 else 1.0))
    return codonsRel


def calculate_rel_codon_usage_invert(codons, ref_codons):
    """Calculate inverted relative global codon usage (RCU) of the pre-calculated absolute codon usage(s)
    (total codon abundances, normalized in total sum to 1.0)

    Args:
        codons (dict): Mono codon usage (RCU)
        ref_codons (dict): Reference genome mono codon usage (RCU)

    Returns:
        codonsPropRel: Relative synonymous mono codon usage inverted
    """
    # Numpy: No scientific number representation
    np.set_printoptions(suppress=True)

    # Calculate Codon usage - relative values ()
    codons_keys = sorted(codons.keys())

    # Step 01: Translate abs. codon-usages into RCUs
    codonsRel = calculate_rel_codon_usage(codons)
    freqs = codonsRel.values()
    codonsRelRef = calculate_rel_codon_usage(ref_codons)
    ref_freqs = codonsRelRef.values()

    # View rel. codon-usages
    # codonsRel = dict(zip(codons_keys, freqs))
    print("# Target codons")
    print(codonsRel)
    print("sum: %f" % sum(freqs))
    print("min: %0.5f %s" % (min(freqs), codonsRel.keys()[codonsRel.values().index(min(freqs))]))
    print("max: %0.5f %s" % (max(freqs), codonsRel.keys()[codonsRel.values().index(max(freqs))]))
    # View rel. codon-usages (ref_codons)
    # codonsRelRef = dict(zip(codons_keys, ref_freqs))
    print("# Reference codons")
    print(codonsRelRef)
    print("sum: %f" % sum(ref_freqs))
    print("min: %0.5f %s" % (min(ref_freqs), codonsRelRef.keys()[codonsRelRef.values().index(min(ref_freqs))]))
    print("max: %0.5f %s" % (max(ref_freqs), codonsRelRef.keys()[codonsRelRef.values().index(max(ref_freqs))]))

    # Step 02: Calculate synonymous (group-based) codon-usages
    codonsRelSyn = calculate_rel_syn_codon_usage(codonsRel)
    print("# Target codons (syn.)")
    print(codonsRelSyn)
    print("sum: %f" % sum(codonsRelSyn.values()))
    codonsRelRefSyn = calculate_rel_syn_codon_usage(codonsRelRef)
    print("# Reference codons (syn.)")
    print(codonsRelRefSyn)
    print("sum: %f" % sum(codonsRelRefSyn.values()))
    freqs_syn = []; ref_freqs_syn = []
    # Garantuee same order of keys and values
    for c in codons_keys:
        freqs_syn.append(codonsRelSyn[c])
        ref_freqs_syn.append(codonsRelRefSyn[c])

    # Step 03: Mirror CUs
    ref_freqs_syn = np.asarray(ref_freqs_syn)
    freqs_syn = np.asarray(freqs_syn)
    weight_lambda = 1.40
    freqs_invert_syn = ref_freqs_syn*np.exp(weight_lambda * np.log(ref_freqs_syn/freqs_syn))
    # Collect inf/nan
    freqs_invert_syn[np.isnan(freqs_invert_syn)] = 0.0
    print("# Diff freqs")
    print(codons_keys)
    # print(diff_syn)
    print("# Inverted freqs")
    print(freqs_invert_syn)

    # Normalize to 1.0
    freqsInvertDict = dict(zip(codons_keys, freqs_invert_syn))
    freqsInvertSynDict = calculate_rel_syn_codon_usage(freqsInvertDict)

    print("# Normed RSCU dict with inverted freqs")
    print(freqsInvertSynDict)
    print("sum: %f" % sum(freqsInvertSynDict.values()))
    print("Min: %f" % min(freqsInvertSynDict.values()))
    print("Max: %f" % max(freqsInvertSynDict.values()))

    # Normalize with group proportions
    codonsRelInvertDict = calculate_rel_prop_codon_usage(freqsInvertSynDict, codonsRel)
    print("# Normed Proportional RCU with inverted freqs")
    print(codonsRelInvertDict)
    print("sum: %f" % sum(codonsRelInvertDict.values()))

    return codonsRelInvertDict


def convert_RCU_to_RSCU(codons_rel):
    """Convert relative (global) codon usage (RCU) to relative synonymous codon usage (RSCU)

    Args:
        codons_rel (dict): Relative mono codon usage (RCU)

    Returns:
        codonsRSCU: Relative synonymous codon usage (RSCU)
    """
    # Calculate Codon usage - relative synonymous values ()
    codonsRSCU = {}
    # Calc. relative synonymous codon usage
    for codon, value in codons_rel.iteritems():
        # Define all altenative codons for amino acid
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        # Build sum over all alternative codons of one group
        sum_alt_codons = sum([codons_rel[c] for c in alt_codons])
        # Calc. RSCU values
        codonsRSCU[codon] = np.divide(codons_rel[codon], float(sum_alt_codons if sum_alt_codons > 0.0 else 1.0))
    return codonsRSCU


def calculate_rel_dicodon_usage(dicodons):
    """Calculate relative di-codon-usage (RCU) of the pre-calculated absolute di-codon-usage (total di-codon abundances)
    (normalized in total to 1.0)

    Args:
        dicodons (dict): Dictionary of the absolute di-codon-usages

    Returns:
        dicodonsRel: Relative di-codon usage (RCU)
    """
    # IMPORTANT: Assure that pairings are the same as before
    # Calculate Codon usage - relative values ()
    codons_keys = sorted(dicodons.keys())
    dicodonsRel = dicodons.copy()
    # Calculate Di-Codon usage - relative values (normalized to each line = first codon fixed)
    # Total sum over all codons (normalized to whole matrix)
    sums = sum(freq for fc, freqs in dicodons.items() for sc, freq in freqs.items())
    for first_codon, sec_codons in dicodons.items():
        sec_codons_values = []
        for c in codons_keys:
            sec_codons_values.append(sec_codons[c])
        freqs = np.divide(list(sec_codons_values), float(sums if sums > 0.0 else 1.0))
        dicodonsRel[first_codon] = dict(zip(codons_keys, freqs))
    return dicodonsRel


def calculate_rel_dicodon_usage_based_rel_codon_usage(dicodons, codons):
    """Calculate relative di-codon-usage (RCU) of the pre-calculated absolute di-codon-usage (total di-codon abundances)
    weighted by mono-codon-usage (RCU) (normalized in total to 1.0)

    Args:
        dicodons (dict): Dictionary of the absolute di-codon-usages
        codons (dict): Dictionary of the absolute mono-codon-usages

    Returns:
        dicodonsRel: Relative di-codon usage (RCU)
    """
    # IMPORTANT: Assure that pairings are the same as before
    # Calculate Codon usage - relative values ()
    codons_keys = sorted(dicodons.keys())
    dicodonsRel = dicodons.copy()
    # Calculate Di-Codon usage - relative values (normalized to each line = first codon fixed)
    # Total sum over all codons (normalized to whole matrix)
    sums = sum(freq for fc, freqs in dicodons.items() for sc, freq in freqs.items())
    for first_codon, sec_codons in dicodons.items():
        sec_codons_values = []
        for c in codons_keys:
            sec_codons_values.append(sec_codons[c])
        freqs = np.divide(list(sec_codons_values), float(codons[first_codon] if codons[first_codon] > 0.0 else 1.0))
        dicodonsRel[first_codon] = dict(zip(codons_keys, freqs))
    return dicodonsRel


def calculate_rel_dicodon_usage_invert(dicodons_acu, ref_dicodons_acu, codons_acu):
    """Calculate inverted relative di-codon-usage (RCU) of the pre-calculated absolute di-codon-usage and reference
    ACU (total di-codon abundances, normalized in total to 1.0)

    Args:
        dicodons_acu (dict): Dictionary of the absolute di-codon-usages (dataset to invert)
        ref_dicodons_acu (dict): Dictionary of the absolute di-codon-usages (reference-dataset to mirror by)
        codons_acu (dict): Dictionary of the absolute di-codon-usages

    Returns:
        dicodonsRelInvertRCU: Inverted relative di-codon usage (RCU)
    """

    # Define result dictionary
    dicodonsRelInvert = {}
    codons_keys = sorted(dicodons_acu.keys())
    # Convert ACU and ACU_Ref into RCUs
    dicodonsRel = calculate_rel_dicodon_usage(dicodons_acu)
    dicodonsRelRef = calculate_rel_dicodon_usage(ref_dicodons_acu)

    # Invert columns of DiCodonMatrix (RCU) by traversing row-by-row
    for first_codon, sec_codons in dicodonsRel.items():

        # print(first_codon)
        # Step 01 : Calculate synonymous (group-based) codon-usages
        codonsRelSyn = calculate_rel_syn_codon_usage(dicodonsRel[first_codon])
        print("sum: %f" % sum(codonsRelSyn.values()))
        codonsRelRefSyn = calculate_rel_syn_codon_usage(dicodonsRelRef[first_codon])
        print("sum: %f" % sum(codonsRelRefSyn.values()))

        freqs_syn = []; ref_freqs_syn = []
        # Garantuee same order of keys and values
        for c in codons_keys:
            freqs_syn.append(codonsRelSyn[c])
            ref_freqs_syn.append(codonsRelRefSyn[c])

        # Step 02: Mirror CUs
        ref_freqs_syn = np.asarray(ref_freqs_syn)
        freqs_syn = np.asarray(freqs_syn)
        diff_syn = ref_freqs_syn - freqs_syn

        # Method 01: Simple difference # freqs_invert = freqs - ( 2*(freqs-ref_freqs) )
        freqs_invert_syn = ref_freqs_syn + diff_syn
        # Correct values less than 0 by splitting the negative part to corresponding codons
        min_val = 0.0001
        for i, freq in enumerate(freqs_invert_syn):
            if freq < 0:
                codon = codons_keys[i]
                alt_codons = list(codonMapping['backward'][codonMapping['forward'][codon]])
                alt_codons.pop(alt_codons.index(codon))
                # Correct group freqs weighted
                alt_freqs = {}
                for c in alt_codons:
                    c_i = codons_keys.index(c)
                    alt_freqs.__setitem__(c, freqs_invert_syn[c_i])
                alt_freqs_rcu = calculate_rel_codon_usage(alt_freqs)
                # Correct all alternative codons of one group
                for c in alt_codons:
                    c_i = codons_keys.index(c)
                    correct_val = abs(freq - min_val) * alt_freqs_rcu[c]
                    if freqs_invert_syn[c_i] < ref_freqs_syn[c_i]:
                        freqs_invert_syn[c_i] += correct_val
                    else:
                        freqs_invert_syn[c_i] -= correct_val
                # Set initial codon-freq to minimum
                freqs_invert_syn[i] = min_val

                # => AA-Groups (inverted) need to sum again up to 1.0
                # TODO: Check

        # Merge into dict.
        freqsInvertSynDict = dict(zip(codons_keys, freqs_invert_syn))

        # Step 03: Normalize rows to RCU with previous group proportions
        freqsInvertDict = calculate_rel_prop_codon_usage(freqsInvertSynDict, sec_codons)
        print(freqsInvertDict)
        dicodonsRelInvert[first_codon] = freqsInvertDict

    # Store
    dicodonsRelInvertRCU = dicodonsRelInvert

    # Step 04:  Normalize to RCU with previous group proportions
    print('Step 04: Before Conversion to RCU')
    sums_01 = sum(freq for fc, freqs in dicodonsRelInvert.items() for sc, freq in freqs.items())
    print(sums_01)

    codons, dicodonsRelDoubleInvert = initialize_codons()
    # Invert rows of DiCodonMatrix (RCU) by traversing col-by-col
    for first_codon, sec_codons in dicodonsRelInvert.items():

        # Transpose col-entries to row
        codonsRel = {}; codonsRelRef = {}; codonsRelAGroup = {}; codonsRelRefAGroup = {};
        for c in codons_keys:
            # Boxes of occuring aminoacids
            codonsRel[c] = dicodonsRelInvert[c][first_codon]
            codonsRelRef[c] = dicodonsRelRef[c][first_codon]

        # Step 01 : Calculate synonymous (group-based) codon-usages
        codonsRelSyn = calculate_rel_syn_codon_usage(codonsRel)
        codonsRelRefSyn = calculate_rel_syn_codon_usage(codonsRelRef)

        freqs_syn = []; ref_freqs_syn = []
        # Garantuee same order of keys and values
        for c in codons_keys:
            freqs_syn.append(codonsRelSyn[c])
            ref_freqs_syn.append(codonsRelRefSyn[c])

        # Step 02: Mirror CUs
        ref_freqs_syn = np.asarray(ref_freqs_syn)
        freqs_syn = np.asarray(freqs_syn)
        diff_syn = ref_freqs_syn - freqs_syn

        # Method 01: Simple difference # freqs_invert = freqs - ( 2*(freqs-ref_freqs) )
        freqs_invert_syn = ref_freqs_syn + diff_syn
        # Correct values less than 0 by splitting the negative part to corresponding codons
        min_val = 0.0001
        for i, freq in enumerate(freqs_invert_syn):
            if freq < 0:
                codon = codons_keys[i]
                alt_codons = list(codonMapping['backward'][codonMapping['forward'][codon]])
                alt_codons.pop(alt_codons.index(codon))
                # Correct group freqs weighted
                alt_freqs = {}
                for c in alt_codons:
                    c_i = codons_keys.index(c)
                    alt_freqs.__setitem__(c, freqs_invert_syn[c_i])
                alt_freqs_rcu = calculate_rel_codon_usage(alt_freqs)
                # Correct all alternative codons of one group
                for c in alt_codons:
                    c_i = codons_keys.index(c)
                    correct_val = abs(freq - min_val) * alt_freqs_rcu[c]
                    if freqs_invert_syn[c_i] < ref_freqs_syn[c_i]:
                        freqs_invert_syn[c_i] += correct_val
                    else:
                        freqs_invert_syn[c_i] -= correct_val
                # Set initial codon-freq to minimum
                freqs_invert_syn[i] = min_val

        # Merge into dict.
        freqsInvertSynDict = dict(zip(codons_keys, freqs_invert_syn))

        # Step 03: Normalize rows to RCU with previous group proportions
        freqsInvertDict = calculate_rel_prop_codon_usage(freqsInvertSynDict, codonsRel)

        # Transpose back
        for c in codons_keys:
            dicodonsRelDoubleInvert[c][first_codon] = freqsInvertDict[c]

    # Store
    dicodonsRelInvertRCU = dicodonsRelDoubleInvert

    print('Step 04: After Conversion to RCU')
    sums_02 = sum(freq for fc, freqs in dicodonsRelInvertRCU.items() for sc, freq in freqs.items())
    print(sums_02)

    return dicodonsRelInvertRCU


# Method: Calculate/Convert relative Di-Codon usage of the pre-calculated
#         absolute Di-Codon usage
#
# Parameter:
# dicodons[dict]   Dictionary of the absolute di-codon-usages
#
def calculate_groupwise_rel_dicodon_usage_invert(dicodons_acu, ref_dicodons_acu):
    """Calculate groupwise inverted relative di-codon-usage (RCU) of the pre-calculated absolute di-codon-usage and
    reference ACU (total di-codon abundances, normalized in total to 1.0)

    Args:
        dicodons_acu (dict): Dictionary of the absolute di-codon-usages (dataset to invert)
        ref_dicodons_acu (dict): Dictionary of the absolute di-codon-usages (reference-dataset to mirror by)

    Returns:
        dicodonsRelInvertRCU: Inverted relative di-codon usage (RCU)
    """
    dicodonsRelInvert = {}  # Result object
    dicodons_rcu = calculate_rel_dicodon_usage(dicodons_acu)
    ref_dicodons_rcu = calculate_rel_dicodon_usage(dicodons_acu)

    # Calculate Di-Codon usage - relative values (normalized to each line = first codon fixed)
    for first_codon, sec_codons in dicodons_rcu.items():
        # Calculate Codon usage - relative values ()

        # Step 01 : Calculate synonymous (group-based) codon-usages
        codonsRelSyn = calculate_rel_syn_codon_usage(dicodons_rcu[first_codon])
        codonsRelRefSyn = calculate_rel_syn_codon_usage(ref_dicodons_rcu[first_codon])

        # Step 02: Mirror CUs
        # Invert groupwise by flipping values in AA-groups
        codonsRelSynInvert = deepcopy(codonsRelSyn)
        for aa, syn_codons in codonMapping['backward'].iteritems():
            # Collect group values
            ref_group_rcu = dict((id, codonsRelRefSyn[id]) for id in syn_codons)
            ref_avg = mean(ref_group_rcu.values())
            for id in syn_codons:
                codonsRelSynInvert[id] = codonsRelSyn[id] - 2*(codonsRelSyn[id] - ref_avg)

        # Step 04:  Normalize to RCU with previous group proportions
        codonsRelInvertDict = calculate_rel_prop_codon_usage(codonsRelSynInvert, dicodons_rcu[first_codon])
        dicodonsRelInvert[first_codon] = codonsRelInvertDict

    # Step 04:  Normalize to RCU with previous group proportions
    codons_rcu = {}
    for first_codon, sec_codons in dicodons_rcu.iteritems():
        codons_rcu[first_codon] = sum(sec_codons.values())

    dicodonsRelInvertRCU = dicodonsRelInvert

    return dicodonsRelInvertRCU


def generate_rscu(codons):
    """Calculate relative synonymous codon usage (RSCU) of the pre-calculated absolute codon usage (total codon abundances)

    Implemented after: Li et. all (1987) The codon adaptation index - a measure of directional synonymous codon usage
    bias, and its potential applications. NAR. 15(3)

    Args:
        codons (dict): Dictionary of the absolute mono codon-usages

    Returns:
        dict, dict: RSCU, relative Adaptiveness
    """

    codons_abs = codons.copy()
    rscu = codons.copy()
    rel_adapt = codons.copy()

    # Calc. relative synonymous codon usage
    for codon, value in codons_abs.iteritems():
        # Define all altenative codons for amino acid
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        n = len(alt_codons)
        # Build sum over all alternative codons of one group
        sum_alt_codons = sum([codons_abs[c] for c in alt_codons])
        # Calc. RSCU values
        rscu[codon] = np.divide(codons_abs[codon], sum_alt_codons/float(n))

    # Calc. relative Adaptiveness
    for codon in rscu:
        # Define all alternative codons for amino acid
        alt_codons = codonMapping['backward'][codonMapping['forward'][codon]]
        max_rscu_codon = max([rscu[c] for c in alt_codons])
        # Calc. Relative Adaptiveness
        rel_adapt[codon] = np.divide(rscu[codon], max_rscu_codon)
        # print('{:s} {:03.3f} {:03.3f} {:d} {:d}'.format(codon, rscu[codon], rel_adapt[codon], n, sum_alt_codons))

    return rscu, rel_adapt


def calculate_cai(sequence, rel_adapt):
    """Calculate codon adaptation index (CAI) for a gene as the geometric mean of the RSCU values corresponding to
    each of the codons used in the gene

    Implemented after: Li et. all (1987) The codon adaptation index - a measure of directional synonymous codon usage
    bias, and its potential applications. NAR. 15(3)

    Args:
        sequence (str): Sequence
        rel_adapt (dict): relative Adaptiveness

    Returns:
        float: CAI
    """
    adapt_vals = []; step = 3; CAI = 0.0
    for n in range(0, len(sequence), step):
        codon = sequence[n:n+step]
        if codon not in stop_codons:
            adapt_vals.append(np.log(rel_adapt[codon]))
            CAI = np.exp(np.sum(adapt_vals)/(len(sequence)/step))
        else:
            print("Sequence contains internal stop codons. ", codon)

    return CAI


def calculate_major_codons(codons):
    """Calculate the major codons from the absolute or relative codon usage

    Args:
        codons (dict): ACU or RCU

    Returns:
        dict: Major codons
    """
    major_codons = {}
    # Calculate Major Codons
    for aa, syn_codons in codonMapping['backward'].iteritems():
        # Print all values for AA
        # print [(codons[c], c) for c in syn_codons]
        # Show Major codon for AA
        print (aa, max([(codons[c], c) for c in syn_codons]))
        major_codons[aa] = max([codons[c] for c in syn_codons])

    return major_codons


def calculate_codon_usage_similarity(cu_one, cu_two):
    """Calculate average similarity between two relative codon-usages. Assumes that both dicts are normalized and
    structurally the same / even

    Args:
        cu_one (dict): RCU
        cu_two (dict): RCU

    Returns:
        float: Similarity between 0.0 and 1.0 (1.0 == same)
    """
    similarity = 0.0; differences = {}; sum_up = 0.0
    codons = sorted(cu_one.keys())
    # Check structure of both CUs (length, sorting, same fields)
    if codons == sorted(cu_two.keys()):
        # Calculate similarity
        for codon in codons:
            differences[codon] = abs(cu_one[codon]-cu_two[codon])
            sum_up += abs(cu_one[codon]+cu_two[codon])

        similarity = 1.0 - np.divide(sum([diff for diff in differences.values()]), sum_up)
        return np.round(similarity, 16)
    else:
        print "The passed arguments are structural not even."


def remove_stop_codons(seq):
    """Remove all stop codons (ending)

    Args:
        seq (str):

    Returns:
        str: Modified sequence
    """
    for codon in stop_codons:
        seq = re.sub(codon + "$", "", seq)
    return seq


def distance_from_model_entropy(setup, model_entropy, bioObj):
    """Calculate distance between sample and model entropy

    Args:
        setup (dict): Current setup (mono-, di-codon-usage, Markov-chain)
        model_entropy: Value of model_entropy
        bioObj: FASTA-SeqRecord (BioPython) of reference sequence

    Returns:
        float: Distance
    """
    # Generate random sequences (strings)
    seq_info = generate_probability_enriched_sequence(bioObj, setup['G'], setup['nl2id'])
    seq_distance = calculate_entropy_distance(seq_info.seq_info, model_entropy, setup['G'], setup['nl2id'])

    return abs(seq_distance)


def mean_sample_distance(setup, sample_size):
    """Calculate mean distance from sample to model entropy (sampled average distance measure)

    Args:
        setup (dict): Current setup (mono-, di-codon-usage, Markov-chain)
        sample_size: Number of samples to generate

    Returns:
        float: Mean Distance
    """
    # Calc entropy and distance measure
    entropy = calculate_model_entropy(setup['G'], setup['nl2id'])

    # Generate random sequences (strings)
    sample_distances = []; # samples = [];
    for i in range(1, sample_size):
        length = random.randint(141, 141)
        seq_info = generate_random_sequence(length, setup['G'], setup['nl2id'], setup['nl2codon'], extended=True)
        seq_distance = calculate_entropy_distance(seq_info, entropy, setup['G'], setup['nl2id'])
        sample_distances.append(seq_distance)

    print("- Min. Entropy-distance: %8.15f" % (min(sample_distances)))
    print("- Max. Entropy-distance: %8.15f" % (max(sample_distances)))
    print("- Mean Entropy-distance: %8.15f" % (mean(sample_distances)))
    print("- Std. Entropy-distance: %8.15f" % (np.std(sample_distances)))
    print("- Var. Entropy-distance: %8.15f" % (np.var(sample_distances)))

    return mean(sample_distances)


def mean_sample_entropy(setup, sample_size):
    """Calculate mean sample entropy (sampled average entropy)

    Args:
        setup (dict): Current setup (mono-, di-codon-usage, Markov-chain)
        sample_size: Number of samples to generate

    Returns:
        float: Mean sample entropy
    """
    # Generate random sequences (strings)
    entropies = []
    for i in range(1, sample_size):
        length = random.randint(141, 141)
        seq_info = generate_random_sequence(length, setup['G'], setup['nl2id'], setup['nl2codon'], extended=True)
        sample_entropy = calculate_sample_entropy(seq_info, setup['G'], setup['nl2id'])
        entropies.append(sample_entropy)

    print("- Min. Entropy: %8.15f" % (min(entropies)))
    print("- Max. Entropy: %8.15f" % (max(entropies)))
    print("- Mean Entropy: %8.15f" % (mean(entropies)))
    print("- Std. Entropy: %8.15f" % (np.std(entropies)))
    print("- Var. Entropy: %8.15f" % (np.var(entropies)))

    return mean(entropies)


########################################################################################################################


def identify_restriction_sites(artificial_seqs, restriction_sites):
    """Search and identify passed restrictions-sites in artificially created heterologous sequences

    Args:
        artificial_seqs (list): List of artificial sequences
        restriction_sites (dict): Parsed/Structured REBASE-Bairoch dataset

    Returns:
        dict: Match object with all matching restriction sites
    """
    overall_matches = {'prototypes': [], 'matches': {}}; i = -1
    for record in artificial_seqs:
        local_matches = {'prototypes': [], 'matches': {}}; i += 1
        # Iterate over Single restriction sites
        for key, enzym in restriction_sites.iteritems():
            for site_pattern in enzym['RSS']:
                # Filter restriction sites with size larger than x[=3]
                #    and by "unused" Prototype-restriction sequences
                if len(site_pattern) >= 4:
                    local_matches['prototypes'].append(enzym['PT'])
                    for m in re.finditer(site_pattern, str(record.seq)):
                        # print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0)))
                        try:
                            local_matches['matches'][site_pattern].append([i, m.start(), m.end()])
                        except KeyError:
                            local_matches['matches'][site_pattern] = [[i, m.start(), m.end()]]
                        try:
                            overall_matches['matches'][site_pattern].append([i, m.start(), m.end()])
                        except KeyError:
                            overall_matches['matches'][site_pattern] = [[i, m.start(), m.end()]]

        # Merge results
        overall_matches['prototypes'] = list(set(overall_matches['prototypes'] + local_matches['prototypes']))

    return overall_matches


def check_restriction_sites(artificial_seqs, restriction_sites):
    """Check for passed restrictions-sites in artificially created heterologous sequences

    Args:
        artificial_seqs (list): List of artificial sequences
        restriction_sites (dict): Parsed/Structured REBASE-Bairoch dataset

    Returns:
        bool
    """
    for record in artificial_seqs:
        # Iterate over Single restriction sites
        for key, enzym in restriction_sites.iteritems():
            for site_pattern in enzym['RSS']:
                # Filter restriction sites with size larger than x[=3]
                #    and by "unused" Prototype-restriction sequences
                if len(site_pattern) >= 4:
                    if re.findall(site_pattern, str(record.seq)):
                        return True
    return False


def parse_bairoch_file(input_file):
    """Parse REBASE Bairoch files

    Args:
        input_file (str): Path to REBASE (Bairoch) file

    Returns:
        list: List of parsed REBASE restriction site entries
    """
    bairochArr = []

    if os.path.isfile(input_file):
        # Reading content from specified file
        f = open(input_file, 'r')
        fileContent = f.read()            # Read in one String
        fileArr = fileContent.split("//") # Split in lines

        # Simple grammar to match
        ident = oneOf("ID ET OS PT RS MS CR CM RN RA RL", caseless=False)
        macroDef = ident.setResultsName("name") + pp.empty + restOfLine.setResultsName("value")
        for dataset in fileArr:
            tmp = {}
            for t, s, e in macroDef.scanString(dataset):
                # print(t.name,":", t.value)
                # if ["RN", "RA", "RL"].index(t.name):
                tmp[t.name] = t.value
            if bool(tmp):
                bairochArr.append(tmp)
    else:
        print("Can't find file at the specified location.")
        print(input_file)
    return bairochArr


def bairoch_select_restriction_enzymes(bairochArr, selected):
    """Filter selected enzymes of all available Bairoch restriciton sites

    Args:
        bairochArr (list): List of dicts with parsed REBASE file (Bairoch) information
        selected (list): List of enzyms (IDs)

    Returns:
        dict: Filtered dict with parsed REBASE file (Bairoch) information
    """
    bairochAvailEnzymes = {}; restrictionList = []
    for dataset in bairochArr:
        try:
            # Filter by number of providers
            if dataset['ID'] in selected:
                # Parse RestrictionSites -  Filter empty entries
                restrictionSites = [x.strip().split(',')[0] for x in dataset['RS'].split(";") if x]
                restrictionList.append(restrictionSites)
                dataset["RSS"] = restrictionSites
                bairochAvailEnzymes[dataset['ID']] = dataset

        except KeyError:
            next

    print('Number of filtered restriction sites: ', len(bairochAvailEnzymes))
    print('Restriction sites to look for: ', bairochAvailEnzymes)
    return bairochAvailEnzymes


def bairoch_filter_identical_enzymes(bairochArr, commercial=False):
    """Filter functional identical available enzymes

    Args:
        bairochArr (list): List of dicts with parsed REBASE file (Bairoch) information
        commercial (bool): Additional filter for commercial enzymes

    Returns:
        dict: Filtered list of dicts with parsed REBASE file (Bairoch) information
    """
    bairochAvailEnzymes = {}; restrictionList = []
    for dataset in bairochArr:
        try:
            if commercial:
                # Filter by number of providers
                if dataset['CR'] != '.' and len(dataset['CR']) > 1:
                    # Parse RestrictionSites -  Filter empty entries
                    restrictionSites = [x.strip().split(',')[0] for x in dataset['RS'].split(";") if x]
                    restrictionList.append(restrictionSites)
                    dataset["RSS"] = restrictionSites
                    try:
                        bairochAvailEnzymes[dataset['PT']].append(dataset)
                    except KeyError:
                        bairochAvailEnzymes[dataset['PT']] = []
                        bairochAvailEnzymes[dataset['PT']].append(dataset)
            else:
                # Parse RestrictionSites -  Filter empty entries
                restrictionSites = [x.strip().split(',')[0] for x in dataset['RS'].split(";") if x]
                restrictionList.append(restrictionSites)
                dataset["RSS"] = restrictionSites
                try:
                    bairochAvailEnzymes[dataset['PT']].append(dataset)
                except KeyError:
                    bairochAvailEnzymes[dataset['PT']] = []
                    bairochAvailEnzymes[dataset['PT']].append(dataset)

        except KeyError:
            next

    return bairochAvailEnzymes


def bairoch_filter_commercial_enzymes(bairochArr):
    """Filter commercially available enzymes

    Args:
        bairochArr (list): List of dicts with parsed REBASE file (Bairoch) information

    Returns:
        dict: Filtered list of dicts with parsed REBASE file (Bairoch) information
    """
    bairochAvailEnzymes = []; restrictionList = []
    for dataset in bairochArr:
        try:
            # Filter by number of providers
            if dataset['CR'] != '.' and len(dataset['CR']) > 1:
                # Parse RestrictionSites -  Filter empty entries
                restrictionSites = [x.strip().split(',')[0] for x in dataset['RS'].split(";") if x]
                restrictionList.append(restrictionSites)
                dataset["RSS"] = restrictionSites
                bairochAvailEnzymes.append(dataset)
        except KeyError:
            next
    print('Number of filtered restriction sites: ', len(bairochAvailEnzymes))
    return bairochAvailEnzymes


########################################################################################################################


def parse_paxdb(input_file):
    """Parse PAX-DB files

    Args:
        input_file (str): Filename of the PAX-DB dataset of specific species

    Returns:
        dict: Parsed (grammar) PAX-DB file information
    """
    # Read file-content in one string
    f = open(input_file, 'r')
    fileContent = f.read()

    # Commentary parser
    key = Word(alphas)
    commentary = Literal("#") + key.setResultsName("key") + ":" + restOfLine.setResultsName("comment")
    for t, s, e in commentary.scanString( fileContent ):
        print(t.key, ":", t.comment)

    # Simple grammar to extract information in hash
    local_ident = Word(nums)
    ncbi_ident  = Word(alphanums + ".-")
    macroDef = pp.LineStart() + local_ident.setResultsName("id") + pp.empty + ncbi_ident.setResultsName("ncbi") + pp.empty + restOfLine.setResultsName("abundance")

    # Execute grammar
    i = 0; paxDbArr = {}
    for t, s, e in macroDef.scanString( fileContent ):
        # print(t.id, ":", t.ncbi.split('.')[1], ":", t.abundance)
        paxDbArr[str(t.ncbi.split('.')[1])] = [t.abundance, t.ncbi.split('.')[0]]

    return paxDbArr


def pull_ncbi_dna_data(paxDbArr, organism, outfile, begin):
    """Download the referenced nucleotid PAX-DB files

    Args:
        paxDbArr (dict): Dictionary parsed PAX-DB files information
        organism (str): Currently used species (relevant for specification of NCBI requests)
        outfile (str): Output filename with project path
        begin (int): Start index for downloading records

    Returns:
        list: List of parsed FASTA-BioSeq-objects
    """
    Entrez.email = "dominic.simm@mpibpc.mpg.de"     # Always tell NCBI who you are
    nucleotides = []; success = []; failure = []; step = 10; batch_size = 5
    start = int(begin/step) if begin > 0 else 0  # 0 | int(3450/step)
    for times in range(start, int(len(paxDbArr)/step)):
        print("Main: Going to download record %i to %i" % (times*step, (times+1)*step))

        search_term = "(" + "[Gene Name] OR ".join(list(paxDbArr.keys())[times*step:(times+1)*step]) + "[Gene Name])"
        print(search_term + " AND " + organism + "[Organism]")

        try:
            search_handle = Entrez.esearch(db="gene", term=search_term + " AND " + organism + "[Organism]", usehistory="y")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            print(search_results)

            count = int(search_results["Count"])
            print("Found %i results" % count)

            # Store session information
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Begin downloading
            for start in range(0, count, batch_size):
                end = min(count, start+batch_size)
                print("- Sub: Going to download record %i to %i" % (start+1, end))
                fetch_handle = Entrez.efetch(db="gene", rettype="gb", retmode="xml",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
                data = Entrez.read(fetch_handle)
                fetch_handle.close()
                # print(data)

                # Get relevant gene information
                for i in range(0, len(data)):
                    ncbi_gene_id = ''; ncbi_ref_genomic_id = ''; ncbi_ref_genomic_acc = ''; ncbi_ref_mrna_acc = '';
                    ncbi_ref_protein_acc = ''; paxdb_locus_tag = '';
                    exception = False;
                    print(data[i])
                    # Get the information
                    try:
                        # Requested origin_id should be same as in the first request!
                        ncbi_gene_id = data[i]['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
                        ncbi_ref_genomic_id = data[i]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-int']
                        ncbi_ref_genomic_acc = data[i]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
                        paxdb_locus_tag = data[i]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str2']
                        ncbi_ref_mrna_acc = data[i]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_accession']
                        version = data[i]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_version']
                    except KeyError as e:
                        print("KeyError", sys.exc_info())
                        exception = True
                        next  # Skip this entry and jump to next one

                    # Fetch the file - DNA sequence
                    if ncbi_ref_mrna_acc.startswith('NM_') and exception == False:
                        record = None
                        try:
                            print("%s.%s" % (ncbi_ref_mrna_acc, version))
                            seq_handle = Entrez.efetch(db="nucleotide", id="%s.%s"%(ncbi_ref_mrna_acc,version), rettype="fasta", retmode="text")
                            record = SeqIO.read(seq_handle, "fasta")
                            # Customize header
                            record.id += "|" + ncbi_gene_id + "|" + paxdb_locus_tag
                            record.name += "#"
                            nucleotides.append(record)
                            success.append([ncbi_gene_id, ncbi_ref_genomic_id, ncbi_ref_genomic_acc, paxdb_locus_tag, ncbi_ref_mrna_acc, version])

                        except:
                            print("Unexpected error:", sys.exc_info())
                            failure.append([ncbi_gene_id, ncbi_ref_genomic_id, ncbi_ref_genomic_acc, paxdb_locus_tag, ncbi_ref_mrna_acc, version])
                            # failure.append(data[i])

                    else:
                        print("No mRNA found or unexpected error ...")
                        failure.append([ncbi_gene_id, ncbi_ref_genomic_id, ncbi_ref_genomic_acc, paxdb_locus_tag, ncbi_ref_mrna_acc, version])
                        # failure.append(data[i])
        except:
            print("Unexpected error:", sys.exc_info())
            return nucleotides

    # Write downloaded content in a batch
    SeqIO.write(nucleotides, outfile, "fasta")
    with open('%s/%s.failure.csv' % (os.path.dirname(outfile), os.path.basename(outfile)), 'wb') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerows(failure)
    with open('%s/%s.success.csv' % (os.path.dirname(outfile), os.path.basename(outfile)), 'wb') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerows(success)

    return nucleotides, failure


def pull_ncbi_protein_data(paxDbArr, organism, organism_id, begin):
    """Download the referenced protein PAX-DB files

    Args:
        paxDbArr (dict): Dictionary parsed PAX-DB files information
        organism (str): Currently used species (relevant for specification of NCBI requests)
        organism_id (int): Taxonomy identifier of currently used species (relevant for specification of NCBI requests)
        begin (int): Start index for downloading records

    Returns:
        list: List of parsed FASTA-BioSeq-objects
    """
    Entrez.email = "dominic.simm@mpibpc.mpg.de"     # Always tell NCBI who you are
    proteins = {}; step = 2; batch_size = 3
    for times in range(0, int(len(paxDbArr)/step)):
        print("\n" + "Going to download PaxDB record %i to %i" % (times*step, (times+1)*step))

        search_arr = list(paxDbArr.keys())[times*step:(times+1)*step]
        search_term_fp = "(" + ("-" + str(organism_id) + " OR ").join(search_arr) + "-" + str(organism_id) + ")"
        search_term = search_term_fp + " AND " + organism + "[Organism]"
        print(search_term)

        # Request
        search_handle = Entrez.esearch(db="protein", term=search_term, usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()

        count = int(search_results["Count"])
        print("Found %i results" % count)

        # Store session information
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Begin downloading
        for start in range(0, count, batch_size):
            end = min(count, start+batch_size)
            print("Check and download requested record %i to %i" % (start+1, end))
            fetch_handle = Entrez.efetch(db="protein", rettype="gb", retmode="xml",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
            dataset = Entrez.read(fetch_handle) # Read and parse XML response
            for data in dataset:
                paxDb_id = None
                for feat in data['GBSeq_feature-table']:
                    if feat['GBFeature_key'] == 'gene':
                        for qual in feat['GBFeature_quals']:
                            if qual['GBQualifier_name'] == 'locus_tag':
                                print(qual['GBQualifier_value']),
                                paxDb_id = qual['GBQualifier_value']
                try:
                    if paxDb_id in search_arr:
                        proteins[paxDb_id] = {'primary-accession': data['GBSeq_primary-accession'], 'sequence': data['GBSeq_sequence'].upper(), 'keywords': data['GBSeq_keywords']}
                        print(" ......... ok"),
                    else:
                        print(" ......... fail"),
                except KeyError:
                    print("Keyerror:", sys.exc_info())
            fetch_handle.close()
    # Return
    print(proteins)
    return proteins


def create_paxdb_compatible_dataset(paxDbArr, input_file):
    """Prepare parsed PAX-DB information to store it easily into the database

    Args:
        paxDbArr (dict): Parsed Nucleotide information (PAX-DB/NCBI)
        input_file (str): Downloaded nucleotide content (NCBI)

    Returns:
        dict: Compatible PAX-DB dataset
    """
    # Make a copy of original paxDB-information
    paxDBdataset = deepcopy(paxDbArr)

    # Read downloaded content in a batch
    for seq_record in SeqIO.parse(input_file, "fasta"):
        # Parse FASTA header line
        info = seq_record.id.split('|')
        paxDbId = info[2]

        # Set up and prepare data
        dummy = {'paxdbDatasetId' : 0, 'paxdbId' : 0, 'uniprotId' : '', 'dnaStart' : 0, 'dnaEnd' : 0, 'dnaStrand' : '', 'dnaSequence' : '', 'abundance' : 0.0}

        dummy['paxdbId'] = paxDbId
        dummy['uniprotId'] = info[1]
        dummy['dnaSequence'] = str(seq_record.seq)
        try:
            dummy['paxdbDatasetId'] = int(paxDbArr[paxDbId][1])
            dummy['abundance'] = float(paxDbArr[paxDbId][0])
        except KeyError:
            dummy['paxdbDatasetId'] = 0
            dummy['abundance'] = 0.0
            print(paxDbId)

        # Store information
        paxDBdataset[paxDbId] = dummy

    print(paxDBdataset)
    return paxDBdataset


def store_paxdb_information(paxDBdataset, proteins):
    """Store downloaded PAX-DB (NCBI) into Database (Postgres)

    Args:
        paxDbDataset (dict): Parsed Nucleotide information (PAX-DB/NCBI)
        proteins (dict): Parsed Protein information (PAX-DB/NCBI)

    Returns:
        bool: Indicator whether storing succeeded or nor
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()
    print("Connected!\n")

    # tell postgres to use more work memory
    work_mem = 2048

    # by passing a tuple as the 2nd argument to the execution function our
    # %s string variable will get replaced with the order of variables in
    # the list. In this case there is only 1 variable.
    # Note that in python you specify a tuple with one item in it by placing
    # a comma after the first variable and surrounding it in parentheses.
    #
    # Iterate over all prepared DNA and protein datasets
    cnt = 0
    for id in paxDBdataset.keys():
        if isinstance(paxDBdataset[id], list):
            continue
        cnt += 1
        try:
            locus_tag = id
            d = paxDBdataset[locus_tag]
            p = proteins[locus_tag]
            cursor.execute('INSERT INTO paxdb (paxdb_dataset_id,locus_tag,uniprot_id,protein_sequence,abundance,gb_keywords,dna_sequence,geninfo_id,dna_start,dna_end,dna_strand) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)',
                           (d['paxdbDatasetId'],locus_tag,p['primary-accession'],p['sequence'],d['abundance'],p['keywords'],d['dnaSequence'],d['uniprotId'],d['dnaStart'],d['dnaEnd'],d['dnaStrand']))
        except KeyError:
            locus_tag = id
            d = paxDBdataset[locus_tag]
            bioObj = load_fasta_data(unicode('>\n'+d['dnaSequence']))[0]
            translation = bioObj.seq.translate()
            cursor.execute('INSERT INTO paxdb (paxdb_dataset_id,locus_tag,uniprot_id,protein_sequence,abundance,gb_keywords,dna_sequence,geninfo_id,dna_start,dna_end,dna_strand) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)',
                           (d['paxdbDatasetId'],locus_tag,'',str(translation[0:-1]),d['abundance'],'',d['dnaSequence'],d['uniprotId'],d['dnaStart'],d['dnaEnd'],d['dnaStrand']))

        # Save data in chunks to DB
        if cnt%500 == 0:
            # Make the changes to the database persistent
            conn.commit()

    # Make the changes to the database persistent
    conn.commit()
    print("Make the changes to the database persistent")

    # Close communication with the database
    cursor.close()
    conn.close()

    return True


def select_paxdb_information(organism):
    """Request downloaded PAX-DB (NCBI) from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB

    Returns:
        dict, dict: Entries and keywords
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    keywords = {}; paxdb_id = host_organisms[organism]['paxdb_id']
    cursor.execute('SELECT paxdb_dataset_id, locus_tag, abundance, gb_keywords FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ' ORDER BY locus_tag ASC')
    entries = [dict(paxdb_id=row[0], locus_tag=row[1], abundance=row[2], keywords=row[3]) for row in cursor.fetchall()]
    for entry in entries:
        for keyword in re.sub('[{}"]', '', entry['keywords']).split(','):
            try:
                keywords[keyword] += 1;
            except KeyError:
                keywords[keyword] = 1;
    # Remove empty gb_key
    if len(keywords) > 0:
        keywords.pop('')
    # Close communication with the database
    cursor.close()
    conn.close()

    return entries, keywords


def get_paxdb_portion(organism, min, max):
    """Request PAX-DB portion for range from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB
        min (int): Range lower bound
        max (int): Range upper bound

    Returns:
        float, float, float: portion, log_min, log_max
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    paxdb_id = host_organisms[organism]['paxdb_id']
    sql = 'SELECT MAX(abundance) FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ';'
    cursor.execute(sql)
    number = cursor.fetchone()[0]

    # Logarithmic scaling
    # Position will be between 0 and 100
    minp = 0; maxp = 100;
    # The result should be between 0 an SQL:COUNT()
    minv = 0 # np.log(1)
    maxv = np.log(number)
    # Calculate adjustment factor
    scale = (maxv-minv) / (maxp-minp)
    log_min = np.exp(minv + scale*(min-minp)) - 1.0  # Need to get down to lowest abundance value (0.0)
    log_max = np.exp(minv + scale*(max-minp))

    sql = 'SELECT COUNT(abundance) FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ' AND abundance >= ' + str(log_min) + ' AND abundance <= ' + str(log_max) + ';'
    cursor.execute(sql)
    portion = cursor.fetchone()[0]
    sql = 'SELECT abundance FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ' AND abundance >= ' + str(log_min) + ' AND abundance <= ' + str(log_max) + ' ORDER BY abundance ASC;'
    cursor.execute(sql)
    entries = [row[0] for row in cursor.fetchall()]
    log_min, log_max = entries[0], entries[-1]

    # Close communication with the database
    cursor.close()
    conn.close()

    return {'portion': portion, 'log_min': log_min, 'log_max': log_max}


def get_paxdb_data(organism):
    """Request downloaded PAX-DB (NCBI) from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB

    Returns:
        dict: Entries
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    paxdb_id = host_organisms[organism]['paxdb_id']
    sql = 'SELECT paxdb_id, abundance FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ' ORDER BY abundance ASC;'
    cursor.execute(sql)
    entries = [row[1] for row in cursor.fetchall()]

    # Close communication with the database
    cursor.close()
    conn.close()

    return entries


def paxdb_select_organism_seqs_by_number(organism, most_abundant):
    """Request downloaded PAX-DB (NCBI) from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB
        most_abundant (int): Number of most abundant proteins

    Returns:
        list: Resulting CDS sequences
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    seqs = []; paxdb_id = host_organisms[organism]['paxdb_id']
    sql = "SELECT paxdb_dataset_id, locus_tag, abundance, dna_sequence FROM paxdb WHERE paxdb_dataset_id = " + str(paxdb_id) + " AND dna_sequence != '' ORDER BY abundance DESC LIMIT " + str(most_abundant) + ";"
    cursor.execute(sql)
    entries = [dict(paxdb_id=row[0], locus_tag=row[1], abundance=row[2], dna_sequence=row[3]) for row in cursor.fetchall()]

    # As single sequences
    for entry in entries:
        seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+entry['dna_sequence'])), "fasta"))
    # # Boost: As one large sequence
    # seq = ''
    # for entry in entries:
    #     seq += entry['dna_sequence']
    # seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+seq)), "fasta"))

    print(len(seqs))
    # Close communication with the database
    cursor.close()
    conn.close()

    return seqs


def paxdb_select_organism_seqs_by_range(organism, min, max):
    """Request PAX-DB organism sequence information for range from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB
        min (int): Range lower bound
        max (int): Range upper bound

    Returns:
        list: Resulting CDS sequences
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    paxdb_id = host_organisms[organism]['paxdb_id']
    sql = 'SELECT MAX(abundance) FROM paxdb WHERE paxdb_dataset_id = ' + str(paxdb_id) + ';'
    cursor.execute(sql)
    number = cursor.fetchone()[0]

    # Logarithmic scaling
    # Position will be between 0 and 100
    minp = 0; maxp = 100;
    # The result should be between 0 an SQL:COUNT()
    minv = 0 # np.log(1)
    maxv = np.log(number)
    # Calculate adjustment factor
    scale = (maxv-minv) / (maxp-minp)
    log_min = np.exp(minv + scale*(min-minp)) - 1.0  # Need to get down to lowest abundance value (0.0)
    log_max = np.exp(minv + scale*(max-minp))

    seqs = []
    sql = "SELECT paxdb_dataset_id, locus_tag, abundance, dna_sequence FROM paxdb WHERE paxdb_dataset_id = " + str(paxdb_id) + " AND dna_sequence != '' AND abundance >= " + str(log_min) + " AND abundance <= " + str(log_max) + " ORDER BY abundance DESC;"
    cursor.execute(sql)
    entries = [dict(paxdb_id=row[0], locus_tag=row[1], abundance=row[2], dna_sequence=row[3]) for row in cursor.fetchall()]

    # As single sequences
    for entry in entries:
        seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+entry['dna_sequence'])), "fasta"))
    # # Boost: As one large sequence
    # seq = ''
    # for entry in entries:
    #     seq += entry['dna_sequence']
    # seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+seq)), "fasta"))

    print(len(seqs))

    # Close communication with the database
    cursor.close()
    conn.close()

    return seqs


def paxdb_select_organism_seqs_by_keywords(organism, gb_keywords):
    """Request PAX-DB organism sequence information by keywords from Database (Postgres)

    Args:
        organism (str): Organism to search information in DB
        gb_keywords (list): List of GB keywords

    Returns:
        list: Resulting CDS sequences
    """
    # Define our connection string
    conn_string = "host='%s' dbname='HeterologousProteinExpression' user='%s' password='%s' port='%s'" % (DB_HOST, DB_USER, DB_PASS, DB_PORT)

    # get a connection, if a connect cannot be made an exception will be raised here
    conn = psycopg2.connect(conn_string)

    # conn.cursor will return a cursor object, you can use this cursor to perform queries
    cursor = conn.cursor()

    seqs = []; paxdb_id = host_organisms[organism]['paxdb_id']
    sql = "SELECT paxdb_dataset_id, locus_tag, abundance, dna_sequence FROM paxdb WHERE paxdb_dataset_id = " + str(paxdb_id) + " AND dna_sequence != '' AND gb_keywords LIKE '%" + str(gb_keywords[0]) + "%' ORDER BY abundance;"
    cursor.execute(sql)
    entries = [dict(paxdb_id=row[0], locus_tag=row[1], abundance=row[2], dna_sequence=row[3]) for row in cursor.fetchall()]

    # As single sequences
    for entry in entries:
        seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+entry['dna_sequence'])), "fasta"))
    # # Boost: As one large sequence
    # seq = ''
    # for entry in entries:
    #     seq += entry['dna_sequence']
    # seqs.append(SeqIO.read(StringIO(unicode(">"+entry['locus_tag']+"\n"+seq)), "fasta"))

    # Close communication with the database
    cursor.close()
    conn.close()

    return seqs


def load_fasta_file(input_file):
    """Load and parse FASTA-file with BioPython methods

    Args:
        input_file (str): Filename of the FASTA-file

    Returns:
        list: List of parsed FASTA-BioSeq-objects
    """
    nucleotides = []
    if os.path.exists(input_file):
        for record in SeqIO.parse(input_file, "fasta"):
            nucleotides.append(record)
        # Print number of parsed/found records
        print(len(nucleotides))
    return nucleotides


def load_fasta_data(input_string):
    """Load and parse FASTA-formatted string with BioPython methods

    Args:
        input_string (str): FASTA-formatted string with sequences

    Returns:
        list: List of parsed FASTA-BioSeq-objects
    """
    nucleotides = []
    if input_string:
        for record in SeqIO.parse(StringIO(input_string), "fasta"):
            nucleotides.append(record)
    return nucleotides


def store_fasta_data(input_obj, path):
    """Store FASTA-information to file

    Args:
        input_obj (object): FASTA-BioSeq objects
        path (str): Path to store FASTA

    Returns:
        bool
    """
    SeqIO.write(input_obj, path, "fasta")
    return True


def execute_rnafold(seq):
    """Execute RNAfold for passed DNA sequence

    Args:
        seq (str): FASTA-BioSeq object

    Returns:
        dict: RNAfold information
    """
    # Create working directory
    uniq = uniq_id(); workdir = BASE_APP_PATH + TMP_IMAGE_PATH + uniq
    os.chdir(TMP_IMAGE_PATH)
    os.makedirs(uniq)
    os.chdir(uniq)

    # Create temp sequence-file
    in_w = open('temp_rnafold_sequence.fasta', 'w')
    in_w.write(seq)
    in_w.close()

    # Call RNAfold
    # RNAfold -p < temp_rnafold_sequence.fasta > temp_rnafold_ouputs.txt
    process = subprocess.Popen(["export PATH='/usr/local/bin:/usr/bin:$PATH' && RNAfold -p < temp_rnafold_sequence.fasta > temp_rnafold_ouputs.txt"], shell=True)
    process.wait()

    # Read and parse result result for MFE-value
    f = open('temp_rnafold_ouputs.txt', 'r'); mfe_val = 0
    fileContent = f.read()            # Read in one String
    for m in re.finditer('\((-[0-9.]+)\)', fileContent):
        mfe_val = m.group(1)

    # Convert PS to PNG [ImageMagick]
    process = subprocess.Popen(["export PATH='/usr/local/bin:/usr/bin:$PATH' && convert -flatten -density 320 -colorspace RGB rna.ps rna.png"], shell=True)
    process.wait()
    process = subprocess.Popen(["export PATH='/usr/local/bin:/usr/bin:$PATH' && convert -flatten -density 400 -colorspace RGB dot.ps dot.png"], shell=True)
    process.wait()
    # Change back to base app path
    os.chdir(BASE_APP_PATH)

    return {'MFE': mfe_val, 'img': '/'.join([WEB_TMP_IMAGE_PATH + uniq, 'rna.png']), 'dot_img': '/'.join([WEB_TMP_IMAGE_PATH + uniq, 'dot.png']), 'dot_mat': '/'.join([workdir, 'dot.ps'])}


def parse_rnafold_base_pairing(path):
    """Parse RNAfold result for base pairing search

    Args:
        path (str): RNAFold result file

    Returns:
        list, dict: selection and matches (base pairing regions)
    """
    # This file contains the square roots of the base pair probabilities in the form
    # i  j  sqrt(p(i,j)) ubox
    # path = BASE_PATH_BIO + 'local_variables/RNAfold/gi|330443482:125124-126118|YBL050W_dp.ps'
    filepath = os.path.split(path); input_file = '/'.join(filepath)
    rnaFoldArr = []; maxVec = []; limitArr = []; sumVec = []

    # If file exists
    if os.path.isfile(input_file):
        # Reading content from specified file
        f = open(input_file, 'r')
        fileContent = f.read()            # Read in one String
        fileArr = fileContent.split("\n") # Split in lines

        # Define Simple grammar to match 'RNAfold base pairing probability data'
        intNum = Word(nums); floatNum = Word(nums + ".")
        macroDef = intNum.setResultsName("x") + pp.empty + intNum.setResultsName("y") + pp.empty + floatNum.setResultsName("bpp") + pp.empty + restOfLine.setResultsName("rest")
        # Parse RNAfold base pairing probability data
        for line in fileArr:
            tmp = {}
            for t, s, e in macroDef.scanString( line ):
                if t.rest == 'ubox':
                    tmp['x'] = int(t.x); tmp['y'] = int(t.y); tmp['bpp'] = float(t.bpp)
            if bool(tmp):
                rnaFoldArr.append(tmp)

        # Organize parsed data
        tmpArr = {}
        for dataset in rnaFoldArr:
            try:
                tmpArr[dataset['x']].append(dataset)
            except:
                tmpArr[dataset['x']] = []
                tmpArr[dataset['x']].append(dataset)

        # Get position probabilities by sum
        for key, subset in tmpArr.items():
            if subset:
                maxNum = np.max([dataset['bpp'] for dataset in subset])
                maxVec.append(maxNum)
        # print(maxVec)
        # Get position probabilities by limit / threshold
        for key, subset in tmpArr.items():
            if subset:
                limitVec = [dataset for dataset in subset if dataset['bpp'] > 0.9]
                limitArr.append(limitVec[0]) if limitVec != [] else None
        # print(limitArr)
        # Get position probabilities by max
        for key, subset in tmpArr.items():
            if subset:
                sumNum = np.sum([dataset['bpp'] for dataset in subset])
                sumVec.append(sumNum)
        # print(sumVec)

        # Identify pairing regions longer >= 3 nucleotides
        selection = []; matches = {'matches': {'rnafold': []}}
        predecessor = -1; row = []
        for idx, nuc in enumerate(limitArr):
            if predecessor == -1:
                predecessor = nuc
            elif nuc['y']+1 != predecessor['y'] or idx == len(limitArr)-1:
                if len(row) >= 3:
                    yVec = [dataset['y'] for dataset in row]
                    xVec = [dataset['x'] for dataset in row]
                    color = random_color()
                    # Prepare output
                    matches['matches']['rnafold'].append([0, min(yVec), max(yVec), color])
                    matches['matches']['rnafold'].append([0, min(xVec), max(xVec), color])
                    selection.append(row)
                row = []
                predecessor = -1
            elif nuc['y']+1 == predecessor['y']:
                row.append(nuc)
                predecessor = nuc

    return selection, matches


def seq_statistics(seqObj, codon_usage, rnafold_switch):
    """Generate sequence statistics to artificial sequence

    Args:
        seqObj (object): Sequence to generate statistics for
        codon_usage (dict): Mono Codon usage (ACU or RCU)
        rnafold_switch (bool): Turn RNAFold execution on/off

    Returns:
        dict: Statistics (CAI, GC, etc.)
    """
    # Cast SeqObj into normal String
    seq = str(seqObj.seq)
    seq_len = float(len(seq))

    # Calculate nucleotide frequencies
    nuc_freq = ['{:s}: {:03.3f}'.format(i, seq.count(i)/seq_len) for i in list("ACGT")]
    nuc_freq = ', '.join(nuc_freq)

    # CAI calculation
    rscu, relAdaptiveness = generate_rscu(codon_usage)
    cai_custom = calculate_cai(remove_stop_codons(seq), relAdaptiveness)

    cai = CodonUsage.CodonAdaptationIndex()
    # cai.set_cai_index(relAdaptiveness)
    cai_num = cai.cai_for_gene(remove_stop_codons(seq))
    # print(CodonUsage.SharpEcoliIndex)
    # print(relAdaptiveness)
    # print('{:03.3f} {:03.3f}'.format(cai_custom, cai_num))

    # Begin RNAFold analysis
    seq_highlight_matches = ''; rnafold_data = {'MFE': '', 'img': '', 'dot_img': '', 'dot_mat': ''}
    if rnafold_switch:
        start_time = datetime.datetime.now()
        rnafold_data = execute_rnafold(seq)
        regions, matches = parse_rnafold_base_pairing(rnafold_data['dot_mat'])
        if matches:
            seq_highlight_matches = make_alignment_matches([seqObj], matches, False)
        print('> Calc. time RNAfold: ' + str(millis(start_time)/1000) + ' sec.')

    return {'CAI': '{:03.3f}'.format(cai_custom), 'GC': '{:03.3f}'.format(GC(seq)), 'nuc_freq': nuc_freq, 'MFE': rnafold_data['MFE'], 'rnafold_img': rnafold_data['img'], 'rnafold_dot_img': rnafold_data['dot_img'], 'bpp_sites': seq_highlight_matches}


def make_nuc_record(protein_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = protein_record.seq.translate(cds=True),
                     id = "trans_" + protein_record.id,
                     description = "translation of CDS, using default table")


def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_record.seq.translate(cds=True),
                     id = "trans_" + nuc_record.id,
                     description = "translation of CDS, using default table")


def make_cu_comparison_csv(setup_one, setup_two, filename):
    """Write CSV into TMP_PATH_DATA with relative-codon values (CUs)"""

    relCodonsOne = get_rel_codon_usage(setup_one['G'], setup_one['nl2codon'])
    relCodonsTwo = get_rel_codon_usage(setup_two['G'], setup_two['nl2codon'])
    synRelCodonsOne = calculate_rel_syn_codon_usage(relCodonsOne)
    synRelCodonsTwo = calculate_rel_syn_codon_usage(relCodonsTwo)

    csv_data = []
    for k, v in synRelCodonsOne.iteritems():
        if math.isnan(v): v = float(0)
        if math.isnan(synRelCodonsTwo[k]): synRelCodonsTwo[k] = float(0)
        csv_data.append([k, v, synRelCodonsTwo[k], codonMapping['forward'][k]])
    sorted_csv = sorted(csv_data, key=lambda el: el[0])
    with open('%s%s%s.csv' % (APP_TMP_PATH, 'img/', filename), 'wb') as f:
        writer = csv.writer(f, delimiter=';')  # dialect='excel')
        writer.writerows(sorted_csv)

    return None


def make_cu_comparison_gnu(csv_file, x_title, y_title):
    """Generate GNUplots into TMP_PATH_DATA - CU comparison"""
    gnu_template = ''; image_path = "%s%s" % (APP_TMP_PATH, 'img/')
    infilepath = "%s%s" % (image_path, 'codon_usage_versus.gnu_tmpl')
    if os.path.isfile(infilepath):
        f = open(infilepath, 'r')
        gnu_template = f.read()  # Read in one String
        f.close()

    gnu_template = re.sub(r"###OUTPUT###", "%s%s" % (image_path, csv_file), gnu_template)
    gnu_template = re.sub(r"###FILE###",   "%s%s" % (image_path, csv_file), gnu_template)
    gnu_template = re.sub(r"###XTITLE###", x_title, gnu_template)
    gnu_template = re.sub(r"###YTITLE###", y_title, gnu_template)

    outfilepath = "%s%s" % (image_path, 'codon_usage_versus.gnu')
    in_w = open(outfilepath, 'w')
    in_w.write(gnu_template)
    in_w.close()

    process = subprocess.Popen(["export PATH='/usr/local/bin:/usr/bin:$PATH' && %scodon_usage_versus.gnu " % image_path], shell=True)
    process.wait()
    # process = subprocess.Popen(["export PATH='/usr/local/bin:/usr/bin:$PATH' && convert -density 1200 -resize 200x200 %s%s.svg %s%s.png" % (image_path,csv_file,image_path,csv_file)], shell=True)
    # process.wait()

    return None
