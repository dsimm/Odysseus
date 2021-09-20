#!/usr/bin/env python
#
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
#

__author__ = 'dsimm'

# General Inclusions
from difflib import SequenceMatcher

# Local Inclusions
from .bio_utilities import *
from .help_utilities import *
from .bio_seq_encoder import *
from .markov_chain import *


def prepare_requirements():
    """
    Prepares MC parameter-sets (var. CUs) for all ref. CDS-data (FASTA) in hosts

    Returns:
        dict: Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G
    """

    hashes = {}
    # Host organisms
    profiles = OrderedDict([
        ('e.coli', host_organisms['e.coli']['genome']),
        ('s.cerevisiae', host_organisms['s.cerevisiae']['genome']),
        ('s.cerevisiae_high', host_organisms['s.cerevisiae_high']['genome']),
        ('s.cerevisiae_mid', host_organisms['s.cerevisiae_mid']['genome']),
        ('s.cerevisiae_low', host_organisms['s.cerevisiae_low']['genome']),
        ('s.cerevisiae_high_weighted', host_organisms['s.cerevisiae_high']['genome']),
        ('s.cerevisiae_mid_weighted', host_organisms['s.cerevisiae_mid']['genome']),
        ('s.cerevisiae_low_weighted', host_organisms['s.cerevisiae_low']['genome']),
        ('s.cerevisiae_high_weighted_invert', host_organisms['s.cerevisiae_high']['genome']),
        ('s.cerevisiae_mid_weighted_invert', host_organisms['s.cerevisiae_mid']['genome']),
        ('s.cerevisiae_low_weighted_invert', host_organisms['s.cerevisiae_low']['genome']),
        ('a.thaliana', host_organisms['a.thaliana']['genome']),
        # ('h.sapiens', host_organisms['h.sapiens']['genome']),
    ])
    for profile, file_path in profiles.items():
        # Initialize variables
        start_time = datetime.datetime.now()
        # Begin calculation
        print('Host: ' + profile + ' - ' + file_path)
        # Create Codon-Usages (Mono-, Di-Codon-Usage)
        print("* Create absolute Codon-Usages")
        if profile.find('weight') == -1:
            codons_tmp, dicodons_tmp = calculate_abs_codon_usages(BASE_PATH_BIO + file_path)
        else:
            codons_tmp, dicodons_tmp = calculate_abs_codon_usages_weighted(BASE_PATH_BIO + file_path)
        if profile == host_organisms[profile]['ref']:
            print("* Create relative Codon-Usages")
            codonsRel = calculate_rel_codon_usage(codons_tmp)
            dicodonsRel = calculate_rel_dicodon_usage(dicodons_tmp)
            # Create Markov Chain
            print("* Create Markov Chain")
            G, nl2codon, nl2id = create_markov_chain(all_codons, codonsRel, dicodonsRel)
            # Prepare hash-object (dict) of variable-copies (avoid value-overwriting)
            hashes[profile] = {'codons': codons_tmp.copy(), 'dicodons': deepcopy(dicodons_tmp), 'G': G.copy(),
                               'nl2codon': nl2codon.copy(), 'nl2id': nl2id.copy()}
            # Measure time
            print('> Calc. time: ' + str(millis(start_time) / 1000) + ' sec.')
        else:
            if profile.find('invert') == -1:
                hashes[profile] = recalculate_requirements(
                    codons_tmp, dicodons_tmp,
                    hashes[host_organisms[profile]['ref']]['codons'], hashes[host_organisms[profile]['ref']]['dicodons']
                )
            else:
                hashes[profile] = recalculate_requirements_invert(
                    codons_tmp, dicodons_tmp,
                    hashes[host_organisms[profile]['ref']]['codons'], hashes[host_organisms[profile]['ref']]['dicodons']
                )

    return hashes


def recalculate_requirements(codons_acu, dicodons_acu, ref_codons_acu, ref_dicodons_acu):
    """
    Prepares MC parameter-sets (various CUs) for passed absolute codon-usage (mono, di)
    Helper: Fills up empty CU-fields from reference CUs (intended from entire genome CUs)

    Args:
        codons_acu (dict): Absolute mono-codon-usage
        dicodons_acu (int): Absolute di-codon-usage
        ref_codons_acu (dict): Reference absolute mono-codon-usage
        ref_dicodons_acu (dict): Reference absolute di-codon-usage

    Returns:
        dict: Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G
    """
    # Initialize variables
    start_time = datetime.datetime.now()

    print("* Create relative Codon-Usages")
    codons_rcu = calculate_rel_codon_usage(codons_acu)
    dicodons_rcu = calculate_rel_dicodon_usage(dicodons_acu)

    # # Simulate Zero-cases
    # pre_codons_rcu['ATG'] = 0.0;
    # pre_codons_rcu['GCT'] = 0.0; pre_codons_rcu['GCC'] = 0.0;
    # pre_codons_rcu['GCA'] = 0.0; pre_codons_rcu['GCG'] = 0.0;
    # Check: Null frequencies and fill up with Genome-RCU
    if not check_codons(codons_rcu, dicodons_rcu):
        ref_codons_rcu = calculate_rel_codon_usage(ref_codons_acu)
        ref_dicodons_rcu = calculate_rel_dicodon_usage(ref_dicodons_acu)
        codons_rcu = calculate_rel_codon_usage(codons_rcu_fill_up(codons_rcu, ref_codons_rcu))
        dicodons_rcu = calculate_rel_dicodon_usage(dicodons_rcu_fill_up(dicodons_rcu, ref_dicodons_rcu))

    # Create Markov Chain
    print("* Create Markov Chain")
    G, nl2codon, nl2id = create_markov_chain(all_codons, codons_rcu, dicodons_rcu)
    # Prepare hash-object (dict) of variable-copies (avoid value-overwriting)
    setup_hash = {'codons': codons_acu.copy(), 'dicodons': deepcopy(dicodons_acu), 'G': G.copy(),
                  'nl2codon': nl2codon.copy(), 'nl2id': nl2id.copy()}
    # Measure time
    print('> Calc. time: ' + str(millis(start_time) / 1000) + ' sec.')

    return setup_hash


def recalculate_requirements_invert(codons_acu, dicodons_acu, ref_codons_acu, ref_dicodons_acu):
    """
    Prepares MC parameter-sets inverted (var. CUs) for passed abs. codon-usage (mono, di)
    Helper: Fills up empty CU-fields from reference CUs (intended from entire genome CUs)

    Args:
        codons_acu (dict): Absolute mono-codon-usage
        dicodons_acu (int): Absolute di-codon-usage
        ref_codons_acu (dict): Reference absolute mono-codon-usage
        ref_dicodons_acu (dict): Reference absolute di-codon-usage

    Returns: Setup-Dictionary
        dict: Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G
    """
    # Initialize variables
    start_time = datetime.datetime.now()

    print("* Create relative Codon-Usages")
    codons_rcu = calculate_rel_codon_usage(codons_acu)

    # Check: Null frequencies and fill up with Genome-RCU
    if not check_codons_rcu(codons_rcu):
        ref_codons_rcu = calculate_rel_codon_usage(ref_codons_acu)
        codons_rcu = calculate_rel_codon_usage(codons_rcu_fill_up(codons_rcu, ref_codons_rcu))
    # Invert
    codons_rcu_invert = calculate_rel_codon_usage_invert(codons_rcu, ref_codons_rcu)
    dicodons_rcu_invert = calculate_rel_dicodon_usage_invert(dicodons_acu, ref_dicodons_acu, codons_acu)

    # Create Markov Chain
    print("* Create Markov Chain")
    G, nl2codon, nl2id = create_markov_chain(all_codons, codons_rcu_invert, dicodons_rcu_invert)
    # Prepare hash-object (dict) of variable-copies (avoid value-overwriting)
    setup_hash = {'codons': codons_acu.copy(), 'dicodons': deepcopy(dicodons_acu), 'G': G.copy(),
                  'nl2codon': nl2codon.copy(), 'nl2id': nl2id.copy()}
    # Measure time
    print('> Calc. time: ' + str(millis(start_time) / 1000) + ' sec.')

    return setup_hash


def create_heterologous_sequences(filepath, nos, setup):
    """
    Wrapper-function:
        Generate synonymous gene sequences for the passed sequence (FASTA) with pre-compiled setup-dict.

    Args:
        filepath (str): filepath to source CDS.
        nos (int): number of sequences.
        setup (dict): Setup holding model configuration (Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G)

    Returns: List
        list: BioPython sequence objects (FASTA formatted)
    """
    # Initialize
    nos = 10 if (int(nos) < 0) else nos

    # Main part | Specification DNA sequence (protein) of interest [Insulin, Hexokinase, etc.]
    record = SeqIO.read(filepath, "fasta")

    # Begin iterative calculation of new heterologous sequences based on MC
    print("###############################")
    print("Alternative nucleotide sequence")
    sequences = []
    for n in range(0, nos):
        # Generate synonymous DNA sequence
        seq_mc = generate_alternative_dna_sequence(record, setup['G'], setup['nl2id'], setup['nl2codon'])

        # Calc entropy and distance measure
        # seq_test = generate_random_sequence(100, setup['G'], setup['nl2id'], setup['nl2codon'])
        # entropy = calculate_model_entropy(setup['G'], setup['nl2id'])
        # calculate_distance(seq_mc, entropy, setup['G'], setup['nl2id'])

        seq_mc_combined = ''.join(list(v[0] for v in seq_mc.values()))
        seq_fasta = unicode(">" + record.description + ("|alternative:%s\n" % n) + seq_mc_combined)
        # Create BioPython FASTA sequence
        bio_seq = SeqIO.read(StringIO(seq_fasta), "fasta")
        bio_seq.description = record.description + ("|alternative:%s\n" % n)
        sequences.append(bio_seq)

    return sequences


def create_heterologous_sequences_from_seq_record(record, nos, setup):
    """
    Wrapper-function:
        Generate synonymous gene sequences for the passed sequence object (BioPython) with pre-compiled setup-dict.

    Args:
        record (object): BioPython sequence object.
        nos (int): number of sequences.
        setup (dict): Setup holding model configuration (Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G)

    Returns:
        list: BioPython sequence objects (FASTA formatted)
    """
    # Initialize
    nos = 10 if (int(nos) < 0) else nos
    sequences = []
    for n in range(0, nos):
        # Generate synonymous DNA sequences
        seq_mc = generate_alternative_dna_sequence(record, setup['G'], setup['nl2id'], setup['nl2codon'])
        seq_mc_combined = ''.join(list(v[0] for v in seq_mc.values()))
        seq_fasta = unicode(">" + record.description + ("|alternative:%s\n" % n) + seq_mc_combined)
        # Create BioPython FASTA sequence
        bio_seq = SeqIO.read(StringIO(seq_fasta), "fasta")
        bio_seq.description = record.description + ("|alternative:%s\n" % n)
        sequences.append(bio_seq)

    return sequences


def search_restricted_areas(artificial_seqs):
    """
    Wrapper-function: Filter for sequences by commercially available enzyme sites

    Args:
        artificial_seqs (list): BioPython sequence objects

    Returns: List
        list: Filtered (commercially) list of regex-matches
    """
    bairoch_arr = parse_bairoch_file(BAIROCH_FILE_PATH)
    bairoch_avail_enzymes = bairoch_filter_commercial_enzymes(bairoch_arr)
    # Results from Restriction Site identification
    print("\n#############################################\nResults from Restriction Site identification")
    return identify_restriction_sites(artificial_seqs, bairoch_avail_enzymes)


def matrix_distance_measure(sequence, setup):
    """
    Calculate matrix distance measure from generated sequence object

    Args:
        sequence (object): Sequence
        setup (dict): Setup holding model configuration (Mono-ACU, Di-ACU, MC-graph G, codon-translation lists of G)
    """
    # Generate custom codon usage from sequenceObj
    print('Estimate codon usages:')
    seq_obj_arr = [sequence]
    codons_custom, dicodons_custom = calculate_abs_codon_usages_string(seq_obj_arr)

    # Create Markov-Chain from Sequence
    setup_custom = recalculate_requirements(codons_custom, dicodons_custom)

    # Calculate target value: dSeq
    dicodonsRel = calculate_rel_dicodon_usage(setup['dicodons'])
    dSeq = calculate_matrix_distance(dicodonsRel, calculate_rel_dicodon_usage(setup_custom['dicodons']))

    # Generate random sequences (strings)
    sample_size = 10000;
    sample_distances = []
    for i in range(1, sample_size):
        length = random.randint(100, 1000)
        seq_info = generate_random_sequence(length, setup_custom['G'], setup_custom['nl2id'], setup_custom['nl2codon'],
                                            extended=True)
        seq = ''.join(list(v[0] for v in seq_info.values()))
        seqfasta = unicode((">sample:%s\n" % i) + seq)
        seq_obj_arr = [SeqIO.read(StringIO(seqfasta), "fasta")]
        codons_sample, dicodons_sample = calculate_abs_codon_usages_string(seq_obj_arr)
        # Calc distances from Genome
        sample_distances.append(calculate_matrix_distance(dicodonsRel, calculate_rel_dicodon_usage(dicodons_sample)))

    # F-Test: More than 95% of Samples lie beyond target-value dSeq
    print("Target value: dSeq = %f" % dSeq)
    print(sorted(sample_distances, reverse=True))
    position = bisect(sorted(sample_distances), dSeq)
    threshold = position / float(sample_size)
    print("%d - %f" % (position, 1.0 - threshold))
    print("F-Test: Outcome: %2.2f perc. samples lie under reference sequence." % (1.0 - threshold))
    print("Passed" if 1.0 - threshold < 95.0 else "Failed")

    return -1
