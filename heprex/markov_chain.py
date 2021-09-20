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
# Created by dsimm1 on 04/08/15.
#

# Create Markov-Chain in form of a fully-interconnected graph
# Support for: Nodes, Edges with stationary and transitions probabilities

# Inclusions
from bio_utilities import *
from networkx.drawing.nx_pydot import write_dot
import networkx as nx
import random
from bisect import bisect
from scipy.spatial import distance


def create_markov_chain(all_codons, codonsRel, dicodonsRel):
    """Create the Markov Chain model on basis of the codon_set, mono-codon-usage (RCU) and di-codon-usage (RCU)
    - State (node) probabilities: Mono-Codon-Usage (RCU)
    - Transition probabilities: Node A -> Node B
                    transition: Di-Codon-Usage (RCU)
             normed_transition:

    Args:
        all_codons (list): List of all valid codons plus stop codons
        codonsRel (dict): Mono-Codon-Usage (RCU)
        dicodonsRel (dict): Di-Codon-Usage (RCU)

    Returns:
        DiGraph, dict, dict: Graph-object G, translation-dicts: nl2codon & nl2id
    """
    G = nx.DiGraph()
    node_num = 0

    # Create nodes corresp. to codons [61]
    for codon in all_codons:
        node_num += 1
        G.add_node(node_num, codon=codon, prob=codonsRel[codon])
    # print(G.nodes(data=True))

    # Create mapping [codon: id] (= node_list)
    node_list = nx.get_node_attributes(G, 'codon')
    node_list = dict((v, k) for k, v in node_list.items())

    # Set edges with transition probabilities (full net)
    for n in G:
        for codon, freq in dicodonsRel[G.node[n]['codon']].items():
            G.add_edge(n, node_list[codon], transition=freq, normed_transition=0.0)
    # print(G.edges(data=True))

    # Add normalized transitions for edges
    normed_transitions = {}
    for src_codon, src_aa in codonMapping['forward'].items():
        for end_aa, end_codons in codonMapping['backward_extended'].items():
            # Get numerical ids for alternative codons
            freqs = dict(((node_list[src_codon], node_list[end_codon]), dicodonsRel[src_codon][end_codon]) for end_codon in end_codons)
            # General normalization over Mono-Codon-Usage
            sums = sum(list(G.get_edge_data(node_list[src_codon], num)['transition'] for num in range(1, 65)))

            # Normalize the actual possible codon-frequences (Conditional probability)
            freqs_norm = freqs.copy()
            for el in freqs:
                freqs_norm[el] /= float(sums if sums > 0.0 else 1.0)
            normed_transitions.update(freqs_norm)
    # print normed_transitions
    nx.set_edge_attributes(G, 'normed_transition', normed_transitions)

    # Create mapping [codon: node_id] (= node_list)
    node_list = nx.get_node_attributes(G, 'codon')
    nl2codon = node_list
    node_list = dict((v, k) for k, v in node_list.items())
    nl2id = node_list

    # Feature: Write Markov Chain to .dot-file
    # nx.draw_graphviz(G)
    # nx.write_dot(G, BASE_APP_PATH + BASE_PATH_BIO + 'mc.dot')

    # Return values
    return G, nl2codon, nl2id


def get_rel_codon_usage(G, nl2codon):
    """Get Mono-RCU from passed MC-model G

    Args:
        G (DiGraph): DiGraph object of current Markov chain
        nl2codon (dict): Mapping of node-ids to codon-names

    Returns:
        dict: Relative mono codon usage (stationary probs. of nodes)
    """
    relCodons = {}
    # Start with first codon of sequence; Codon like 'ATG', but random:
    node_probs = nx.get_node_attributes(G, 'prob')
    for i, prob in node_probs.iteritems():
        relCodons.__setitem__(nl2codon[i], prob)
    return relCodons


def get_rel_dicodon_usage(G, nl2codon):
    """Get Di-RCU from passed MC-model G

    Args:
        G (DiGraph): DiGraph object of current Markov chain
        nl2codon (dict): Mapping of node-ids to codon-names

    Returns:
        dict: Relative di-codon usage (transition probs. between nodes/egdes)
    """
    relCodons, relDiCodons = initialize_codons()
    edge_probs = nx.get_edge_attributes(G, 'transition')
    for pair, prob in edge_probs.iteritems():
        relDiCodons[nl2codon[pair[0]]][nl2codon[pair[1]]] = prob
        # relDiCodons[nl2codon[pair[0]]].__setitem__([nl2codon[pair[1]]], prob)
    return relDiCodons


def weighted_choice(choices):
    """Select one out of a list of weighted elements

    Args:
        choices (dict): Codon usages / Probabilities

    Returns:
        str: Key of choices dict (codon)
    """
    total = sum(w for c, w in choices.items())
    if total > 0:
        r = random.uniform(0, total)
        upto = 0
        # print(choices)
        for c, w in choices.items():
            if upto + w > r:
                return c
            upto += w
    # Can only occur, when all entries are Zero
    else:
        return random.choice(choices.keys())

    assert False, "Shouldn't get here"


########################################################################################################################
# Information theory approach: Entropy-based
########################################################################################################################

def calculate_model_entropy(G, nl2id):
    """Calculate entropy of given MC model (parameters)
    - Relative codon-usage
    - Relative dicodon-usage

    Args:
        G (DiGraph): DiGraph object of current Markov chain
        nl2id (dict): Mapping of codon-names to node-ids

    Returns:
        float: Model entropy
    """
    sum_list = []
    for src_codon, src_aa in codonMapping['forward'].items():
        # sum(list(G.get_edge_data(nl2id[src_codon], num)['transition'] * log2(1 / float(G.get_edge_data(nl2id[src_codon], num)['normed_transition'])) for num in range(1, 64)))
        for num in range(1, 64):
            edge = G.get_edge_data(nl2id[src_codon], num)
            try:
                sum_list.append( edge['transition'] * log2(1/float(edge['normed_transition'])) )
            except ZeroDivisionError:
                sum_list.append(float(0))
    # Return summed entropy
    return sum(sum_list)


def calculate_entropy_distance(sequence, entropy, G, nl2id):
    """Calculate distance measure from model entropy

    Args:
        sequence (dict): DNA sequence as SeqObj with DiCodon-probability-information
        entropy (float): MC model entropy
        G (DiGraph): DiGraph object of current Markov chain
        nl2id (dict): Mapping of codon-names to node-ids

    Returns:
        float: Distance (value over 0.0)
    """
    start_codon = sequence[0][0]
    seq_prob = log2(1/G.node[nl2id[start_codon]]['prob']) + np.sum(list(log2(1/v[1]) if v[1] > 0 else 1.0 for v in sequence.values()))
    avg_seq_prob = seq_prob/len(sequence)
    distance = entropy - avg_seq_prob

    return abs(distance)


def calculate_sample_entropy(sequence, G, nl2id):
    """Calculate sample entropy

    Args:
        sequence (dict): DNA sequence as SeqObj with DiCodon-probability-information
        G (DiGraph): DiGraph object of current Markov chain
        nl2id (dict): Mapping of codon-names to node-ids

    Returns:
        float: Sample entropy (value over 0.0)
    """
    start_codon = sequence[0][0]
    seq_prob = log2(1/G.node[nl2id[start_codon]]['prob']) + np.sum(list(log2(1/v[1]) if v[1] > 0 else 1.0 for v in sequence.values()))
    avg_seq_prob = seq_prob/len(sequence)

    return avg_seq_prob


########################################################################################################################
# Mathematical approach: Matrix-based comparison
########################################################################################################################

def calculate_matrix_distance(matrixA, matrixB):
    """Calculate distance between to matrices A and B

    Args:
        matrixA (dict):
        matrixB (dict):

    Returns:
        float: Sample entropy (value over 0.0)
    """
    # Y = spatial.distance.cdist(dict2matrix(matrixA), dict2matrix(matrixB), lambda u, v: np.sqrt(((u - v) ** 2).sum().sum()))
    Y = distance.cdist(dict2matrix(matrixA), dict2matrix(matrixB), 'euclidean')

    return np.trace(Y)


def dict2matrix(matrix):
    """Convert a dict into a matrix

    Args:
        matrix (dict):

    Returns:
        list: List of lists (vectors)
    """
    sort_matrix = []; sort_list = sorted(matrix)
    for key in sort_list:
        sort_vector = []
        for sec_key in sort_list:
            sort_vector.append(float(matrix[key][sec_key]))
        sort_matrix.append(sort_vector)
    return sort_matrix


########################################################################################################################
# Sequence-Generators with Probability-information: Alternative DNA-sequence, fully-random DNA-Sequence
########################################################################################################################

def generate_alternative_dna_sequence(seq_record, G, nl2id, nl2codon, extended = True):
    """Generate on basis of a predefined Markov-chain model synonymous DNA sequences for heterologous protein expression

    Args:
        seq_record (SeqRecord): Aminoacid (AA) sequence; FASTA parsed Bio.Seq object
        G (DiGraph): Markov-chain in form of full graph
        nl2id (dict): Codon to Graph-ID-Mapping
        nl2codon (dict): Graph-ID to Codon-Mapping

    Returns:
        dict: Sequence (alternative) with probability information
    """
    alternative_seq = {}
    # Start with first codon of sequence and its frequency
    node_probs = nx.get_node_attributes(G, 'prob')
    states = codonMapping['backward_extended'][seq_record.seq[0]]
    codon_ids = list(nl2id[codon] for codon in states)
    codon_freqs = {id: node_probs[id] for id in codon_ids}
    codon_id = weighted_choice(codon_freqs)
    codon = nl2codon[codon_id]; active_state = codon
    alternative_seq[0] = (active_state, node_probs[codon_id])
    # alternative_seq[0] = (active_state, 1.0)
    for n in range(1, len(seq_record), 1):
        # Determine current amino acid
        aa = seq_record.seq[n]
        # Determine encoded codons
        codons = codonMapping['backward'][aa]
        # Get numerical ids for alternative codons
        codon_ids = list(nl2id[codon] for codon in codons)
        freqs_norm = dict((id, G.get_edge_data(nl2id[active_state], id)['normed_transition']) for id in codon_ids)
        codon_id = weighted_choice(freqs_norm)
        # print codon_id
        codon = nl2codon[codon_id]
        # print codon
        alternative_seq[n] = (codon, freqs_norm[codon_id])
        active_state = codon

    return alternative_seq


def generate_random_sequence(length, G, nl2id, nl2codon, extended=True):
    """Generate on basis of a predefined Markov-chain model completely random protein sequences for heterologous protein
    expression

    Args:
        length (int): Length of output sequence
        G (DiGraph): Markov-chain in form of full graph
        nl2id (dict): Codon to Graph-ID-Mapping
        nl2codon (dict): Graph-ID to Codon-Mapping

    Returns:
        dict: Sequence (alternative) with probability information
    """
    random_seq = {}

    # Get numerical ids for codons
    codon_ids = list(nl2id[codon] for codon in all_codons)

    # Start with first codon of sequence; Codon like 'ATG', but random:
    node_probs = nx.get_node_attributes(G, 'prob')
    codon_id = weighted_choice(node_probs)
    codon = nl2codon[codon_id]; active_state = codon
    random_seq[0] = (active_state, node_probs[codon_id])

    # Generate sequence
    for n in range(1, length, 1):
        # Gather transition frequences of MC model for active state [amino acid]
        freqs_norm = dict((id, G.get_edge_data(nl2id[active_state], id)['normed_transition']) for id in codon_ids)
        # Choose follow-up codon
        codon_id = weighted_choice(freqs_norm)
        codon = nl2codon[codon_id]
        random_seq[n] = (codon, freqs_norm[codon_id])
        active_state = codon

    return random_seq


def generate_probability_enriched_sequence(seq_record, G, nl2id, cu='di'):
    """Add to passed DNA-Sequence (SeqRecord) probability information from a predefined Markov-chain model
    => Mono- | Di-Codon-Usage

    Args:
        seq_record (SeqRecord): DNA sequence; FASTA parsed Bio.Seq object
        G (DiGraph): Markov-chain in form of full graph
        nl2id (dict): Codon to Graph-ID-Mapping
        cu (str): Codon-usage selector ['mono'|'di']

    Returns:
        SeqRecord: The passed SeqRecord object enriched with probability information
    """

    # Add additional field to BioSeq:SeqRecord
    seq_record.seq_info = {}
    # Start with first codon of sequence; Codon like 'ATG', but random:
    start_codon = seq_record.seq[0:3]
    node_probs = nx.get_node_attributes(G, 'prob')
    codon_id = nl2id[start_codon]
    seq_record.seq_info[0] = (start_codon, node_probs[codon_id])

    # Di-Codon-Usage
    if cu == 'di':
        # Generate sequence info object
        i = 1; step = 3;   # 1 (all 3 frames) | 3 (specified "coding" frame)
        for n in range(0, len(seq_record)-step, step):
            first_codon = seq_record.seq[n:n + 3]
            sec_codon = seq_record.seq[n + 3:n + 6]
            freq_norm = G.get_edge_data(nl2id[first_codon], nl2id[sec_codon])['normed_transition']
            seq_record.seq_info[i] = (sec_codon, freq_norm)
            i += 1
    # Mono-Codon-Usage
    else:
        # Generate sequence info object
        i = 1;
        step = 3;  # 1 (all 3 frames) | 3 (specified "coding" frame)
        for n in range(0, len(seq_record) - step, step):
            # first_codon = seq_record.seq[n:n + 3]
            sec_codon = seq_record.seq[n + 3:n + 6]
            freq_norm = node_probs[nl2id[sec_codon]]
            seq_record.seq_info[i] = (sec_codon, freq_norm)
            i += 1

    return seq_record

