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
# Created by dsimm on 28/01/16.
#

# Inclusions
# Load MatplotLib package
from matplotlib.pylab import *

# Load System packages
import random, re

# Local Inclusions
from .globals import *


def make_matrix(alignment):
    """Convert a list of SeqRecords into codon-list-matrix

    Args:
        alignment (list): List of sequences

    Returns:
        list, int: Matrix and length
    """

    matrix = []; sequence_length = 0; step = 3
    for idx, sequence in enumerate(alignment):
        for n in range(0, len(sequence), step):
            codon = str(sequence.seq[n:n+step])
            try:
                matrix[idx].append(codon)
            except IndexError:
                matrix.append([codon])
        sequence_length = len(matrix[idx])

    # Combine and return data
    return matrix, sequence_length


def add_aa_colors(sequence):
    """Highlight codons (by amino acids) of a SeqRecord (add surrounding div-tags to the letters)

    Args:
        sequence (str): Sequence

    Returns:
        str: Coloured sequence (HTML)
    """
    title = ''; coloured_sequence = ''; step = 3
    for n in range(0, len(sequence), step):
      codon = sequence.seq[n:n+step]
      aa = codonMapping['forward'][codon]
      bordered = 'unbordered' if re.match('[a-zA-Z]', aa) == None else 'bordered'
      coloured_sequence += "<div class='ui-sortable-element sequence nuc-" + aa.capitalize() + " " + bordered + "' id='' style='left: 0px; z-index: 100;'>" + str(codon).upper() + "</div>"

    # Combine and return data
    return title + coloured_sequence


def add_diff_colors(sequence, consensus_seq):
    """Highlight codons (by differences) of a SeqRecord (add surrounding div-tags to the letters)

    Args:
        sequence (object): Sequence (SeqRecord)
        consensus_seq (str): Consensus sequence of alignment

    Returns:
        str: Coloured sequence (HTML)
    """
    colored_sequence = ''; step = 3
    for n in range(0, len(sequence), step):
      codon = str(sequence.seq[n:n+step])
      cons_codon = consensus_seq[n:n+step]
      aa = 'Z' if codon == cons_codon else codon[2]
      bordered = 'unbordered' if re.match('[a-zA-Z]', aa) == None else 'bordered'
      colored_sequence += "<div class='ui-sortable-element sequence nuc-" + aa.lower() + " " + bordered + "' id='' style='left: 0px; z-index: 100;'>" + str(codon).upper() + "</div>"

    # Combine and return data
    return colored_sequence


def add_match_colors(sequence, matches):
    """Highlight codons (by matches) of a SeqRecord (add surrounding div-tags to the letters)

    Args:
        sequence (object): Sequence (SeqRecord)
        matches (list): List of match pairs e.g. base pairings

    Returns:
        str: Coloured sequence (HTML)
    """
    colored_sequence = ''; step = 1
    for n in range(0, len(sequence), step):
        aa = 'Z'; color = ''
        for match in matches:
            if match[1] <= n < match[2]:
                if len(match) > 3:
                    color = ' background-color: ' + str(match[3])
                else:
                    aa = 'Q'
        # Color atom (nucleotide, aminoacid)
        atom = str(sequence.seq[n:n+step])
        bordered = 'bordered' if aa == 'Z' else 'bordered-match'
        colored_sequence += "<div class='ui-sortable-element sequence nuc-" + aa.capitalize() + " " + bordered + "' id='' style='left: 0px; z-index: 100;" + color + "'>" + atom.upper() + "</div>"

    # Return data
    return colored_sequence


def make_alignment_aa(alignment):
    """Generate a MSA view with codons colored by aminoacids

    Args:
        alignment (list): List of sequences (SeqRecords)

    Returns:
        str: Coloured sequences (HTML)
    """
    lines = []
    for sequence in alignment:
        lines.append(add_aa_colors(sequence))

    # Join highlighted HTML into one string
    data = '</div><div class="oneliner">'.join(lines)
    data = '<div class="oneliner">' + data + '</div>'
    html_code = "<div style='padding: 4px;'>" + data + "</div>"
    return html_code


def make_alignment_diffs(alignment):
    """Generate a MSA view with highlighted codons colored by difference

    Args:
        alignment (list): List of sequences (SeqRecords)

    Returns:
        str: Coloured sequences (HTML)
    """
    lines = []; seq_num = len(alignment); consensus_seq = ''
    matrix, seq_len = make_matrix(alignment)
    for m in range(0, seq_len, 1):
        codons = {}
        for n in range(0, seq_num, 1):
            try:
                codons[matrix[n][m]] += 1
            except KeyError:
                codons[matrix[n][m]] = 1
        max = {'key': '', 'val': 0}
        for key, val in codons.items():
            if val > max['val']:
                max['key'] = key
                max['val'] = val
        consensus_seq += max['key']

    for sequence in alignment:
        lines.append(add_diff_colors(sequence, consensus_seq))

    # Join highlighted HTML into one string
    data = '</div><div class="oneliner">'.join(lines)
    data = '<div class="oneliner">' + data + '</div>'
    html_code = "<div style='padding: 4px;'>" + data + "</div>"
    return html_code


def make_alignment_matches(alignment, matches, horizontal=True):
    """Generate a MSA view with highlighted regions (matches)

    Args:
        alignment (list): List of sequences (SeqRecords)
        matches (dict): Entries with seq_id, start, end, [color]
        horizontal (bool): Orientation

    Returns:
        str: Coloured sequences (HTML)
    """
    lines = []
    for i, sequence in enumerate(alignment):
        seq_matches = []
        for site, coords in matches['matches'].iteritems():
            for coord in coords:
                if coord[0] == i:
                    seq_matches.append(coord)
        # Highlight the matching regions of the sequence
        lines.append(add_match_colors(sequence, seq_matches))

    # Join highlighted HTML into one string
    if horizontal:
        data = '</div><div class="oneliner">'.join(lines)
        data = '<div class="oneliner">' + data + '</div>'
        html_code = "<div style='padding: 4px;'>" + data + "</div>"
    else:
        data = '<br>'.join(lines)
        html_code = "<div style='padding: 4px;'>\
     <pre id='sortable' class='alignment' style='padding: 0px;'>" + data + "</pre>\
        </div>"
    return html_code


def get_possible_sequence_codon_probabilities(seqObj, G, nl2id, cu='mono'):
    """Generate height profile sequence-codon-probabilities

    Args:
        seqObj (object): FASTA-SeqRecord (BioPython) of reference sequence
        G (DiGraph): Markov-chain with current values (mono-, di-codon-usage)
        nl2id (list): Mapping from codons to node-ids (G)
        cu (str): Codon-usage selector ['mono'|'di']

    Returns:
        dict: Frequencies
    """

    freqs = {}
    # Generate codon opportunities for sequence info object
    # Di-Codon-Usage
    if cu == 'di':
        # Start with first codon of sequence
        start_codon = seqObj[0][0]
        start_aa = codonMapping['forward'][start_codon]
        syn_codons = codonMapping['backward'][start_aa]
        freqs[0] = {}
        for syn_codon in syn_codons:
            freqs[0][syn_codon] = G.node[nl2id[syn_codon]]['prob']
        # Rest of sequence
        step = 1
        for n in range(0, len(seqObj) - step, step):
            first_codon = seqObj[n][0]
            sec_codon = seqObj[n + 1][0]
            sec_aa = codonMapping['forward'][sec_codon]
            syn_codons = codonMapping['backward'][sec_aa]
            freqs[n + 1] = {}
            for syn_codon in syn_codons:
                freqs[n + 1][syn_codon] = G.get_edge_data(nl2id[first_codon], nl2id[syn_codon])['normed_transition']
    # Mono-Codon-Usage
    else:
        # Start with first codon of sequence
        start_codon = seqObj[0][0]
        start_aa = codonMapping['forward'][start_codon]
        syn_codons = codonMapping['backward'][start_aa]
        freqs[0] = {}
        for syn_codon in syn_codons:
            freqs[0][syn_codon] = G.node[nl2id[syn_codon]]['prob']
        # Rest of sequence
        step = 1
        for n in range(0, len(seqObj) - step, step):
            sec_codon = seqObj[n + 1][0]
            sec_aa = codonMapping['forward'][sec_codon]
            syn_codons = codonMapping['backward'][sec_aa]
            freqs[n + 1] = {}
            for syn_codon in syn_codons:
                freqs[n + 1][syn_codon] = G.node[nl2id[syn_codon]]['prob']

    return freqs


def make_sequence_codon_usage_profile(seqObj, G, nl2id, cu='mono', scale='relative_adapted', display='normal'):
    """Helper method to generate height profile sequence-codon-probabilities

    Args:
        seqObj (object): FASTA-SeqRecord (BioPython) of reference sequence
        G (DiGraph): Markov-chain with current values (mono-, di-codon-usage)
        nl2id (list): Mapping from codons to node-ids (G)
        cu (str): Codon-usage selector ['mono'|'di']
        scale (str): Profile scale for 'di'-mode: ['relative_adapted'|'uniform']
        display (str): Display mode: normal (one-line) with maxima or 'extended' (10-step profile)

    Returns:
        str: HTML-code of sequence codon-usage profile
    """

    # Generate codon opportunities for sequence info object
    freqs = get_possible_sequence_codon_probabilities(seqObj.seq_info, G, nl2id, cu)

    # Prepare HTML output
    r = {0: ['<tr>'], 1: ['<tr>'], 2: ['<tr>'], 3: ['<tr>'], 4: ['<tr>'], 5: ['<tr>'], 6: ['<tr>'], 7: ['<tr>'],
         8: ['<tr>'], 9: ['<tr>'], 10: ['<tr>'], 11: ['<tr>'], 12: ['<tr>']}
    nil_cell = '<td class="cell"></td>'
    min_cell = '<td class="cell min-cell"></td>'
    mid_cell = '<td class="cell mid-cell"></td>'
    max_cell = '<td class="cell max-cell"></td>'
    for index, value in seqObj.seq_info.items():
        codon = value[0]; prob = value[1];
        minVal = min(freqs[index].values())
        maxVal = max(freqs[index].values())
        minInd = prob == minVal and len(freqs[index]) > 1
        maxInd = prob == maxVal and len(freqs[index]) > 1

        # Profile: Scaling of frequencies (uniform|relative Adaptiveness)
        if scale == 'uniform':
            # Discrete | Uniform
            prop_pos = (sorted(freqs[index].values()).index(prob)+1)/float(len(freqs[index]))
            rAInd = 11-int(ceil(round(prop_pos / 0.1, 1)))
        else:
            # Relative Adaptiveness
            # prop_rA = prob/float(max(freqs[index].values()))
            prop_rA = (prob-minVal)/float(maxVal-minVal) if (prob-minVal) > 0 and (maxVal-minVal) > 0 else 0.1
            rAInd = 11-int(ceil(round(prop_rA / 0.1, 1)))
        # print rAInd

        # Switch display mode
        if display == 'normal':
            min_cell = '<td class="cell block-min-cell"></td>'
            max_cell = '<td class="cell block-max-cell"></td>'
            if minInd:
                r[1].append(min_cell)
            elif maxInd:
                r[1].append(max_cell)
            else:
                r[1].append(nil_cell)
            r[0].append('<td class="cell head-cell">%d</td>' % (index + 1))
            r[2].append('<td class="cell head-cell">%s</td>' % codon)
            r[3].append('<td class="cell foot-cell">%s</td>' % codonMapping['forward'][codon])
        elif display == 'extended':
            if rAInd == 10 or minInd:
                r[rAInd].append(min_cell)
            elif rAInd == 1 or maxInd:
                r[rAInd].append(max_cell)
            else:
                r[rAInd].append(mid_cell)
            for i in range(1, 11):
                if i != rAInd:
                    r[i].append(nil_cell)
            r[0].append('<td class="cell head-cell">%d</td>' % (index + 1))
            r[11].append('<td class="cell head-cell">%s</td>' % codon)
            r[12].append('<td class="cell foot-cell">%s</td>' % codonMapping['forward'][codon])

    data = ''
    for i, row in r.items():
        data = data + "".join(row) + '</tr>'
    html_code = "<h5>" + seqObj.description + "</h5><table style='width: inherit; min-width: 100%; table-layout: fixed;'>" + data + "</table>"

    # Return markup
    return html_code


def random_color():
    """Helper method to return a random color from a HTML color palette"""
    colors = ['#77dd88', '#99ee66', '#55bb33', '#55bb33', '#9999ff',
              '#77dd88', '#5555ff', '#66bbff', '#ffcc77', '#66bbff',
              '#66bbff', '#55bb33', '#eeaaaa', '#55bb33', '#ffcc77',
              '#ff4455', '#ff4455', '#66bbff', '#9999ff', '#9999ff']
    # Return one color
    return random.sample(colors, 1)[0]


def embed_html_into_svg(code):
    """Image helper method: Embed a given HTML string into a standalone SVG-image"""
    header = "data:image/svg+xml;utf8,"
    pre_code = "<svg xmlns='http://www.w3.org/2000/svg' width='540' height='344'>\
            <foreignObject width='100%' height='100%'>\
            <style type='text/css'>.cu_table { padding: 2px; border: 1px solid #CCC; font-family: Consolas, Courier, 'Courier New'; font-size: 14px; }</style>\
            <div xmlns='http://www.w3.org/1999/xhtml'>";
    suf_code = "</div></foreignObject></svg>";
    html_code = header + pre_code + code + suf_code
    # Return SVG-code
    return html_code
