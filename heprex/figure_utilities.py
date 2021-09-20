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
# Load BioPython packages
from Bio import Entrez, SeqIO

# Switching between X and no-X environments
import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
# Import the pyplot interface to matplotlib:
from matplotlib import pyplot as plt

# Load MatplotLib package
from matplotlib.pylab import *
from matplotlib import gridspec
# Import Collections
from collections import OrderedDict

# Local Inclusions
from .globals import *
from .help_utilities import *


def make_codon_usage_table(host, codons, title=''):
    """Prepare and store an overview plot from observed codon-usages

    Args:
        host (str): Name of host organism
        codons (dict): Mono Codon usage (RCU)
        title (str): Title for table

    Returns:
        str: HTML code for Codon usage table
    """
    # Make vector (of values) of Codon-Usage and sort values by amino-acids [codonMapping]
    codonsSorted = OrderedDict(sorted(codons.items(), key=lambda x: codonMapping['forward'][x[0]]))

    # Relative values - Set interval [0.0; 1.0]
    codonsRelProbVec = np.divide(list(codons.values()), float(np.sum(codons.values())))
    rel_codons = dict(zip(codons.keys(), codonsRelProbVec))

    cu = \
"<table class='cu_table'>"\
"<tr><th colspan='4'>" + title + "</th></tr>"\
"<tr><td>TTT F ( %(TTT)03.3f )</td><td>TCT S ( %(TCT)03.3f )</td><td>TAT Y ( %(TAT)03.3f )</td><td>TGT C ( %(TGT)03.3f )</td></tr>"\
"<tr><td>TTC F ( %(TTC)03.3f )</td><td>TCC S ( %(TCC)03.3f )</td><td>TAC Y ( %(TAC)03.3f )</td><td>TGC C ( %(TGC)03.3f )</td></tr>"\
"<tr><td>TTA L ( %(TTA)03.3f )</td><td>TCA S ( %(TCA)03.3f )</td><td>TAA * ( %(TAA)03.3f )</td><td>TGA * ( %(TGA)03.3f )</td></tr>"\
"<tr><td>TTG L ( %(TTG)03.3f )</td><td>TCG S ( %(TCG)03.3f )</td><td>TAG * ( %(TAG)03.3f )</td><td>TGG W ( %(TGG)03.3f )</td></tr>"\
"<tr><td>CTT L ( %(CTT)03.3f )</td><td>CCT P ( %(CCT)03.3f )</td><td>CAT H ( %(CAT)03.3f )</td><td>CGT R ( %(CGT)03.3f )</td></tr>"\
"<tr><td>CTC L ( %(CTC)03.3f )</td><td>CCC P ( %(CCC)03.3f )</td><td>CAC H ( %(CAC)03.3f )</td><td>CGC R ( %(CGC)03.3f )</td></tr>"\
"<tr><td>CTA L ( %(CTA)03.3f )</td><td>CCA P ( %(CCA)03.3f )</td><td>CAA Q ( %(CAA)03.3f )</td><td>CGA R ( %(CGA)03.3f )</td></tr>"\
"<tr><td>CTG L ( %(CTG)03.3f )</td><td>CCG P ( %(CCG)03.3f )</td><td>CAG Q ( %(CAG)03.3f )</td><td>CGG R ( %(CGG)03.3f )</td></tr>"\
"<tr><td>ATT I ( %(ATT)03.3f )</td><td>ACT T ( %(ACT)03.3f )</td><td>AAT N ( %(AAT)03.3f )</td><td>AGT S ( %(AGT)03.3f )</td></tr>"\
"<tr><td>ATC I ( %(ATC)03.3f )</td><td>ACC T ( %(ACC)03.3f )</td><td>AAC N ( %(AAC)03.3f )</td><td>AGC S ( %(AGC)03.3f )</td></tr>"\
"<tr><td>ATA I ( %(ATA)03.3f )</td><td>ACA T ( %(ACA)03.3f )</td><td>AAA K ( %(AAA)03.3f )</td><td>AGA R ( %(AGA)03.3f )</td></tr>"\
"<tr><td>ATG M ( %(ATG)03.3f )</td><td>ACG T ( %(ACG)03.3f )</td><td>AAG K ( %(AAG)03.3f )</td><td>AGG R ( %(AGG)03.3f )</td></tr>"\
"<tr><td>GTT V ( %(GTT)03.3f )</td><td>GCT A ( %(GCT)03.3f )</td><td>GAT D ( %(GAT)03.3f )</td><td>GGT G ( %(GGT)03.3f )</td></tr>"\
"<tr><td>GTC V ( %(GTC)03.3f )</td><td>GCC A ( %(GCC)03.3f )</td><td>GAC D ( %(GAC)03.3f )</td><td>GGC G ( %(GGC)03.3f )</td></tr>"\
"<tr><td>GTA V ( %(GTA)03.3f )</td><td>GCA A ( %(GCA)03.3f )</td><td>GAA E ( %(GAA)03.3f )</td><td>GGA G ( %(GGA)03.3f )</td></tr>"\
"<tr><td>GTG V ( %(GTG)03.3f )</td><td>GCG A ( %(GCG)03.3f )</td><td>GAG E ( %(GAG)03.3f )</td><td>GGG G ( %(GGG)03.3f )</td></tr>"\
"</table>" % rel_codons

    return cu


def make_codon_usage_plots(host, codons, dicodons, title=''):
    """Prepare and store an overview plot from observed codon-usages

    Args:
        host (str): Name of host organism
        codons (dict): Mono Codon usage (RCU)
        dicodons (dict): Di-Codon usage (RCU)
        title (str): Title for table

    Returns:
        str: Path to rendered Codon usage plots
    """

    labl_font = {'family' : 'normal',     # normal, sans-serif, monospace
                 'weight' : 'normal',     # normal, bold, ...
                 'size'   :  8 }
    axis_font = {'family' : 'monospace',  # normal, sans-serif, monospace
                 'weight' : 'normal',     # normal, bold, ...
                 'size'   :  8 }

    # Fix the defaut saveconfig dpi
    matplotlib.rcParams['savefig.dpi'] = 160
    matplotlib.rcParams.update({'font.size': 8})

    # Make matrix (of values) of Dicodon-Usage
    # Sort 2-dimensional array [sorted by amino-acids]
    codons_sorted_keys = sorted(codons.keys(), key=lambda x: codonMapping['forward'][x])
    codons_sorted_values = []; dicodonsMat = []
    # Garantuee same order of keys and values
    for c in codons_sorted_keys:
        codons_sorted_values.append(codons[c])
        dicodonsMat.append(list(dicodons[c][d] for d in codons_sorted_keys))
    # Relative values - Set interval [0.0; 1.0]
    codonsVec = np.divide(list(codons_sorted_values), float(np.amax(codons_sorted_values)))
    # Hack: case matshow can't plot vectors
    duoCodonsVec = [codonsVec] + [codonsVec]
    # codonsVecNorm =  np.log(codonsVec)

    # Labels
    # alpha = sorted(all_codons_marked) # ['AAA', 'ACA', ..., 'TTT']
    alpha = sorted(all_codons, key=lambda x: codonMapping['forward'][x])  # ['AAA', 'ACA', ..., 'TTT']

    # Plotting settings
    fig = plt.figure(figsize=(8,9), dpi=240)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])

    # Codon-Usage
    ax0 = fig.add_subplot(gs[0])
    cax = ax0.matshow(duoCodonsVec, interpolation='nearest')
    fig.colorbar(cax, orientation="horizontal", pad=0.05, shrink=2.0)
    ax0.set_xticks(range(0,64,1))
    ax0.set_yticks(range(0,2,1))
    ax0.set_xticklabels(alpha, rotation=70, **axis_font)  #, rotation_mode="anchor")
    ax0.set_yticklabels(['Codon-Usage'], **axis_font)

    # DiCodon-Usage
    ax1 = fig.add_subplot(gs[1])
    ax1.matshow(dicodonsMat, interpolation='nearest')
    ax1.set_xticks(range(0,64,1))
    ax1.set_yticks(range(0,64,1))
    ax1.set_xticklabels(alpha, rotation=70, **axis_font)  #, rotation_mode="anchor")
    ax1.set_yticklabels(alpha, **axis_font)
    ax1.set_ylabel('First codon', fontsize=12, color='#000000')
    ax1.set_xlabel('Second codon', fontsize=12, color='#000000')
    ax1.xaxis.set_label_coords(0.5, 1.085)

    # Figure description
    if len(title) == 0:
        record = SeqIO.read(BASE_PATH_BIO + host_organisms[host]['genbank'], "genbank")
        annot = record.annotations
        title = annot['organism'] + " v." + str(annot['sequence_version']) + " gi:" + annot['gi'] + "\n" + str(annot['keywords']) + " " + str(annot["references"][0].journal)
    fig.suptitle(title, fontsize=12)

    # Optimize figure arrangement
    fig.set_tight_layout(True)
    # gs.tight_layout(fig) # iPython version for interactive plotting

    # Store in tempfile
    base_filename = "codon_usages_colormap_" + host + "_" + uniq_id() + ".png"
    plt.savefig(TMP_IMAGE_PATH + base_filename)
    return WEB_TMP_IMAGE_PATH + base_filename


def make_histograms(host, codons):
    """Prepare and store an histogram plot from observed codon-usages

    Args:
        host (str): Name of host organism
        codons (dict): Mono Codon usage (RCU)

    Returns:
        str: Path to rendered histogram plot
    """

    # Definition of Matplotlib-fonts:
    axis_font = {'family' : 'monospace',  # normal, sans-serif, monospace
                 'weight' : 'normal',     # normal, bold, ...
                 'size'   :  4 }

    # Fix the defaut saveconfig dpi
    matplotlib.rcParams['savefig.dpi'] = 160
    matplotlib.rcParams.update({'font.size': 4})

    # Make vector (of values) of Codon-Usage and sort values
    codonsSorted = OrderedDict(sorted(codons.items(), key=lambda x: x[0]))
    codonsVec = list(codonsSorted.values())

    # Labels
    alpha = sorted(all_codons_marked)  # ['AAA', 'ACA', ..., 'TTT']
    # Presets for the figure
    ind = np.arange(64)  # the x locations for the groups
    width = 0.5          # the width of the bars

    # Plotting settings
    fig = plt.figure(figsize=(3, 3), dpi=160)
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])

    # Bar Chart from Codon-Usage
    ax1 = fig.add_subplot(gs[0])
    ax1.bar(ind, np.multiply(codonsVec, float(0.75)), width, color='#000000')
    ax1.set_xticks(ind+width)
    ax1.set_xticklabels(alpha, rotation=90, **axis_font)
    # ax1.set_yticks(range(0,90000,10000))
    plt.margins(0.02, 0)

    # Bar Chart from Codon-Usage
    ax1 = fig.add_subplot(gs[1])
    ax1.bar(ind, np.multiply(codonsVec, float(1.00)), width, color='#000000')
    ax1.set_xticks(ind+width)
    ax1.set_xticklabels(alpha, rotation=90, **axis_font)
    # ax1.set_yticks(range(0,90000,10000))
    plt.margins(0.02, 0)

    # Bar Chart from Codon-Usage
    ax1 = fig.add_subplot(gs[2])
    ax1.bar(ind, np.multiply(codonsVec, float(1.25)), width, color='#000000')
    ax1.set_xticks(ind+width)
    ax1.set_xticklabels(alpha, rotation=90, **axis_font)
    # ax1.set_yticks(range(0,90000,10000))
    plt.margins(0.02, 0)

    # Optimize figure arrangement
    # plt.axis('tight') # No margins in plot
    fig.set_tight_layout(True)

    # Store in tempfile
    base_filename = "codon_usages_hist_" + host + "_" + uniq_id() + ".png"
    plt.savefig(TMP_IMAGE_PATH + base_filename)
    return WEB_TMP_IMAGE_PATH + base_filename


def make_histogram(host, data, lower_bound, upper_bound):
    """Prepare and store an ranged histogram plot from observed codon-usages

    Args:
        host (str): Name of host organism
        data (dict): Data such as Mono Codon usage (RCU)
        lower_bound (float): lower bound of range
        upper_bound (float): upper bound of range

    Returns:
        str: Path to rendered histogram plot
    """

    # Definition of Matplotlib-fonts:
    axis_font = {'family' : 'monospace', # normal, sans-serif, monospace
                 'weight' : 'normal',    # normal, bold, ...
                 'size'   :  4 }

    # Fix the defaut saveconfig dpi
    matplotlib.rcParams['savefig.dpi'] = 480
    matplotlib.rcParams.update({'font.size': 4})

    # Plotting settings
    fig = plt.figure(figsize=(3, 1), dpi=480)
    gs = gridspec.GridSpec(1, 1)

    # Bar Chart from Codon-Usage
    ax1 = fig.add_subplot(gs[0])
    # Sort and Split into 3 partitions
    np.msort(data)
    low_ind, up_ind = data.index(lower_bound), data.index(upper_bound)
    splitted_data = np.split(data, [low_ind, up_ind])
    # Make values scalable by log10(): Need to start with values > 1.0 to prevent negative values
    data = [np.log10(el) for el in np.add(splitted_data, 0.01)]
    # n, bins, patches = ax1.hist(data, 100, facecolor='green', alpha=0.5)
    n, bins, patches = ax1.hist(data, 100)
    for patch in patches[0]: patch.set_facecolor('lightblue'); patch.set_alpha(0.33);
    for patch in patches[1]: patch.set_facecolor('lightgreen')
    for patch in patches[2]: patch.set_facecolor('crimson'); patch.set_alpha(0.33);

    # Configurations
    x_ticks = range(-2, 6, 2)
    labels = [0.01, 1, 100, 10000]
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(labels, rotation=0.0, **axis_font)
    # ax1.set_xscale('symlog')
    ax1.set_xlabel('Abundance', fontsize=8.0, **axis_font)
    ax1.set_ylabel('Sequences', fontsize=6.0, **axis_font)
    plt.margins(0.02, 0.05)

    # Optimize figure arrangement
    # plt.axis('tight') # No margins in plot
    fig.set_tight_layout(True)
    # Store in tempfile
    base_filename = "hist_" + host + "_" + uniq_id() + ".png"
    plt.savefig(TMP_IMAGE_PATH + base_filename)
    return WEB_TMP_IMAGE_PATH + base_filename


# Image helper method: Embed a given HTML string into a standalone SVG-image data string
def embed_html_into_svg(code):
    """Image helper method: Embed a given HTML string into a standalone SVG-image"""
    # header = "data:image/svg+xml;utf8,"
    pre_code = "<svg xmlns='http://www.w3.org/2000/svg' width='540' height='344'>\
            <foreignObject width='100%' height='100%'>\
            <style type='text/css'>.cu_table { background-color: #FFF; padding: 2px; border: 1px solid #CCC; font-family: Consolas, Courier, 'Courier New'; font-size: 14px; }</style>\
            <div xmlns='http://www.w3.org/1999/xhtml'>"
    suf_code = "</div></foreignObject></svg>"
    # html_code = header + pre_code + code + suf_code
    html_code = pre_code + code + suf_code
    # Return SVG-code
    return html_code


# Image helper method: Embed a given HTML string into a standalone SVG-image file
def embed_html_into_svg_file(code):
    """Image helper method: Embed a given HTML string into a standalone SVG-image"""
    pre_code = "<svg xmlns='http://www.w3.org/2000/svg' width='540' height='344'>\
            <foreignObject width='540' height='344'>\
            <style type='text/css'>.cu_table { background-color: #FFF; padding: 2px; margin: 4px; border: 1px solid #CCC; font-family: Consolas, Courier, 'Courier New'; font-size: 14px; }</style>\
            <div xmlns='http://www.w3.org/1999/xhtml'>"
    suf_code = "</div></foreignObject></svg>"
    html_code = pre_code + code + suf_code

    # Store in tempfile
    base_filename = "svg_" + uniq_id() + ".svg"
    in_w = open(TMP_IMAGE_PATH + base_filename, 'w')
    in_w.write(html_code)
    in_w.close()

    return WEB_TMP_IMAGE_PATH + base_filename
