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
from Bio.Data import CodonTable
# Load System packages
import os.path

###########################################################################
# Definitions / Constants
###########################################################################

# PATHS
BASE_APP_PATH      =  str(os.path.dirname(os.path.realpath(__file__))) + '/../'
"""str: Path to base destination of the HePrEx project folder."""
BASE_PATH_BIO      = 'bio_data/'
"""str: Subfolder to all biological source data."""
SYSTEM_TMP_PATH    = '/tmp/heprex/'
"""str: Main temporary folder."""
APP_TMP_PATH       = 'web/tmp/'
"""str: Project specific path to temporary folder."""
TMP_IMAGE_PATH     = 'web/static/tmp/img/'
"""str: Path to temporary image folder."""
WEB_TMP_IMAGE_PATH = 'static/tmp/img/'
"""str: Path to public visible image folder."""
BAIROCH_FILE_PATH  = BASE_PATH_BIO + 'local_variables/REBASE/rebase_v506_bairoch.txt'
"""str: Path to the Bairoch source file (specification of restriction sites)."""

# Postgres-Database-Connection
DB_HOST = os.getenv('DB_HOST', 'localhost')
DB_PORT = os.getenv('DB_PORT', 5433)
DB_USER = os.getenv('DB_USER', 'postgres')
DB_PASS = os.getenv('DB_PASS', 'databix')

# Caching-Server: Redis
CACHE_HOST = os.getenv('CACHE_HOST', 'localhost')
CACHE_PORT = os.getenv('CACHE_PORT', 6379)

aminoacids = {'A': 'Ala', 'C': 'Cys', 'E': 'Glu', 'D': 'Asp',
              'G': 'Gly', 'F': 'Phe', 'I': 'Ile', 'H': 'His',
              'K': 'Lys', 'M': 'Met', 'L': 'Leu', 'N': 'Asn',
              'Q': 'Gln', 'P': 'Pro', 'S': 'Ser', 'R': 'Arg',
              'T': 'Thr', 'W': 'Trp', 'V': 'Val', 'Y': 'Tyr'}
"""dict: Dictionary of Single-Letter Aminoacids mapped to 3-Letter abbreviations."""

codonMapping = { 'forward': {}, 'backward': {} }
codonMapping['forward'] = { 'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                            'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                            'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
                            'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                            'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                            'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
                            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                            'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                            'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
codonMapping['backward'] = {
                             'A': ["GCT", "GCC", "GCA", "GCG"],
                             'D': ["GAT", "GAC"],
                             'F': ["TTT", "TTC"],
                             'H': ["CAT", "CAC"],
                             'K': ["AAA", "AAG"],
                             'M': ["ATG"],                # START
                             'P': ["CCT", "CCC", "CCA", "CCG"],
                             'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                             'T': ["ACT", "ACC", "ACA", "ACG"],
                             'W': ["TGG"],
                             'C': ["TGT", "TGC"],
                             'E': ["GAA", "GAG"],
                             'G': ["GGT", "GGC", "GGA", "GGG"],
                             'I': ["ATT", "ATC", "ATA"],
                             'L': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
                             'N': ["AAT", "AAC"],
                             'Q': ["CAA", "CAG"],
                             'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                             'V': ["GTT", "GTC", "GTA", "GTG"],
                             'Y': ["TAT", "TAC"],
                             '*': ["TAA", "TGA", "TAG"]  # STOP
}
codonMapping['backward_extended'] = {
                             'A': ["GCT", "GCC", "GCA", "GCG"],
                             'B': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "AAT", "AAC"],  # R or N
                             'D': ["GAT", "GAC"],
                             'F': ["TTT", "TTC"],
                             'H': ["CAT", "CAC"],
                             'K': ["AAA", "AAG"],
                             'M': ["ATG"],                # START
                             'P': ["CCT", "CCC", "CCA", "CCG"],
                             'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                             'T': ["ACT", "ACC", "ACA", "ACG"],
                             'W': ["TGG"],
                             'C': ["TGT", "TGC"],
                             'E': ["GAA", "GAG"],
                             'G': ["GGT", "GGC", "GGA", "GGG"],
                             'I': ["ATT", "ATC", "ATA"],
                             'L': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
                             'J': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA"],  # L or I
                             'N': ["AAT", "AAC"],
                             'Q': ["CAA", "CAG"],
                             'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                             'V': ["GTT", "GTC", "GTA", "GTG"],
                             'Y': ["TAT", "TAC"],
                             'X': [],  # Unknown
                             'Z': ["GAA", "GAG", "CAA", "CAG"],  # E or Q
                             '*': ["TAA", "TGA", "TAG"]  # STOP
}
host_organisms = {  'e.coli': {
                        'genome': 'host_organisms/e.coli_k12_mg1655/NC_000913.3.nucleotide.fasta.txt',
                        'ref': 'e.coli',
                        'ref_genome': 'host_organisms/e.coli_k12_mg1655/NC_000913.3.nucleotide.fasta.txt',
                        'genbank': 'host_organisms/e.coli_k12_mg1655/NC_000913.3.gb',
                        'fullname': 'Escherichia coli str. K-12 substr. MG1655',
                        'paxdb_id': 511145,
                    },
                    's.cerevisiae': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'Saccharomyces cerevisiae S288c',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_low': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_low_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - low expressed proteins',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_low_weighted': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_low_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - low expressed proteins - weighted',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_low_weighted_invert': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_low_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - low expressed proteins - weighted, inverted',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_mid': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_mid_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - mid expressed proteins',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_mid_weighted': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_mid_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - mid expressed proteins - weighted',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_mid_weighted_invert': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_mid_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - mid expressed proteins - weighted, inverted',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_high': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_high_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - high expressed proteins',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_high_weighted': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_high_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - high expressed proteins - weighted',
                        'paxdb_id': 4932,
                    },
                    's.cerevisiae_high_weighted_invert': {
                        'genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerev_high_expr_dna_seqs.fasta.txt',
                        'ref': 's.cerevisiae',
                        'ref_genome': 'host_organisms/s.cerevisiae_paxdb_4932/s.cerevisiae.paxdb.postgres.all.exp.fasta',
                        'genbank': 'host_organisms/s.cerevisia_s288c/NC_001133.9.gb',
                        'fullname': 'S. cerevisiae - high expressed proteins - weighted, inverted',
                        'paxdb_id': 4932,
                    },
                    'a.thaliana': {
                        'genome': 'host_organisms/a.thaliana/a.thaliana_complete_genome.fasta.txt',
                        'ref': 'a.thaliana',
                        'ref_genome': 'host_organisms/a.thaliana/a.thaliana_complete_genome.fasta.txt',
                        'genbank': 'host_organisms/a.thaliana/NC_003070.9.gb',
                        'fullname': 'Arabidopsis thaliana',
                        'paxdb_id': 3702,
                    },
                    'h.sapiens': {
                        'genome': 'host_organisms/h.sapiens/rna.fa',
                        'ref': 'h.sapiens',
                        'ref_genome': 'host_organisms/h.sapiens/rna.fa',
                        'genbank': '',
                        'fullname': 'Homo Sapiens',
                        'paxdb_id': 9606,
                    },
                    'a.niger': None,
                    'c.elegans': None,
                    'd.melanogaster': None,
                  }

# Define various codon sets
standard_table    = CodonTable.unambiguous_dna_by_name["Standard"]
valid_codons      = list(standard_table.forward_table.keys())
stop_codons       = standard_table.stop_codons
all_codons        = valid_codons + stop_codons
all_codons_marked = valid_codons + ['TAA*', 'TAG*', 'TGA*']

# Responsible project leaders
responsibles = [
    {   'prename': 'Dominic',
        'surname': 'Simm',
        'title': 'M. Sc.',
        'tel_tech': '+4955139172060',
        'tel_human': '+ 49 (0)551 39 172060',
        'fax': '+ 49 (0)551 39 14693',
        'mail': 'dominic.simm@cs.uni-goettingen.de',
        'inst': 'Institute of Computer Science',
        'oa': 'Department: Theoretical Computer Science and Algorithmic Methods',
        'img': 'static/img/profile_simm.jpg'},
    {   'prename': 'Martin',
        'surname': 'Kollmar',
        'title': 'Dr.',
        'tel_tech': '+495512012260',
        'tel_human': '+ 49 (0)551 201 2260',
        'fax': '+ 49 (0)551 201 2202',
        'mail': 'mako@nmr.mpibpc.mpg.de',
        'inst': 'Max-Planck-Institute for Biophysical Chemistry',
        'oa': 'Department: NMR-based Structural Biology; Group: Systems Biology of Motor Proteins',
        'img': 'static/img/profile_kollmar.png'},
    {   'prename': 'Stefan',
        'surname': 'Waack',
        'title': 'Prof.',
        'tel_tech': '+4955139172011',
        'tel_human': '+ 49 (0)551 39 172011',
        'fax': '+ 49 (0)551 39 14693',
        'mail': 'waack@informatik.uni-goettingen.de',
        'inst': 'Institute of Computer Science',
        'oa': 'Department: Theoretical Computer Science and Algorithmic Methods',
        'img': 'static/img/profile_waack.jpg'},
]


# Make Codon dictionaries with initial values
def initialize_codons():
    """Helper function to fill codon dictionaries (codons, dicodons) with initial values

    ``{'GAG': 0, 'AGC': 0, ...}``

    Returns:
        dict, dict: The both dictionaries (codons, dicodons).
    """

    # Make Single Codon dictionary with initial values
    # 1-dimensional dict: {'GAG':0, 'AGC':0, ..}
    codons = dict((c, 0) for c in all_codons)
    # To prevent (divide by zero) side-effects, see all codons at least one time?
    # codons = dict((c, 1) for c in all_codons)

    # Make Di-Codon dictionary with initial values
    # 2-dimensional dict: {'GAG: {'GAG':0, 'AGC':0, ..}, 'AGC': {..}, ..}
    dicodons = dict((c, codons.copy()) for c in all_codons)

    return codons, dicodons


# Call on first inclusion
codons_rel, dicodons = initialize_codons()
