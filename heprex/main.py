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
# Created by dsimm1 on 04/08/15.
#

__author__ = 'dsimm'

# General Inclusions
from difflib import SequenceMatcher
# Local Inclusions
from .bio_seq_encoder import *
from .bio_utilities import *
from .markov_chain import *

# Calls initial methods
def main():
    # Create Codon-Usages
    dicodons = calculate_abs_codon_usages("host_organisms/e.coli_k12_mg1655/NC_000913.3.nucleotide.fasta.txt")
    dicodonsRel = calculate_rel_dicodon_usage(dicodons)
    # Create Markov Chain
    G, nl2codon, nl2id = create_markov_chain(all_codons, dicodonsRel)

if __name__ == '__main__':
    main()

