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

# Import Flask packages
from flask import Flask, render_template, request
from werkzeug.contrib.cache import RedisCache

# Import Standard python packages
import base64, json
from collections import OrderedDict

# Import HePrEx package
import heprex
from heprex.bio_utilities import *
from heprex.help_utilities import *

########################################################################################################################

# Define routes
app = Flask(__name__)


@app.route("/")
def index():
    # Parse enzymes from REBASE source
    enzymes_arr = cache.get('setup')['restriction']
    enzymes = heprex.bairoch_filter_identical_enzymes(enzymes_arr, True)

    return render_template("_heprex_gui.html", enzymes=enzymes, site='index')


@app.route("/contact")
def contact():
    return render_template("_contact.html", responsibles=responsibles, site='contact')


@app.route("/documentation")
def documentation():
    return render_template("_documentation.html", site='doc')


@app.route("/faqs")
def faqs():
    return render_template("_faqs.html", site='faqs')


@app.route("/guide")
def guide():
    return render_template("_guide.html", site='guide')


@app.route('/get_codon_usage_data', methods=['POST', 'GET'])
def get_codon_usage_data():
    # Get [base-encoded] form-data (host_organism)
    organism = str(request.args.get('organism'))
    # Get prepared setup
    # From Cache
    setup = cache.get('setup')[organism]
    # # From Backup-file
    # setup_dump = heprex.APP_TMP_PATH + "setup.dat"
    # if os.path.isfile(setup_dump):
    #     setup = pickle.load(open(setup_dump, "r"))[host_organism]

    # Make histograms
    plot_path_hist = heprex.make_histograms(organism, setup['codons'])

    # Request PaxDB data
    protein_abundances, keywords = heprex.select_paxdb_information(organism)

    return render_template("_codon_usage_mod.html", organism=organism, protein_abundances=protein_abundances,
                           keywords=keywords, plot_path=plot_path_hist)


@app.route('/get_restriction_data', methods=['POST', 'GET'])
def get_restriction_data():
    # Get [base-encoded] form-data (host_organism)
    commercial = True if request.args.get('commercial') == 'true' else False
    # Parse enzymes from REBASE source
    enzymes_arr = cache.get('setup')['restriction']
    enzymes = heprex.bairoch_filter_identical_enzymes(enzymes_arr, commercial)

    return render_template("_restriction_enzymes_mod.html", enzymes=enzymes)


@app.route('/get_protein_abundance', methods=['POST', 'GET'])
def get_protein_abundance():
    # Get [base-encoded] form-data (host_organism)
    organism = str(request.args.get('organism'))
    lower_bound = int(request.args.get('min'))
    upper_bound = int(request.args.get('max'))
    # Filter enzymes from REBASE source
    resultDict = heprex.get_paxdb_portion(organism, lower_bound, upper_bound)
    entries = heprex.get_paxdb_data(organism)
    resultDict['img_path'] = heprex.make_histogram(organism, entries, resultDict['log_min'], resultDict['log_max'])
    resultDict['log_min'] = '{:.1f}'.format(resultDict['log_min'])
    resultDict['log_max'] = '{:.1f}'.format(resultDict['log_max'])
    print(resultDict)

    return json.dumps(resultDict)


@app.route('/get_dna_translation', methods=['POST', 'GET'])
def get_dna_translation():
    # Get [base64-encoded] form-data (search-term, host_organism, etc.)
    base64_code = request.args.get('search_protein')
    search_word = unicode(base64.decodestring(base64_code))

    # Parse uploaded sequence into BioSeq-Object
    if len(search_word) > 0:
        fastaObj = heprex.load_fasta_data(search_word)
    else:
        return '{"error": "empty_sequence"}'
    # Determine type of uploaded sequence (dna/protein)
    # ... and convert if necessary into protein sequence
    bioObj = fastaObj[0]; bioSeq = ''
    try:
        bioObj.seq.alphabet = IUPAC.unambiguous_dna
        bioObj.seq = bioObj.seq.translate()
        bioSeq = str(bioObj.seq)
    except CodonTable.TranslationError:
        bioObj.seq.alphabet = IUPAC.protein
        # TODO: Identify Extended_IUPAC
        # If some of the following abbrev. (B,X,Z,J,U,O) are found
        if bioObj.seq.count('B') > 0 or bioObj.seq.count('X') > 0 or bioObj.seq.count('Z') > 0 or \
                bioObj.seq.count('J') > 0 or bioObj.seq.count('U') > 0 or bioObj.seq.count('O') > 0:
            return '{"error": "extended_iupac"}'

    return '{"trans_seq": "' + bioSeq + '"}'


@app.route('/get_prediction_data', methods=['POST', 'GET'])
def get_prediction_data():

    # Get [base64-encoded] form-data (search-term, host_organism, etc.)
    host_organism = str(request.args.get('host_organism'))
    # Set Standard Host organism
    if host_organism not in heprex.host_organisms.keys() or heprex.host_organisms[host_organism] is None:
        host_organism = 'e.coli'
    number = int(request.args.get('seq_number'))
    # Decode safe transmitting encoding (Base64)
    base64_code = request.args.get('search_protein')
    search_word = unicode(base64.decodestring(base64_code))
    # Get restriction sites for search as array
    restriction_site_ids = json.loads(request.args.get('restriction_sites'))
    most_abundant = int(json.loads(request.args.get('cu_most_abundant_proteins')))
    gb_keywords = json.loads(request.args.get('cu_gb_keywords'))
    abundance_range = json.loads(request.args.get('cu_abundance_range'))
    rnafold = True if request.args.get('rnafold_switch') == 'on' else False

    # Parse uploaded sequence into BioSeq-Object
    if len(search_word) > 0:
        fastaObj = heprex.load_fasta_data(search_word)
    else:
        return '{"error": "empty_sequence"}'
    # Determine type of uploaded sequence (DNA/Protein) declare type
    bioObj = fastaObj[0]; bioDNAObj = None; additional_data = ''
    try:
        bioObj.seq.alphabet = IUPAC.unambiguous_dna
        bioDNAObj = deepcopy(bioObj)
        bioObj.seq = bioObj.seq.translate()
        additional_data = '"trans_seq": "' + str(bioObj.seq) + '"'
    except CodonTable.TranslationError:
        bioDNAObj = None
        bioObj.seq.alphabet = IUPAC.protein
        # TODO: Identify Extended_IUPAC
        # If some of the following abbrev. (B,X,Z,J,U,O) are found
        if bioObj.seq.count('B') > 0 or bioObj.seq.count('X') > 0 or bioObj.seq.count('Z') > 0 or \
                bioObj.seq.count('J') > 0 or bioObj.seq.count('U') > 0 or bioObj.seq.count('O') > 0:
            return '{"error": "extended_iupac"}'

    # Control output
    print(fastaObj)
    print("Seq: %s" % bioObj.seq)
    print("Type: %s\n" % repr(bioObj.seq.alphabet))

    # Passed FASTA sequence has been successful parsed
    if len(fastaObj) > 0 and len(bioObj.seq) > 0:
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        heprex.store_fasta_data(bioObj, tmp_file.name)

        # Get prepared HePrEx setup from Cache
        setup = cache.get('setup')

        # Generate custom codon usages
        # Case 1: Most abundant proteins | '-1': code for all sequences
        if most_abundant > 0 and most_abundant != -1:
            # Get prepared HePrEx setup from Cache
            try:
                setup_custom = cache.get('setup')['model_%s_%d' % (host_organism, most_abundant)]
            except KeyError:
                print('Abundance list:')
                paxdb_seqs = heprex.paxdb_select_organism_seqs_by_number(host_organism, most_abundant)
                codons_custom, dicodons_custom = heprex.calculate_abs_codon_usages_string(paxdb_seqs)
                # Quality check of codons
                if not heprex.check_codons(codons_custom, dicodons_custom):
                    return '{"error": "The selected sequence information is not sufficient to build a reliable model of the markov chains."}'
                setup_custom = heprex.recalculate_requirements(codons_custom, dicodons_custom)
                setup['model_%s_%d' % (host_organism, most_abundant)] = setup_custom
            # Create heterologous sequences
            seqs = heprex.create_heterologous_sequences(tmp_file.name, number, setup_custom)
            codons_custom = setup_custom['codons']; dicodons_custom = setup_custom['dicodons']
            codons_suggest, dicodons_suggest = heprex.calculate_abs_codon_usages_string(seqs)

        # Case 2: Genbank keywords
        elif isinstance(gb_keywords, list) and len(gb_keywords) > 0:
            try:
                # Get prepared HePrEx setup from Cache
                setup_custom = cache.get('setup')['model_%s_%s' % (host_organism, ','.join(gb_keywords))]
            except KeyError:
                # Calculate custom codon-usage
                paxdb_seqs = heprex.paxdb_select_organism_seqs_by_keywords(host_organism, gb_keywords)
                codons_custom, dicodons_custom = heprex.calculate_abs_codon_usages_string(paxdb_seqs)
                if not heprex.check_codons(codons_custom, dicodons_custom):
                    return '{"error": "The selected sequence information is not sufficient to build a reliable model of the markov chains."}'
                setup_custom = heprex.recalculate_requirements(codons_custom, dicodons_custom)
                setup['model_%s_%s' % (host_organism, ','.join(gb_keywords))] = setup_custom
            # Create heterologous sequences
            print('Keywords:')
            seqs = heprex.create_heterologous_sequences(tmp_file.name, number, setup_custom)
            codons_custom = setup_custom['codons']; dicodons_custom = setup_custom['dicodons']
            codons_suggest, dicodons_suggest = heprex.calculate_abs_codon_usages_string(seqs)

        # Case 3: Abundance range (PaxDB)
        elif isinstance(abundance_range, list) and len(abundance_range) == 2 and \
                (int(abundance_range[0]) >= 0 and int(abundance_range[1]) < 100 or int(abundance_range[0]) > 0 and int(abundance_range[1]) <= 100):
            min = int(abundance_range[0]); max = int(abundance_range[1])
            try:
                # Get prepared HePrEx setup from Cache
                setup_custom = cache.get('setup')['model_%s_%d-%d' % (host_organism, min, max)]
            except KeyError:
                # Calculate custom codon-usage
                paxdb_seqs = heprex.paxdb_select_organism_seqs_by_range(host_organism, min, max)
                codons_custom, dicodons_custom = heprex.calculate_abs_codon_usages_string(paxdb_seqs)
                if not heprex.check_codons(codons_custom, dicodons_custom):
                    return '{"error": "The selected sequence information is not sufficient to build a reliable model of the markov chains."}'
                setup_custom = heprex.recalculate_requirements(codons_custom, dicodons_custom)
                setup['model_%s_%d-%d' % (host_organism, min, max)] = setup_custom
            # Create heterologous sequences
            print('Abundance Range:')
            seqs = heprex.create_heterologous_sequences(tmp_file.name, number, setup_custom)
            codons_custom = setup_custom['codons']; dicodons_custom = setup_custom['dicodons']
            codons_suggest, dicodons_suggest = heprex.calculate_abs_codon_usages_string(seqs)

        # Case 4: Precalculated whole Genome codon usages
        else:
            # Get prepared HePrEx setup from Cache
            setup_custom = cache.get('setup')[host_organism]
            setup_custom['restriction'] = cache.get('setup')['restriction']
            # Generate predictions for user input
            print('General:')
            seqs = heprex.create_heterologous_sequences(tmp_file.name, number, setup_custom)
            codons_custom = False
            codons_suggest, dicodons_suggest = heprex.calculate_abs_codon_usages_string(seqs)

        # Store data in Redis and Pickle (disk)
        cache.set('setup', setup, timeout=30*60*60*24)  # One month (in seconds);
        with open(setup_dump, "w+") as f:
            pickle.dump(setup, f)
        # Get prepared HePrEx setup from Cache
        setup = cache.get('setup')[host_organism]

        #####################################################################################
        # Preparation of heterologous genes
        #####################################################################################

        # Search for specified restriction sites
        matches = {}
        if restriction_site_ids:
            enzymes_arr = cache.get('setup')['restriction']
            restriction_sites = heprex.bairoch_select_restriction_enzymes(enzymes_arr, restriction_site_ids)
            matches = heprex.identify_restriction_sites(seqs, restriction_sites)
            print("\n#############################################\nRestriction Site matches:")
            print(matches)

        # Prepare return data
        msa_highlight_aa = heprex.make_alignment_aa(seqs)
        msa_highlight_diff = heprex.make_alignment_diffs(seqs)
        if restriction_site_ids:
            msa_highlight_matches = heprex.make_alignment_matches(seqs, matches)
        else:
            msa_highlight_matches = ''

        # Generate codon-usage plots - Color-Maps and Tables
        figures = {}; cu_table = {}
        rel_codons = heprex.get_rel_codon_usage(setup['G'], setup['nl2codon'])
        rel_dicodons = heprex.get_rel_dicodon_usage(setup['G'], setup['nl2codon'])
        # Display codon-usages of Markov model (normal, weighted, inverted)
        figures['colormap_genome'] = heprex.make_codon_usage_plots(host_organism, rel_codons, rel_dicodons)
        cu_table['genome'] = heprex.embed_html_into_svg_file(heprex.make_codon_usage_table(host_organism, rel_codons, heprex.host_organisms[host_organism]['fullname']))
        figures['colormap_suggest'] = heprex.make_codon_usage_plots(host_organism, codons_suggest, dicodons_suggest, 'Codon usage of suggestions')
        cu_table['suggest'] = heprex.embed_html_into_svg_file(heprex.make_codon_usage_table(host_organism, codons_suggest, 'Codon usage of suggestions'))
        if codons_custom:
            figures['colormap_custom'] = heprex.make_codon_usage_plots(host_organism, codons_custom, dicodons_custom, 'Codon usage of customized protein set')
            cu_table['custom'] = heprex.embed_html_into_svg_file(heprex.make_codon_usage_table(host_organism, codons_custom, 'Codon usage of customized protein set'))
        else:
            figures['colormap_custom'] = ''
            cu_table['custom'] = ''

        # Single sequence information - Stats, RSCU
        print("\n#############################################\nSingle sequence information calculation:")
        cu_table_single = {}; stats_single = {}
        for idx, seqObj in enumerate(seqs):
            codons_single, dicodons_single = heprex.calculate_abs_codon_usages_string([seqObj])
            cu_table_single[idx] = heprex.embed_html_into_svg_file(heprex.make_codon_usage_table(host_organism, codons_single, "Codon usage of sequence " + str(idx + 1)))
            stats_single[idx] = heprex.seq_statistics(seqObj, setup['codons'], rnafold)

        # Render data into template
        # Sequences - Alignment section
        div_alignments = render_template("_result_alignments.html", seqs=seqs, hha=msa_highlight_aa, hhd=msa_highlight_diff, hhm=msa_highlight_matches)
        # Sequences - Single view
        div_sequences = render_template("_result_sequence_singles.html", seqs=seqs, cu_table=cu_table_single, stats=stats_single, rnafold=rnafold)
        # Codon-Usages
        div_figures = render_template("_result_figures.html", figures=figures, cu_table=cu_table)

        # Return pages
        return div_alignments + div_sequences + div_figures

    else:
        return '{"error": "error"}'


if __name__ == "__main__":
    # Initialize memcache
    cache = RedisCache(host=heprex.CACHE_HOST, port=heprex.CACHE_PORT, password=None, db=0, default_timeout=300)
    # Setup tmp-path environment
    if not os.path.exists(heprex.WEB_TMP_IMAGE_PATH):
        os.makedirs(heprex.WEB_TMP_IMAGE_PATH)
    # Setup heprex environment
    setup_dump = heprex.APP_TMP_PATH + "setup.dat"

    # Force reload of profiles into cache
    setup = cache.get('setup')
    cache.delete('setup')
    # # For re-creation of profiles and write new 'setup' dump-file
    # if os.path.isfile(setup_dump): os.unlink(setup_dump)

    # Profiles not stored in Cache
    if cache.get('setup') is None:
        if os.path.isfile(setup_dump):
            print('Restore cache from .dat-file:')
            setup = pickle.load(open(setup_dump, "r"))

            # # Parse enzymes from REBASE source
            # # SLOWs web-app start down! Can not be commented out on Global Restart
            # setup['restriction'] = heprex.parse_bairoch_file(BAIROCH_FILE_PATH)

            # Modify or add profiles
            profiles = OrderedDict([
               # ('s.cerevisiae_high', host_organisms['s.cerevisiae_high']['genome']),
               # ('s.cerevisiae_mid', host_organisms['s.cerevisiae_mid']['genome']),
               # ('s.cerevisiae_low', host_organisms['s.cerevisiae_low']['genome']),
               # ('s.cerevisiae_high_weighted', host_organisms['s.cerevisiae_high']['genome']),
               # ('s.cerevisiae_mid_weighted', host_organisms['s.cerevisiae_mid']['genome']),
               # ('s.cerevisiae_low_weighted', host_organisms['s.cerevisiae_low']['genome']),
               # ('s.cerevisiae_high_weighted_invert', host_organisms['s.cerevisiae_high']['genome']),
               # ('s.cerevisiae_mid_weighted_invert', host_organisms['s.cerevisiae_mid']['genome']),
               # ('s.cerevisiae_low_weighted_invert', host_organisms['s.cerevisiae_low']['genome']),
            ])
            for profile, file_path in profiles.items():
                print('Recreating: %s' % profile)
                if profile.find('weight') == -1:
                    codons_acu, dicodons_acu = heprex.calculate_abs_codon_usages(BASE_PATH_BIO + file_path)
                else:
                    codons_acu, dicodons_acu = heprex.calculate_abs_codon_usages_weighted(BASE_PATH_BIO + file_path)
                if profile.find('invert') == -1:
                    setup_custom = heprex.recalculate_requirements(
                        codons_acu, dicodons_acu,
                        setup[host_organisms[profile]['ref']]['codons'], setup[host_organisms[profile]['ref']]['dicodons']
                    )
                else:
                    setup_custom = heprex.recalculate_requirements_invert(
                        codons_acu, dicodons_acu,
                        setup[host_organisms[profile]['ref']]['codons'], setup[host_organisms[profile]['ref']]['dicodons']
                    )
                # Store into cache
                setup[profile] = setup_custom
                # Set creation date
                try:
                    setup.pop(setup['timestamp'], None)
                    setup['timestamp'] = heprex.timestamp()
                    setup[setup['timestamp']] = ''
                except KeyError:
                    setup['timestamp'] = heprex.timestamp()

            # Store data in Redis
            cache.set('setup', setup, timeout=30*60*60*24)  # One month (in seconds);
            print('Profiles updated in Cache: ' + repr(setup.keys()))
            # # Store pre-computed data in Pickle
            # with open(setup_dump, "w+") as f:
            #     pickle.dump(setup, f)

        else:
            print('Recreate profiles from scratch:')
            setup = heprex.prepare_requirements()

            # Parse enzymes from REBASE source
            # SLOWs web-app start down! Can not be commented out on Global Restart
            enzymes_arr = heprex.parse_bairoch_file(BAIROCH_FILE_PATH)
            setup['restriction'] = enzymes_arr

            # Store data in Redis and Pickle (disk)
            cache.set('setup', setup, timeout=30*60*60*24)  # One month (in seconds);
            print('Profiles written to Cache: ' + repr(setup.keys()))
            with open(setup_dump, "w+") as f:
                pickle.dump(setup, f)

    # Profiles stored in Cache
    else:
        setup = cache.get('setup')
        print('Found cached profiles: ' + repr(setup.keys()))
        # # Helper to store to Cache to file
        # with open(setup_dump, "w+") as f:
        #     pickle.dump(setup, f)

    # Start flask web-application
    print("Start up Flask web-server on port 2020 ...")
    app.run(host='0.0.0.0', port=2020, debug=True)
