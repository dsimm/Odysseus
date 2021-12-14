/**
 * Created by dsimm on 21/10/15.
 */

// Constants
var modal_id = '#info_modal_lg';

// Copy selected codon usage in hidden form field
function select_host_organism(host_code, host_text) {
    $('input[id="host_organism"]').val(host_code);
    $('button[id="button_codon_usage"]').html(host_text + ' <span class="caret"></span>');
    $('span[id="active_organism"]').html(host_text);
    // Get some data
    request_codon_usage_data(host_code);
    $('div[id="codon_usage_mod"]').toggle(true);
}

function change_codon_usage_state(id) {
    $('label[for^=cu_scaling_]').removeClass('active');
    $('label[for=' + id + ']').addClass('active');
}

// Copy selected codon usage in hidden form field
function set_seq_number(seq_number) {
    $('input[id="seq_number"]').val(seq_number);
    $('button[id="button_number"]').html(seq_number + ' <span class="caret"></span>');
}

// Request codon usage data from server
function request_codon_usage_data(organism) {
    $.ajax({
        url: "/get_codon_usage_data",
        data: {
            organism: organism,
        },
        beforeSend: function () {
            $("#loading_modal").modal('show');
        }
    }).done(function (data) {
        $("#loading_modal").modal('hide');
        $('div[id="codon_usage_mod"]').html(data);
        // Activate JS elements
        activate_codon_usage_elements(organism);
    });
}

// Request restriction data from server
function request_restriction_data(commercial) {
    $.ajax({
        url: "/get_restriction_data",
        data: {
            commercial: commercial,
        }
    }).done(function (data) {
        $('div[id="restriction_enzymes_content"]').html(data);
        // Activate JS elements
        activate_restriction_elements();
    });
}

// Request restriction data from server
function request_protein_abundance(organism, min, max) {
    $.ajax({
        url: "/get_protein_abundance",
        data: {
            organism: organism,
            min: min,
            max: max,
        }
    }).done(function (data) {
        response = JSON.parse(data);
        $("#cu_amount").val(response.portion + " (" + response.log_min + ":" + response.log_max + ")");
        $("#abundance_histogram").html('<img src="' + response.img_path + '" width="240" />');
    });
}

// Activate JS elements
function activate_initial_elements() {
    // Chosen.js
    $("#restriction_sites").chosen({placeholder_text_multiple: 'Enter restriction sites or enzymes to filter for ...'});

    // Slider.js (jquery-ui.js)
    $("#slider-vertical_gc").slider({
        value: 80,
        min: 0,
        max: 100,
        step: 10,
        slide: function (event, ui) {
            $("#gc_amount").val(ui.value);
        }
    });
    $("#gc_amount").val($("#slider-vertical_gc").slider("value"));

    // Slider.js (jquery-ui.js)
    $("#slider-vertical_rnafold").slider({
        value: 80,
        min: 0,
        max: 100,
        step: 10,
        slide: function (event, ui) {
            $("#rnafold_amount").val(ui.value);
        }
    });
    $("#rnafold_amount").val($("#slider-vertical_rnafold").slider("value"));
}

function activate_codon_usage_elements(organism) {
    // Convert Multiple Select in Chosen.js Box
    $("#cu_most_abundant_proteins").chosen({disable_search: true, width: '66px'});
    $("#cu_gb_keywords").chosen({placeholder_text_multiple: 'Enter enzym types to filter for ...', width: '240px'});

    // Slider.js (jquery-ui.js)
    $("#slider-abundance-range").slider({
        orientation: "horizontal",
        range: true,
        values: [0, 100],
        //slide: function (event, ui) {
        //    request_protein_abundance(organism, ui.values[0], ui.values[1]);
        //},
        change: function (event, ui) {
            request_protein_abundance(organism, ui.values[0], ui.values[1]);
        },
    });
    request_protein_abundance(organism, 0, 100);
}

// Activate JS elements
function activate_restriction_elements() {
    // Convert Multiple Select in Chosen.js Box
    $("#restriction_sites").chosen({placeholder_text_multiple: 'Enter restriction sites or enzymes to filter for ...'});
}

// Activate JS elements
function activate_results_gallery() {
    // FancyBox
    $(".fancybox-thumb").fancybox({
        prevEffect: 'none',
        nextEffect: 'none',
        helpers: {
            title: {
                type: 'inside'
            },
            thumbs: {
                width: 50,
                height: 50
            }
        }, // helpers
        afterLoad : function() {
            this.title = (this.title ? '' + this.title + '<br />' : '') + 'Image ' + (this.index + 1) + ' of ' + this.group.length;
        } // afterLoad
    });
    // FancyBox
    $(".fancybox-thumb-single").fancybox({
        prevEffect: 'none',
        nextEffect: 'none',
        helpers: {
            title: {
                type: 'outside'
            },
            thumbs: {
                width: 50,
                height: 50
            }
        }
    });
}


// Execute ONLY on prediction/start page
if (window.location.pathname.substring(1) == '')    {

    // Run on document ready
    $(document).ready(function () {

        // Set Ajax initial configuration
        $.ajaxSetup({
            cache: false
        });

        // Activate improved tooltips
        $('[data-toggle="tooltip"]').tooltip();

        // Modal definitions
        $(modal_id).on('show.bs.modal', function (event) {
            // If necessary, you could initiate an AJAX request here (and then do the updating in a callback).
            // Update the modal's content. We'll use jQuery here, but you could use a data binding library or other methods instead.
            var modal = $(this)
        });
        // .modal-backdrop classes
        $(".modal-transparent").on('show.bs.modal', function () {
            setTimeout(function () {
                $(".modal-backdrop").addClass("modal-backdrop-transparent");
            }, 0);
        });
        $(".modal-transparent").on('hidden.bs.modal', function () {
            $(".modal-backdrop").addClass("modal-backdrop-transparent");
        });

        $(".modal-fullscreen").on('show.bs.modal', function () {
            setTimeout(function () {
                $(".modal-backdrop").addClass("modal-backdrop-fullscreen");
            }, 0);
        });
        $(".modal-fullscreen").on('hidden.bs.modal', function () {
            $(".modal-backdrop").addClass("modal-backdrop-fullscreen");
        });

        // Make textareas auto-resizeable
        function h(e) {
            $(e).css({'height': 'auto', 'overflow-y': 'auto'}).height(e.scrollHeight);
        }

        $('textarea').each(function () {
            h(this);
        }).on('input', function () {
            h(this);
        }).on('focus', function () {
            h(this);
        });

        // Position caret or mark section
        $.fn.selectRange = function (start, end) {
            if (typeof end === 'undefined') {
                end = start;
            }
            return this.each(function () {
                if ('selectionStart' in this) {
                    this.selectionStart = start;
                    this.selectionEnd = end;
                } else if (this.setSelectionRange) {
                    this.setSelectionRange(start, end);
                } else if (this.createTextRange) {
                    var range = this.createTextRange();
                    range.collapse(true);
                    range.moveEnd('character', end);
                    range.moveStart('character', start);
                    range.select();
                }
            });
        };

        // Convert Multiple Select in Chosen.js Box
        activate_initial_elements();

        // Set standard values
        select_host_organism('s.cerevisiae', 'S. cerevisiae');
        set_seq_number(5);

        //  Define functions, which need existent elements
        function get_translation_data() {
            $.ajax({
                url: "/get_dna_translation",
                data: {
                    search_protein: btoa($('textarea[id="search_protein"]').val()),
                    host_organism: $('input[id="host_organism"]').val(),
                },
                beforeSend: function () {
                    $("#loading_modal").modal('show');
                }
            }).done(function (data) {
                response = $.parseJSON(data);
                // Show translated DNA sequence to app user, if one was uploaded
                if (('trans_seq' in response) && response.trans_seq != '') {
                    $('#translated_protein').html(response.trans_seq);
                    $('#translated_protein_div').show();
                } else {
                    $('#translated_protein_div').hide();
                }
            });
        }

        //  Define functions, which need existent elements
        function get_prediction_data() {
            console.log($('input[name="rnafold_switch"]:checked').val());
            $.ajax({
                url: "/get_prediction_data",
                data: {
                    search_protein: btoa($('textarea[id="search_protein"]').val()),
                    host_organism: $('input[id="host_organism"]').val(),
                    seq_number: $('input[id="seq_number"]').val(),
                    restriction_sites: JSON.stringify($('select[id="restriction_sites"]').val()),
                    cu_most_abundant_proteins: JSON.stringify($('select[id="cu_most_abundant_proteins"]').val()),
                    cu_gb_keywords: JSON.stringify($('select[id="cu_gb_keywords"]').val()),
                    cu_abundance_range: JSON.stringify($('#slider-abundance-range').slider('values')),
                    rnafold_switch: $('input[name="rnafold_switch"]:checked').val(),
                },
                beforeSend: function () {
                    $("#loading_modal").modal('show');
                }
            }).done(function (data) {
                // Deactivate loading animation
                $("#loading_modal").modal('hide');

                // Control returned data
                //console.log(data);
                try {
                    response = $.parseJSON(data);
                    console.log(response);
                    // Show general error messages
                    if (('error' in response)) {
                        $(modal_id).find('p').hide();
                        switch (response.error) {
                            case 'empty_sequence':
                                $(modal_id).find('#empty_sequence').toggle();
                                break;
                            case 'extended_iupac':
                                $(modal_id).find('#extended_iupac').toggle();
                                break;
                            default:
                                console.log(response.error);
                                $(modal_id).find('#general').html(response.error);
                                $(modal_id).find('#general').toggle();

                        }
                        $(modal_id).modal();
                        return 0;
                    }
                } catch (e) {
                    // Show data ...
                    $("#result .panel-body").html(data);
                    activate_results_gallery();

                    // Make data visible and scroll to it
                    $("#result").toggle(true);
                    $('html, body').animate({ scrollTop:$('#result').offset().top }, 'slow');
                }
            });
        }

        //  Define functions, which need existent elements
        function get_js_prediction_data() {
            test = $('select[id="cu_most_abundant_proteins"]').val();
            console.log(test);

            $.ajax({
                url: "/get_js_prediction_data",
                data: {
                    search_protein: btoa($('textarea[id="search_protein"]').val()),
                    host_organism: $('input[id="host_organism"]').val(),
                    seq_number: $('input[id="seq_number"]').val(),
                    restriction_sites: JSON.stringify($('select[id="restriction_sites"]').val()),
                    cu_most_abundant_proteins: JSON.stringify($('select[id="cu_most_abundant_proteins"]').val()),
                    cu_gb_keywords: JSON.stringify($('select[id="cu_gb_keywords"]').val()),
                    cu_amount: $('input[id="cu_amount"]').val(),
                },
                beforeSend: function () {
                    $("#loading_modal").modal('show');
                }
            }).done(function (data) {
                $("#loading_modal").modal('hide');
                if (console && console.log) {
                    // Control returned data
                    //console.log(data);
                    response = $.parseJSON(data);
                    console.log(response);

                    // Show general error messages
                    if (('error' in response)) {
                        $(modal_id).find('p').hide();
                        switch (response.error) {
                            case 'empty_sequence':
                                $(modal_id).find('#empty_sequence').toggle();
                                break;
                            case 'extended_iupac':
                                $(modal_id).find('#extended_iupac').toggle();
                                break;
                            default:
                                console.log(response.error);
                                $(modal_id).find('#general').html(response.error);
                                $(modal_id).find('#general').toggle();

                        }
                        $(modal_id).modal();
                        return 0;
                    }

                    // Prepare embedding results in html page
                    prepend_code = '<div name="heprex_result" id="heprex_result" class="code">';
                    append_code = '</div>';
                    prepend_header = '<h4><span class="label label-default">';
                    append_header = '</span></h4>';
                    prepend_figures = '<div name="heprex_result_fig" id="heprex_result_fig" class="code-figures">';
                    append_figures = '</div>';
                    prepend_img_tag = '<img src="';
                    append_img_tag = '" class="fancybox-gallery" alt="Codon usage" width="240px" style="float: left;" />';
                    prepend_fancybox = '<a class="fancybox-thumb" rel="fancybox-thumb" href="';
                    middle_fancybox = '" title="">';
                    append_fancybox = '</a>';
                    content = '';
                    content_img = '';

                    // Make navigation tabs/pills
                    navigation_header = '<ul class="nav nav-pills">\
                        <li role="presentation" class="active"><a onclick="$(this).parents(\'.panel\').children(\'.panel-body\').hide(); $(this).parents(\'.panel\').children(\'#conf_suggestions\').show(); $(this).parents(\'.nav\').children(\'li\').removeClass(\'active\'); $(this).parent().toggleClass(\'active\')">Suggested Sequences</a></li>\
                        <li role="presentation"><a href="#">Suggestions (Colored by Aminoacids)</a></li>\
                        <li role="presentation"><a href="#">Suggestions (Colored by Differences)</a></li>\
                        </ul>';

                    // Step through sequences
                    // TODO: Numbers sprintf
                    for (var key in response.sequences) {
                        content += key + ' ' + response.sequences[key].seq + '\n';
                    }
                    content = prepend_code + content + append_code;
                    // Step through figures
                    for (var key in response.fig_paths) {
                        content_img += prepend_fancybox + response.fig_paths[key] + middle_fancybox + prepend_img_tag + response.fig_paths[key] + append_img_tag + append_fancybox;
                        console.log(response.fig_paths[key]);
                    }

                    // Testing ...
                    cu_test = atob(response.colormap_cu_test);
                    embedded_figs = make_embedded_svg('', cu_test);

                    content_img = prepend_header + 'Codon usages' + append_header + prepend_figures + content_img + embedded_figs + append_figures;

                    // Decode Base64-encoding of sequence alignment view
                    if (response.html_alignment != '') {
                        content_msa = prepend_header + 'Aligned suggestions' + append_header;
                        content_msa += prepend_code + atob(response.html_alignment) + append_code;
                    } else {
                        content_msa = '';
                    }
                    if (response.html_highlight_alignment != '') {
                        content_highlight_msa = prepend_header + 'Matched restriction sites' + append_header;
                        content_highlight_msa += prepend_code + atob(response.html_highlight_alignment) + append_code;
                    } else {
                        content_highlight_msa = '';
                    }

                    // Show translated DNA sequence to app user, if one was uploaded
                    if (('trans_seq' in response)) {
                        $('#translated_protein').html(response.trans_seq);
                        $('#translated_protein_div').show();
                    } else {
                        $('#translated_protein_div').hide();
                    }

                    // Stitch content parts together
                    $("#result .panel-body").html(navigation_header + content + content_msa + content_highlight_msa + content_img);

                    activate_results_gallery()
                    $("#result").toggle(true);
                }
            });
        }

        function make_embedded_svg_image(new_element_id, html_content) {
            // Create new DOM element
            var fancyBoxElem = $('<a class="fancybox-thumb" rel="fancybox-thumb" href="" title=""></a>');
            var canvasElem = $('<canvas id="canvas_test" class="fancybox-gallery" alt="Codon usage" style="border:1px solid #CCC; border-radius: 4px;" width="800" height="800"></canvas>');
            var canvas = $('canvas', canvasElem).get(0).outerHTML;

            // Use exsting DOM element
            //var canvas = document.getElementById('canvas');
            //console.log(canvas);

            // Create SVG element on Canvas
            var ctx = canvas.getContext("2d");
            var data = "<svg xmlns='http://www.w3.org/2000/svg' width='800' height='800'>" +
                "<foreignObject width='100%' height='100%'>" +
                "<style type='text/css'>.cu_table { padding: 2px; border: 1px solid black; font-family: Consolas, Courier, 'Courier New'; font-size: 14px; } </style>" +
                "<div xmlns='http://www.w3.org/1999/xhtml' style='font-size:40px'>" +
                html_content +
                "</div>" +
                "</foreignObject>" +
                "</svg>";

            // Turn Canvas into BLOB
            var DOMURL = self.URL || self.webkitURL || self;
            var img = new Image();
            var svg = new Blob([data], {type: "image/svg+xml;charset=utf-8"});
            var url = DOMURL.createObjectURL(svg);
            img.onload = function () {
                ctx.drawImage(img, 0, 0);
                DOMURL.revokeObjectURL(url);
            };
            img.src = url;

            // Integrate into FancyBox
            $('a', fancyBoxElem).attr('href', url);  // set the attribute
            fancyBoxElem.append(canvas);
            //$("#result .panel-body").append(fancyBoxElem);
            console.log(fancyBoxElem.get(0).outerHTML);
            return fancyBoxElem.get(0).outerHTML;
        }

        function make_embedded_svg(new_element_id, html_content) {
            // Create new DOM elements
            var fancyBoxElem = $('<a class="fancybox-thumb" rel="fancybox-thumb" href="" title=""></a>');
            var divElem = $('<div id="canvas" class="fancybox-gallery" alt="Codon usage" style="border:1px solid #CCC; border-radius: 4px;"></div>');

            // Append SVG element on Canvas
            var data = "<svg xmlns='http://www.w3.org/2000/svg' width='560px'>" +
                "<foreignObject width='100%' height='100%'>" +
                "<style type='text/css'>.cu_table { padding: 2px; border: 1px solid #CCC; font-family: Consolas, Courier, 'Courier New'; font-size: 14px; }</style>" +
                "<div xmlns='http://www.w3.org/1999/xhtml'>" +
                html_content +
                "</div>" +
                "</foreignObject>" +
                "</svg>";

            // Integrate into FancyBox
            divElem.append(data);
            fancyBoxElem.append(divElem);
            //console.log(fancyBoxElem.get(0).outerHTML);
            // Return text version of jQuery-Object
            return fancyBoxElem.get(0).outerHTML;
        }

        function make_alignment() { /* ... */
        }

        function fill_input_sequence(data) {
            $('input[id="show_protein"]').val(data);
            textArea = $('textarea[id="search_protein"]');
            textArea.html(data).focus();
            textArea.selectRange(0);
            textArea.scrollTop();
        }

        function paste_protseq1() {
            seq = ">gi|56966013|pdb|1SZ2|B Chain B, Crystal Structure Of E. Coli Glucokinase In Complex With Glucose\nXGSSHHHHHHGSTKYALVGDVGGTNARLALCDIASGEISQAKTYSGLDYPSLEAVIRVYLEEHKVEVKDGCIAIACPITGDWVAXTNHTWAFSIAEXKKNLGFSHLEIINDFTAVSXAIPXLKKEHLIQFGGAEPVEGKPIAVYGAGTGLGVAHLVHVDKRWVSLPGEGGHVDFAPNSEEEAIILEILRAEIGHVSAERVLSGPGLVNLYRAIVKADNRLPENLKPKDITERALADSCTDCRRALSLFCVIXGRFGGNLALNLGTFGGVFIAGGIVPRFLEFFKASGFRAAFEDKGRFKEYVHDIPVYLIVHDNPGLLGSGAHLRQTLGHIL";
            fill_input_sequence(seq);
        }

        function paste_nucseq2() {
            seq = ">gi|186431|gb|AH002844.1|SEG_HUMINS0 Human insulin gene\nATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG";
            fill_input_sequence(seq);
        }

        function paste_protseq2() {
            seq = ">gi|386828|gb|AAA59172.1| insulin [Homo sapiens]\nMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN";
            fill_input_sequence(seq);
        }

        function paste_nucseq3() {
            seq = ">gi|116242516|sp|P19367.3|HXK1_HUMAN RecName: Full=Hexokinase-1; AltName: Full=Brain form hexokinase; AltName: Full=Hexokinase type I; Short=HK I\nATGATCGCCGCGCAGCTCCTGGCCTATTACTTCACGGAGCTGAAGGATGACCAGGTCAAAAAGATTGACAAGTATCTCTATGCCATGCGGCTCTCCGATGAAACTCTCATAGATATCATGACTCGCTTCAGGAAGGAGATGAAGAATGGCCTCTCCCGGGATTTTAATCCAACAGCCACAGTCAAGATGTTGCCAACATTCGTAAGGTCCATTCCTGATGGCTCTGAAAAGGGAGATTTCATTGCCCTGGATCTTGGTGGGTCTTCCTTTCGAATTCTGCGGGTGCAAGTGAATCATGAGAAAAACCAGAATGTTCACATGGAGTCCGAGGTTTATGACACCCCAGAGAACATCGTGCACGGCAGTGGAAGCCAGCTTTTTGATCATGTTGCTGAGTGCCTGGGAGATTTCATGGAGAAAAGGAAGATCAAGGACAAGAAGTTACCTGTGGGATTCACGTTTTCTTTTCCTTGCCAACAATCCAAAATAGATGAGGCCATCCTGATCACCTGGACAAAGCGATTTAAAGCGAGCGGAGTGGAAGGAGCAGATGTGGTCAAACTGCTTAACAAAGCCATCAAAAAGCGAGGGGACTATGATGCCAACATCGTAGCTGTGGTGAATGACACAGTGGGCACCATGATGACCTGTGGCTATGACGACCAGCACTGTGAAGTCGGCCTGATCATCGGCACTGGCACCAATGCTTGCTACATGGAGGAACTGAGGCACATTGATCTGGTGGAAGGAGACGAGGGGAGGATGTGTATCAATACAGAATGGGGAGCCTTTGGAGACGATGGATCATTAGAAGACATCCGGACAGAGTTTGACAGGGAGATAGACCGGGGATCCCTCAACCCTGGAAAACAGCTGTTTGAGAAGATGGTCAGTGGCATGTACTTGGGAGAGCTGGTTCGACTGATCCTAGTCAAGATGGCCAAGGAGGGCCTCTTATTTGAAGGGCGGATCACCCCGGAGCTGCTCACCCGAGGGAAGTTTAACACCAGTGATGTGTCAGCCATCGAAAAGAATAAGGAAGGCCTCCACAATGCCAAAGAAATCCTGACCCGCCTGGGAGTGGAGCCGTCCGATGATGACTGTGTCTCAGTCCAGCACGTTTGCACCATTGTCTCATTTCGCTCAGCCAACTTGGTGGCTGCCACACTGGGCGCCATCTTGAACCGCCTGCGTGATAACAAGGGCACACCCAGGCTGCGGACCACGGTTGGTGTCGACGGATCTCTTTACAAGACGCACCCACAGTATTCCCGGCGTTTCCACAAGACTCTAAGGCGCTTGGTGCCAGACTCCGATGTGCGCTTCCTCCTCTCGGAGAGTGGCAGCGGCAAGGGGGCTGCCATGGTGACGGCGGTGGCCTACCGCTTGGCCGAGCAGCACCGGCAGATAGAGGAGACCCTGGCTCATTTCCACCTCACCAAAGACATGCTGCTGGAGGTGAAGAAGAGGATGCGGGCCGAGATGGAGCTGGGGCTGAGGAAGCAGACGCACAACAATGCCGTGGTTAAGATGCTGCCCTCCTTCGTCCGGAGAACTCCCGACGGGACCGAGAATGGTGACTTCTTGGCCCTGGATCTTGGAGGAACCAATTTCCGTGTGCTGCTGGTGAAAATCCGTAGTGGGAAAAAGAGAACGGTGGAAATGCACAACAAGATCTACGCCATTCCTATTGAAATCATGCAGGGCACTGGGGAAGAGCTGTTTGATCACATTGTCTCCTGCATCTCTGACTTCTTGGACTACATGGGGATCAAAGGCCCCAGGATGCCTCTGGGCTTCACGTTCTCATTTCCCTGCCAGCAGACGAGTCTGGACGCGGGAATCTTGATCACGTGGACAAAGGGTTTTAAGGCAACAGACTGCGTGGGCCACGATGTAGTCACCTTACTAAGGGATGCGATAAAAAGGAGAGAGGAATTTGACCTGGACGTGGTGGCTGTGGTCAACGACACAGTGGGCACCATGATGACCTGTGCTTATGAGGAGCCCACCTGTGAGGTTGGACTCATTGTTGGGACCGGCAGCAATGCCTGCTACATGGAGGAGATGAAGAACGTGGAGATGGTGGAGGGGGACCAGGGGCAGATGTGCATCAACATGGAGTGGGGGGCCTTTGGGGACAACGGGTGTCTGGATGATATCAGGACACACTACGACAGACTGGTGGACGAATATTCCCTAAATGCTGGGAAACAAAGGTATGAGAAGATGATCAGTGGTATGTACCTGGGTGAAATCGTCCGCAACATCTTAATCGACTTCACCAAGAAGGGATTCCTCTTCCGAGGGCAGATCTCTGAGACGCTGAAGACCCGGGGCATCTTTGAGACCAAGTTTCTCTCTCAGATCGAGAGTGACCGATTAGCACTGCTCCAGGTCCGGGCTATCCTCCAGCAGCTAGGTCTGAATAGCACCTGCGATGACAGTATCCTCGTCAAGACAGTGTGCGGGGTGGTGTCCAGGAGGGCCGCACAGCTGTGTGGCGCAGGCATGGCTGCGGTTGTGGATAAGATCCGCGAGAACAGAGGACTGGACCGTCTGAATGTGACTGTGGGAGTGGACGGGACACTCTACAAGCTTCATCCACACTTCTCCAGAATCATGCACCAGACGGTGAAGGAACTGTCACCAAAATGTAACGTGTCCTTCCTCCTGTCTGAGGATGGCAGCGGCAAGGGGGCCGCCCTCATCACGGCCGTGGGCGTGCGGTTACGCACAGAGGCAAGCAGCTAA";
            fill_input_sequence(seq);
        }

        function paste_protseq3() {
            seq = ">gi|116242516|sp|P19367.3|HXK1_HUMAN RecName: Full=Hexokinase-1; AltName: Full=Brain form hexokinase; AltName: Full=Hexokinase type I; Short=HK I\nMIAAQLLAYYFTELKDDQVKKIDKYLYAMRLSDETLIDIMTRFRKEMKNGLSRDFNPTATVKMLPTFVRSIPDGSEKGDFIALDLGGSSFRILRVQVNHEKNQNVHMESEVYDTPENIVHGSGSQLFDHVAECLGDFMEKRKIKDKKLPVGFTFSFPCQQSKIDEAILITWTKRFKASGVEGADVVKLLNKAIKKRGDYDANIVAVVNDTVGTMMTCGYDDQHCEVGLIIGTGTNACYMEELRHIDLVEGDEGRMCINTEWGAFGDDGSLEDIRTEFDREIDRGSLNPGKQLFEKMVSGMYLGELVRLILVKMAKEGLLFEGRITPELLTRGKFNTSDVSAIEKNKEGLHNAKEILTRLGVEPSDDDCVSVQHVCTIVSFRSANLVAATLGAILNRLRDNKGTPRLRTTVGVDGSLYKTHPQYSRRFHKTLRRLVPDSDVRFLLSESGSGKGAAMVTAVAYRLAEQHRQIEETLAHFHLTKDMLLEVKKRMRAEMELGLRKQTHNNAVVKMLPSFVRRTPDGTENGDFLALDLGGTNFRVLLVKIRSGKKRTVEMHNKIYAIPIEIMQGTGEELFDHIVSCISDFLDYMGIKGPRMPLGFTFSFPCQQTSLDAGILITWTKGFKATDCVGHDVVTLLRDAIKRREEFDLDVVAVVNDTVGTMMTCAYEEPTCEVGLIVGTGSNACYMEEMKNVEMVEGDQGQMCINMEWGAFGDNGCLDDIRTHYDRLVDEYSLNAGKQRYEKMISGMYLGEIVRNILIDFTKKGFLFRGQISETLKTRGIFETKFLSQIESDRLALLQVRAILQQLGLNSTCDDSILVKTVCGVVSRRAAQLCGAGMAAVVDKIRENRGLDRLNVTVGVDGTLYKLHPHFSRIMHQTVKELSPKCNVSFLLSEDGSGKGAALITAVGVRLRTEASS";
            fill_input_sequence(seq);
        }

        function paste_nucseq4() {
            seq = ">lcl|NG_011486.1_cds_NP_005451.2_1 [gene=SNCAIP] [protein=synphilin-1 isoform 1] [protein_id=NP_005451.2]\nATGGAAGCCCCTGAATACCTTGATTTGGATGAAATTGACTTTAGTGATGACATATCTTATTCAGTCACATCACTCAAGACGATCCCAGAACTGTGCCGAAGATGTGATACGCAAAACGAAGACAGATCAGTTTCTAGCTCTAGCTGGAATTGTGGCATCTCAACTCTTATTACAAACACGCAAAAGCCCACAGGAATCGCTGATGTGTACAGTAAGTTCCGCCCAGTGAAGCGGGTTTCGCCACTGAAACATCAGCCAGAGACTCTGGAGAACAATGAAAGTGATGACCAAAAGAACCAGAAAGTGGTTGAGTACCAGAAAGGGGGTGAGTCTGACCTGGGCCCCCAGCCTCAGGAGCTTGGCCCTGGAGATGGAGTGGGCGGCCCACCAGGTAAGAGCTCTGAGCCCAGCACATCGCTGGGTGAACTGGAGCACTACGACCTCGACATGGATGAGATTCTGGATGTGCCTTATATTAAATCCAGTCAGCAGCTTGCCTCTTTTACCAAGGTGACTTCAGAAAAAAGAATTTTGGGCTTATGCACAACCATCAATGGCCTTTCTGGCAAAGCCTGCTCTACAGGAAGTTCTGAGAGCTCATCATCCAACATGGCACCATTTTGTGTTCTTTCTCCCGTGAAAAGCCCTCACTTGAGAAAAGCATCAGCTGTCATCCACGACCAGCACAAGCTGTCCACTGAAGAAACCGAGATCTCACCTCCTCTGGTTAAATGTGGCTCTGCATATGAGCCTGAAAACCAGAGTAAAGACTTCCTAAACAAGACATTTAGTGATCCTCATGGTCGAAAAGTTGAGAAGACAACACCAGACTGCCAGCTCAGGGCCTTCCACCTACAATCCTCAGCAGCAGAATCCAAACCAGAAGAGCAGGTCAGTGGCCTAAACCGGACCAGCTCCCAAGGCCCAGAAGAAAGGAGTGAGTATCTGAAAAAAGTGAAAAGCATCTTGAACATTGTTAAAGAAGGACAGATCTCTCTCCTGCCACACCTAGCTGCAGACAATCTAGACAAAATTCACGACGAAAATGGAAACAATCTATTACATATTGCGGCGTCACAGGGACACGCAGAGTGTCTACAGCACCTCACTTCTTTGATGGGAGAAGACTGCCTCAATGAGCGCAACACTGAGAAGTTGACTCCAGCAGGCCTGGCCATTAAGAATGGTCAGTTGGAGTGCGTACGCTGGATGGTGAGCGAAACAGAAGCCATTGCAGAACTGAGTTGTTCTAAGGATTTTCCAAGCCTTATTCATTACGCAGGTTGCTATGGCCAGGAAAAGATTCTTCTGTGGCTTCTTCAGTTTATGCAAGAACAGGGCATCTCGTTGGATGAAGTAGACCAGGATGGCAACAGTGCCGTTCACGTAGCCTCACAGCATGGCTACCTTGGATGCATACAGACCTTGGTTGAATATGGAGCAAATGTCACCATGCAGAACCACGCTGGGGAAAAGCCCTCCCAGAGCGCCGAGCGGCAGGGGCACACCCTGTGCTCCAGGTACCTGGTGGTGGTGGAGACCTGCATGTCGCTGGCCTCTCAAGTGGTGAAGTTAACCAAGCAGCTAAAGGAACAAACAGTAGAACGTGTCACGCTGCAGAACCAACTCCAACAATTTCTAGAAGCCCAGAAATCAGAGGGCAAGTCACTCCCTTCTTCACCCAGTTCACCATCCTCACCTGCCTCCAGAAAGTCCCAGTGGAAATCTCCAGATGCAGATGATGATTCTGTAGCCAAAAGCAAGCCAGGAGTCCAAGAGGGGATTCAGGTTCTTGGAAGCCTGTCAGCCTCCAGCCGGGCTAGACCCAAAGCAAAAGATGAAGATTCTGATAAAATCTTACGCCAGTTATTGGGAAAGGAAATCTCAGAAAATGTCTGCACCCAGGAAAAACTGTCCTTGGAATTCCAGGATGCTCAGGCTTCCTCTAGAAATTCTAAAAAGATCCCACTGGAGAAGAGGGAACTGAAGTTAGCCAGGCTGAGACAGCTGATGCAGAGGTCACTGAGTGAGTCTGACACAGACTCCAACAACTCTGAGGACCCCAAGACTACCCCAGTGAGGAAGGCTGACCGACCAAGGCCGCAGCCCATTGTAGAAAGCGTAGAGAGTATGGACAGCGCAGAAAGCCTGCACCTGATGATTAAGAAACACACCTTGGCATCAGGGGGACGCAGGTTTCCTTTCAGCATCAAGGCCTCCAAATCCCTGGATGGCCACAGCCCATCTCCCACCTCAGAGAGCAGCGAACCAGACTTAGAATCCCAGTATCCAGGCTCAGGGAGTATTCCTCCAAACCAGCCCTCTGGTGACCCTCAGCAGCCCAGCCCTGACAGTACTGCTGCCCAGAAAGTTGCCACAAGTCCCAAGAGTGCCCTCAAGTCTCCATCTTCCAAGCGTAGGACATCTCAGAACTTAAAACTGAGAGTTACCTTTGAGGAGCCTGTGGTGCAGATGGAGCAGCCTAGCCTTGAACTGAATGGAGAAAAAGACAAAGATAAGGGCAGGACTCTCCAGCGGACCTCCACAAGTAACGAATCGGGGGATCAACTGAAAAGGCCTTTTGGAGCCTTTCGATCTATCATGGAGACACTAAGTGGCAACCAAAACAATAATAATAACTACCAGGCAGCCAACCAGCTGAAAACCTCTACATTGCCCTTGACCTCACTTGGGAGGAAGACAGATGCCAAGGGAAACCCTGCCAGCTCCGCTAGCAAAGGAAAGAATAAGGCAGCATAA";
            fill_input_sequence(seq);
        }

        function paste_protseq4() {
            seq = ">lcl|NG_011486.1_prot_NP_005451.2_1 [gene=SNCAIP] [protein=synphilin-1 isoform 1] [protein_id=NP_005451.2]\nMEAPEYLDLDEIDFSDDISYSVTSLKTIPELCRRCDTQNEDRSVSSSSWNCGISTLITNTQKPTGIADVYSKFRPVKRVSPLKHQPETLENNESDDQKNQKVVEYQKGGESDLGPQPQELGPGDGVGGPPGKSSEPSTSLGELEHYDLDMDEILDVPYIKSSQQLASFTKVTSEKRILGLCTTINGLSGKACSTGSSESSSSNMAPFCVLSPVKSPHLRKASAVIHDQHKLSTEETEISPPLVKCGSAYEPENQSKDFLNKTFSDPHGRKVEKTTPDCQLRAFHLQSSAAESKPEEQVSGLNRTSSQGPEERSEYLKKVKSILNIVKEGQISLLPHLAADNLDKIHDENGNNLLHIAASQGHAECLQHLTSLMGEDCLNERNTEKLTPAGLAIKNGQLECVRWMVSETEAIAELSCSKDFPSLIHYAGCYGQEKILLWLLQFMQEQGISLDEVDQDGNSAVHVASQHGYLGCIQTLVEYGANVTMQNHAGEKPSQSAERQGHTLCSRYLVVVETCMSLASQVVKLTKQLKEQTVERVTLQNQLQQFLEAQKSEGKSLPSSPSSPSSPASRKSQWKSPDADDDSVAKSKPGVQEGIQVLGSLSASSRARPKAKDEDSDKILRQLLGKEISENVCTQEKLSLEFQDAQASSRNSKKIPLEKRELKLARLRQLMQRSLSESDTDSNNSEDPKTTPVRKADRPRPQPIVESVESMDSAESLHLMIKKHTLASGGRRFPFSIKASKSLDGHSPSPTSESSEPDLESQYPGSGSIPPNQPSGDPQQPSPDSTAAQKVATSPKSALKSPSSKRRTSQNLKLRVTFEEPVVQMEQPSLELNGEKDKDKGRTLQRTSTSNESGDQLKRPFGAFRSIMETLSGNQNNNNNYQAANQLKTSTLPLTSLGRKTDAKGNPASSASKGKNKAA";
            fill_input_sequence(seq);
        }

        function paste_nucseq5() {
            seq = ">AY049786.1 Homo sapiens synuclein alpha (SNCA) mRNA, complete cds\nATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA";
            // seq = ">Alpha-Synuclein (140 amino acids)\nATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA";
            fill_input_sequence(seq);
        }

        function paste_protseq5() {
            seq = ">sp|P42212.1|GFP_AEQVI RecName: Full=Green fluorescent protein\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*";
            // seq = ">GFP (238 amino acids)\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLKCFARYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*";
            fill_input_sequence(seq);
        }

        function paste_protseq6() {
            seq = ">sp|P37840.1|SYUA_HUMAN RecName: Full=Alpha-synuclein; AltName: Full=Non-A beta component of AD amyloid; AltName: Full=Non-A4 component of amyloid precursor; Short=NACP\nMDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA*";
            fill_input_sequence(seq);
        }

        // Define submit button action ('bind click event')
        $('button#run_btn').bind('click', function (event) {
            get_translation_data();
            get_prediction_data();
        });
        $('button#seq1').bind('click', function (event) {
            paste_protseq1();
        });
        $('button#nucseq2').bind('click', function (event) {
            paste_nucseq2();
        });
        $('button#protseq2').bind('click', function (event) {
            paste_protseq2();
        });
        $('button#nucseq3').bind('click', function (event) {
            paste_nucseq3();
        });
        $('button#protseq3').bind('click', function (event) {
            paste_protseq3();
        });
        $('button#nucseq4').bind('click', function (event) {
            paste_nucseq4();
        });
        $('button#protseq4').bind('click', function (event) {
            paste_protseq4();
        });
        $('button#nucseq5').bind('click', function (event) {
            paste_nucseq5();
        });
        $('button#protseq5').bind('click', function (event) {
            paste_protseq5();
        });
    });

}