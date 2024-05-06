Change Log
==========

Version 0.1.0 (20 May, 2022)

Version 0.2.0 (03 August, 2022)
   - Find automatically the mkdb.pm, so scripts can be run from anywhere (providing the path to them)
   - get_subtaxa.pl - new script 
   - download_bold.pl - option to cut up taxa to subtaxa with less than max_record_n specimen records
   - format_db.pl - new option for VTAM format
   - add benchmarking of select_taxa
   
Version 0.2.1 (07 August, 2022)

    - format_db.pl Correct if there is a taxID 0 in taxonomy file
    
Version 0.2.2 (09 January, 2023)

    - select_taxa.pl recognizes merged taxids

Version 0.2.3 (24 January, 2023)

    - format_rdp.pl added to make sequence.tsv and taxonomy.tsv from RDP training files
    
Version 0.2.4 (17 April, 2023)

    - download_bold.pl works if taxon_list has a heading (taxon_name) and also if it has not

Version 0.3.0 (06 May, 2024)

    - format_bold_package.pl replaces download_bold.pl sinc it is quicker and BOLD API do not work any more for large data download
    - reduce_metadata.pl to create a tsv file with BOLD metadata of the BOLD sequences kept in COInr
    - for better traceability of the BOLD sequences, ids are in the following format: BOLD_markercode_processid
