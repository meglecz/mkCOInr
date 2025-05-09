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

Version 0.3.1 (28 Oct, 2024)

    - format_db.pl add sintax option to produce database for SINTAX
    
    
Version 0.4.0 (09 May, 2025)

    - format_bold_package.pl

        - Read input TSV line by line, to reduce memory need
        - Can delete sequneces without BOLD BIN; New argument delete_noBIN [0/1/2] 
        - Sequence IDs are in the following format: BOLD_MARKER_PROCESSID_BIN (BOLD_COI-5P_JRPAA9741-15_BIN:ADQ9721 or BOLD_COI-5P_GBBAC2495-15_BIN:NA)

    - add_taxids.pl

        - avoid using 0 as a taxid
        - if more than one taxid for name

            - takes the one with highest proportion of taxa matching the upwards lineage (as before)
            - then smallest difference in taxlevel
            - then prefer a valid scientific name if choice between synonyms and scientific names
            - taxon name in BOLD mathcing a scientific name of the taxid (for homotypic synomyms it is possible to have different scientific names)

