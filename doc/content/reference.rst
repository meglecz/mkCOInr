Description
===================================

.. _add_taxids_reference:

add_taxids.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each lineage

    - Find the smallest taxon that matches an already existing :ref:`taxID<taxid_glossary>`  in :ref:`taxonomy.tsv<taxonomy_io>`
    - Assign arbitrary, negative taxIDs to taxa under this taxon

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - lineages (:ref:`lineage tsv without taxID<lineage_tsv_without_taxid_io>`; output of :ref:`format_custom.pl<format_custom_reference>` or :ref:`format_bold.pl<format_bold_reference>`)
    - sequences (:ref:`sequence tsv without taxID<sequence_tsv_without_taxid_io>`)
    - :ref:`outdir<outdir_io>`
    - :ref:`taxonomy.tsv<taxonomy_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NONE

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each lineage finds the lowest taxon that matches an already existing :ref:`taxID<taxid_glossary>` :

    - Finds the smallest taxon that matches an already existing taxID (taking into account all taxIDs in the input taxonomy file)
    - Accepts taxID if at least 0.6 of the taxa in the upward input lineages matches the lineage of the taxID (for species level do not count the genus, since it matches necessarily the species name) OR Both the input taxon and taxID have a phylum rank
    - If the match between the input lineage and the taxID lineage is <=0.25, go to the next taxlevel
    - If the match between the input lineage and the taxID lineage is between 0.25 and 0.6, print lineage to ambiguous file, and it should be checked manually

The search for an existing taxID, takes into account 

    - Existing NCBI or previous arbitrary (negative) taxIDs
    - All synonym taxon names in the taxonomy file
    - If there are more than one taxIDs for the taxon name, take the one with the highest proportion of taxa matching the upwards lineage and then the one that has the same taxonomic rank

If there is no taxID for the lowest taxonomic levels in the lineage, assign arbitrary negative taxIDs to them.
New taxIDs are linked as a child to an existing taxID, and the taxonomy file is updated with them.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`sequences_with_taxIDs.tsv<sequence_tsv_with_taxid_io>`
    - :ref:`lineages_with_taxIDs.tsv<lineage_tsv_with_taxid_add_taxids_io>`
    - :ref:`ambiguous_lineages.tsv<ambiguous_lineages_io>` 
    - :ref:`ambiguous_sequences.tsv<sequence_tsv_without_taxid_io>`
    - :ref:`taxonomy_updated.tsv<taxonomy_io>`


.. _dereplicate_reference:

dereplicate.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Taxonomically aware dereplication.

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - :ref:`tsv sequence tsv with taxID<sequence_tsv_with_taxid_io>`
 - :ref:`outdir<outdir_io>`
 - out (name of the output dereplicated sequence tsv file)

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - vsearch_path (path to vsearch executables if not in the PATH)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare sequences of the same taxID. 
Delete sequences that are substring of another sequence (100% identity on the overlapping region, and one sequence covers the other completely).
If more than 10.000 sequences for the same taxID, first, cluster the sequences with 100% of identity using the cluster_fast algorithm of vsearch, than use the substring search for each cluster.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`tsv sequence tsv with taxID<sequence_tsv_with_taxid_io>`


.. _download_bold_reference:

download_bold.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download BOLD data in tsv format for a list of taxa.

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`taxon_list<taxon_list_io>`
    - :ref:`outdir<outdir_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - try_download (integer; Try to download files *try_download* times if some of the downloaded files are incomplete; Default: 3)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

the scripts downloads all sequences and lineages for all taxa on the taxon_list from BOLD. 

The taxon_list file was constructed manually from taxa on https://www.boldsystems.org/index.php/TaxBrowser_Home. Each taxa on the list has less than 500M specimen records on BOLD. 
The taxon_list constructed on 2022-02-24 is available with on github (data/bold_taxon_list_2022-02-24.txt). It contains all taxa available in BOLD. This file might need to be updated later.

Download is done using BOLD's API. First a small stat file is downloaded to access the number of records available for the taxa, then the tsv file is downloaded with sequences and metainfo.
The stript checks if the number of downloaded records corresponds to the expected one (based on stat file).
If error, it removes the file and retry download try_download times.

If the file exists already, the download is skiped. In this way, if the program stops (for example hitting wall time on a server), it can be simply restarted and the taxa with successful downloads will not be rerun.

NOTE: The download for a long list of taxa can take several days since it is not parallelized. 
You can cut up the input list and run each of them on a separate computer and move the output files to the same folder afterwards.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - json file for each taxon with the total number of records
    - tsv file for each taxon with the following columns:  
    processid	sampleid	recordID	catalognum	fieldnum	institution_storing	collection_code	bin_uri	phylum_taxID	phylum_name	class_taxID	class_name	order_taxID	order_name	family_taxID	family_name	subfamily_taxID	subfamily_name	genus_taxID	genus_name	species_taxID	species_name	subspecies_taxID	subspecies_name	identification_provided_by	identification_method	identification_reference	tax_note	voucher_status	tissue_type	collection_event_id	collectors	collectiondate_start	collectiondate_end	collectiontime	collection_notesite_code	sampling_protocol	lifestage	sex	reproduction	habitat	associated_specimens	associated_taxa	extrainfo	notes	lat	lon	coord_source	coord_accuracy	elev	depth	elev_accuracy	depth_accuracy	country	province_state	region	sector	exactsite	image_ids	image_urls	media_descriptors	captions	copyright_holders	copyright_years	copyright_licenses	copyright_institutions	photographers	sequenceID	markercode	genbank_accession	nucleotides	trace_ids	trace_names	trace_links	run_dates	sequencing_centers	directions	seq_primers	marker_codes

.. _download_taxonomy_reference:

download_taxonomy.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the NCBI taxonomy dmp files (https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/) and prepare :ref:`taxonomy.tsv<taxonomy_io>` file.

Input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`outdir<outdir_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - skip_download (0/1; if 1, skips download, only prepares taxonomy file; Default: 0)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Downloads ncbi taxonomy dump files from https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/ to the ncbi_tax subdirectory and decompress them.
Prepares a :ref:`taxonomy.tsv<taxonomy_io>` file with all taxIDs in it.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - ncbi_tax subdirectory with dmp files
    - :ref:`taxonomy.tsv<taxonomy_io>`

.. _format_bold_reference:

format_bold.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prepare a single sequence file and a lineage file from all tsv files downloaded from BOLD.
Clean and orient sequences.

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - download_dir (name of the folder containing the downloaded BOLD tsv files)
    - :ref:`outdir<outdir_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - marker_list (List of markers to be selected; Default:  'COI-5P COI-3P')
    - check_name (0/1; If 1 keeps only taxa with valid Latin name format: Default: 1)
    - max_n (positive integer; eliminates sequences with max_n or more consecutive Ns; Default:5)
    - min_length (positive integer; minimum length of the cleaned sequence; Default:100)
    - max_length (positive integer; maximum length of the cleaned sequence; Default:2000)
    - check_orientation (0/1; if 1, checks the orientation of the sequences; Default: 0)
    - blast_path (Optional; Path to the BLAST executables if they is not in your PATH)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clean downloaded files and pool information to lineage and sequence files

    - Eliminate partial lines (mostly errors in the database)
    - Select sequences for a given marker list
    - Clean sequences
        - Correct sequence IDs
        - Gaps deleted
        - Non-TCGA changed to N
        - External Ns deleted
        - Sequences with more than max_n consecutive Ns are deleted
        - Keep sequences with length in a min_length and max_length range 
    - Clean lineages
         - If check_name, keep only names matching a correct Latin name format (only letters, spaces and -, correct capitalization)
         - Pool identical lineages into one line with the list of valid sequence IDs in the last field
         - Eliminate lines with environmental and metagenomic samples
    - Orient sequence (optional)
        - Count the TAA, TAG STOP codons in each reading frame
        - Choose the orientation where there is no STOP codon
        - If STOP codon in all frames OR stand + and - among the frames without STOP codon, class it as ambiguous
        - Make a small "reference" db form randomly sampled oriented sequences
        - Blast ambiguous sequences to check orientation
        - Write sequences without hit to the bold_ambiguous_orientation.fas

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`bold_sequences.tsv<sequence_tsv_without_taxid_io>`
    - :ref:`bold_lineages.tsv<lineage_tsv_without_taxid_io>`
    - bold_partial_lines.tsv (lines in the input tsv files that did not have sequences)
    - bold_ambiguous_orientation.fas (sequences that could not be oriented in check_orientation option)

.. _format_custom_reference:

format_custom.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a lineage file for each input taxon name. 
Prepare input for add_taxid.pl

The output lineage file should be checked manually 

    - To see if the suggested lineages are plausible
    - To select the correct lineage if there is more than one (1 in homonymy column) for the same taxon name

Before the next step (add_taxids.pl)

    - Correct/delete/complete lines if lineage in not correct. Try to use taxon names compatible with NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy/ )
    - Delete the homonymy column
    
Input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`custom sequneces<custom_sequences_tsv_io>`
    - :ref:`taxonomy.tsv<taxonomy_io>`
    - :ref:`outdir<outdir_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - max_n (positive integer; eliminate sequences with max_n or more consecutive Ns; Default: 5)
    - min_length (positive integer; minimum length of the cleaned sequence; Default: 100)
    - max_length (positive integer; maximum length of the cleaned sequence; Default: 2000 )
    - check_seqid_format (0/1; if 1 check if seqID resemble to bold and ncbi formats, print out warnings, if yes; Default: 1)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Match names to all taxon names in :ref:`taxonomy.tsv<taxonomy_io>` (including synonyms)

    - Write a lineage to all taxon names where name matches to an existing name in taxonomy.tsv
    - If homonymy, create a lineage for each homonym, and write 1 to the homonymy column
    - If taxon name corresponds to a Latin name format (Genus species) and species name is not known, get the lineage for the genus.

If the check_seqid_format option is activated

    - If some of the sequence IDs are not unique, list the duplicates IDs and exit
    - If sequence ID format is similar to accessions used in BOLD and NCBI/EMBL/DDBJ, list IDs but continue
    - A fairly safe format is xxx_xxx####, where x is a letter, # is a digit

Clean sequences

    - gaps deleted
    - non-TCGA changed to N
    - external Ns deleted
    - sequences with more than max_n consecutive Ns are deleted


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`custom_lineages.tsv<custom_lineages_tsv_io>`
    - :ref:`custom_sequences.tsv<sequence_tsv_without_taxid_io>`

.. _format_db_reference:

format_db.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a database in blast, rdp, qiime or full tsv format from the :ref:`sequence tsv<sequence_tsv_with_taxid_io>` and :ref:`taxonomy.tsv<taxonomy_io>` files

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`sequence tsv<sequence_tsv_with_taxid_io>`
    - :ref:`taxonomy.tsv<taxonomy_io>`
    - outfmt (rdp/blast/qiime/full/vtam; choose the format of the database)
    - :ref:`outdir<outdir_io>`
    - out (string for naming the output files)

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - blast_path (Optional; path to the blast executables if it is not in your PATH)
	
Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BLAST db, VTAM

    - Prepare a fasta file with sequneces and the taxIDs (seqID, taxID)
    - Run the *makeblastdb* commande of blast to make indexed files ready to be used as a blast database
    - for VTAM format, prepare taxonomy file as well as the BLAST database. They can be used directly in VTAM.

RDP, QIIME and FULL

    - Prepare a ranked lineage for each taxID
    - Taxon names are concatenated with taxID to avoid homonymy
    - Missing taxonomic levels are completed by using the name of a higher-level taxon concatenated with the taxonomic ranks
    - Prepare a trainseq fasta and a taxon file for :ref:`rdp<rdp_io>` and :ref:`qiime<qiime_io>`
    - Prepare a single tsv file for :ref:`full<full_tsv_io>`

The  trainseq fasta and the taxon files can be used by the *train* command of rdp_classifier or *feature-classifier* of QIIME to train the dataset before classification

The full tsv format is an easy to parse tsv file with :ref:`ranked lineage<ranked_lineage_glossary>` and :ref:`taxID<taxid_glossary>` for each sequence.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BLAST option

    - :ref:`Indexed files<blast_database_files_io>` ready to be used as a BLAST database 
    
VTAM option

    - :ref:`Indexed files<blast_database_files_io>` ready to be used as a BLAST database 
    - :ref:`taxonomy.tsv<taxonomy_io>` adpted to VTAM

RDP option 

    - :ref:`RDP trainseq fasta<rdp_trainseq_fasta_io>` 
    - :ref:`RDP taxon file<rdp_taxon_file_io>`

QIIME option

    - :ref:`QIIME trainseq fasta<qiime_trainseq_fasta_io>`
    - :ref:`QIIME taxon file<qiime_taxon_file_io>`

FULL option

    - :ref:`tsv file<full_tsv_io>`
	
	
.. _format_ncbi_reference:

format_ncbi.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format the CDS and taxID files to a sequence  :ref:`sequence tsv file with taxIDs<sequence_tsv_with_taxid_io>`.
Select and clean sequences.

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - cds (CDS fasta file; output of nsdpy) 
    - taxids (tsv file with the seqID and taxID columns; output of nsdpy) 
    - :ref:`taxonomy.tsv<taxonomy_io>`
    - :ref:`outdir<outdir_io>`

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - check_name (0/1; If one keep only taxa with valid Latin name fomat: Default: 1)
    - max_n (positive integer; eliminate sequences with max_n or more consecutive Ns; Default: 5)
    - min_length (positive integer; minimum length of the cleaned sequence; Default: 100)
    - max_length (positive integer; maximum length of the cleaned sequence; Default: 2000)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - Select sequences
        - Check if gene names and protein names correspond to COI
        - Eliminate genes if they have introns
        - Can have more than one COI gene in the same mitochondrion
        - Accept only if valid taxID, replace old non-valid taxIDs by up-to-date taxIDs
        - Eliminate sequences from environmental and metagenomic samples
        - If check_name is activated, take the taxID of the smallest taxon with a valid Latin name, otherwise keep the original taxID. 
    - Clean sequences
        - Upper case
        - Replace non-ATCG by N
        - Delete gaps and external Ns
        - Delete sequence if more than max_n consecutive Ns
        - Keep sequences if length is between min_length and max_length
        - Sequences are already in a correct orientation in the input file, since that are coming from CDS files

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - :ref:`ncbi_sequences.tsv<sequence_tsv_with_taxid_io>`
 

.. _pool_and_dereplicate_reference:

pool_and_dereplicate.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pool 2 dereplicated sequence tsv files and do a taxonomically-aware dereplication for taxIDs present in both input files

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`tsv1<sequence_tsv_with_taxid_io>`
    - :ref:`tsv2<sequence_tsv_with_taxid_io>`
    - :ref:`outdir<outdir_io>`
    - out (name of the output dereplicated sequence tsv file)

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - vsearch_path (path to vsearch executables if not in the PATH)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pool sequences from two input tsv files.
The dereplication is done only for taxID shared by the two input files, since they have been dereplicated individually.
The algorithm of dereplication is identical to the one used in :ref:`dereplicate.pl<dereplicate_reference>` 

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :ref:`sequence tsv with taxIDs<sequence_tsv_with_taxid_io>`


.. _select_region_reference:

select_region.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select a target region from input sequences. 
As an input, either a primer pair should be given to identify the target region in some sequences by e-pcr, or a fasta file containing taxonomically diverse sequences limited to the target region.
The sequences are then aligned to the target sequences and trimmed according to the alignment

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`sequence tsv with taxIDs<sequence_tsv_with_taxid_io>`
    - :ref:`outdir<outdir_io>`
    - target_region_fas (A small phylogenetically diverse fasta file with sequences already trimmed to the target region; Optional; Can be produced by e-pcr included in the script.)

Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*E-pcr related parameters*

    - e_pcr  (0/1; if 1, identify the target region of the sequences by e-pcr in the first step)
    - fw (optional; if e_pcr is done, the sequence of the forward primer that amplifies the target region)
    - rv (optional; if e_pcr is done, the sequence of the reverse primer that amplifies the target region)
    - trim_error (real [0-1], the proportion of mismatches allowed between the primer and the sequence during the e_pcr; Default : 0.3)
    - min_overlap (the minimum overlap between primer and the sequence during e-pcr; Defalut 10)
    - min_amplicon_length (The minimum length of the amplicon after primer trimming; Default: 100)
    - max_amplicon_length (The minimum length of the amplicon after primer trimming; Default: 2000)
    - cutadapt_path (Optional; Path to cutadapt if it is not in your PATH)

*usearch_global related parameters* for trimming sequneces according to the alignments

    - tcov_hsp_perc (minimum coverage of the target sequence in *usearch_global* hits; Default: 0.5)
    - perc_identity (minimum percentage of identity between the sequence and the target in *usearch_global* hits; Default: 0.7)
    - vsearch_path (Optional; path to the vsearch if it is not in your PATH)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fasta file is prepared from the input tsv sequence file. 

Sequences are aligned to a small pool of target sequences already limited to the target region (target_region_fas). 
The alignment is done by the *usearch_global* command of vsearch which makes global alignments (unlike BLAST). 

The best hit is used for each sequence to orient and trim them to the target region. 
Only hits with a minimum target coverage (tcov_hsp_perc) and percentage of identity (perc_identity) are used.

The target_region_fas file can be either previously prepared by the users and given as an input, 
or be produced by e-pcr using the e-pcr related parameters.

To reduce runtime, sequences in the target_region_fas are clustered by -cluster_fast algorithm of vsearch 
and the centroids are used as a target file for the *usearch_global*

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - cutadapt_trimmed.fas (if e_pcr; fasta file with sequences recognized and trimmed by E-pcr; equivalent of target_region_fas )
    - target_centroids.fas (fasta file of centroids of the clustering of target_region_fas or cutadapt_trimmed.fas)
    - :ref:`trimmed.tsv<sequence_tsv_with_taxid_io>` (sequence tsv files with taxIDs trimmed to the target region) 
    - :ref:`untrimmed.tsv<sequence_tsv_with_taxid_io>` (sequence tsv files with taxIDs where the target region is not found)

.. _select_taxa_reference:

select_taxa.pl
-------------------------------------------------

Aim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select or eliminate  sequences that belong to taxa on a taxon list.
Select sequences with a minimum taxonomic resolution (e.g. assigned at least to genus).

Input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`sequence tsv with taxIDs<sequence_tsv_with_taxid_io>`
    - :ref:`taxon_list<taxon_list_io>`
    - :ref:`taxonomy<taxonomy_io>`
    - :ref:`outdir<outdir_io>`
    - out (name of the output dereplicated sequence tsv file)
 
Parameters/options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - negative_list (1/0; if 1, keeps all taxa except the ones on the taxon list; Default: 0)
    - min_taxlevel (species/genus/family/order/class/phylum/kingdom/root; Default: root)

Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select sequences that belong to the taxa on a taxon list if *negative_list==0* (default). Select sequences that DO not belong to the taxa on in a taxon list, if *negative_list==1*.

Keep only sequences that are assigned to at least to *min_taxlevel* rank.

If taxIDs are not given in the taxon_list file the script uses all taxIDs that matches the taxon name.

A lineage file is written for each taxon in taxon_list and the corresponding taxIDs.

    - It should be checked manually if lineages are coherent with the target taxa
    - Homonymy column indicates if there are more than one taxID for a taxon
    - If there are incoherent lineages, make a new taxon_list file based on the lineage file including taxon names and taxIDs and rerun the script with the new taxon_list.

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    - :ref:`sequence tsv with taxIDs<sequence_tsv_with_taxid_io>`
    - :ref:`lineage file with taxIDs<lineage_tsv_with_taxid_select_taxa_io>`

