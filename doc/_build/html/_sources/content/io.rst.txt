.. _io_formats_io:

Input/Output
============================

This list of arguments and options can be obtained by typing

.. code-block:: bash

	perl name_of_the_script.pl -help


Taxonomy related files or directories
-------------------------------------------------

.. _ncbitax_dir_io:

ncbitax_dir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Directory containing the NCBI taxonomy dmp files (downloaded by :ref:`download_taxonomy.pl<download_taxonomy_reference>`)

.. _taxonomy_io:

taxonomy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns: 

    - tax_id: can be correct NCBI taxID, or arbitrary negative numbers for taxa not in NCBI
    - parent_tax_id: taxiID of the closest parent of tax_id
    - rank: taxonomic rank (e.g. species, genus, subgenus, no_rank)
    - name_txt: Scientifc name of the taxon
    - old_tax_id: taxIDs that have been merged to the tax_id by NCBI; if more than one for a given tax_id, make one line for each old_tax_id
    - :ref:`taxlevel<taxlevel_index_glossary>`
    - synonyms: list of synonyms

.. code-block:: bash

	tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel	synonyms
	1	1	no rank	root		0	
	2	131567	superkingdom	Bacteria		1	Prokaryotae;Prokaryota;Procaryotae
	6	335928	genus	Azorhizobium		7	
	7	6	species	Azorhizobium caulinodans	395	8	Azotirhizobium caulinodans
	-34968	2778801	subfamily	Callithamnioideae		6.5	
	-35035	-35034	species	Fractonotus caelatus		8	
	-35036	-35030	family	Ramazzottidae		6	


.. _rdp_classifier_taxid_file_io:

RDP classifier taxid file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Taxonomy file in RDP Classifier format. Can be downloaded from https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/
Text file with the following columns separated by stars: 

    - tax_id
    - taxon_name
    - parent_tax_id: taxiID of the closest parent of tax_id
    - rank_index: place in the lineage (starting by 0)
    - rank: taxonomic rank (e.g. species, genus)


.. code-block:: bash

	0*Root*-1*0*rootrank
	1*Bacteria*0*1*domain
	2*Actinobacteria*1*2*phylum
	3*Actinobacteria*2*3*class
	4*Acidimicrobiia*2*3*class


.. _taxon_list_io:

taxon_list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with a list of taxa in the first column. Optionally, taxIDs in the second column.

Any other columns can follow, they will not be taken into account.

.. code-block:: bash

	taxon_name	taxID
	Plecoptera	50622
	Lepidoptera
	Diptera


Lineage files
-------------------------------------------------

.. _bold_data_package_io:

BOLD data package TSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    -processid
    -sampleid
    -specimenid
    -museumid
    -fieldid
    -inst
    -bin_uri
    -identification
    -funding_src
    -kingdom
    -phylum
    -class
    -order
    -family
    -subfamily
    -genus
    -species
    -subspecies
    -identified_by
    -voucher_type
    -collectors
    -collection_date
    -collection_date_accuracy
    -life_stage
    -sex
    -reproduction
    -extrainfo
    -notes
    -coord
    -coord_source
    -coord_accuracy
    -elev
    -depth
    -elev_accuracy
    -depth_accuracy
    -country
    -province
    -country_iso
    -region
    -sector
    -site
    -collection_time
    -habitat
    -collection_note
    -associated_taxas
    -associated_specimen
    -species_reference
    -identification_method
    -recordset_code_arr
    -gb_acs
    -marker_code
    -nucraw
    -sequence_run_site
    -processid_minted_date
    -sequence_upload_date
    -identification_rank

The reduced metadata file do not contain the sequence, and statrs with the sequence ID used in COInd (BOLD_marker_processid format)

.. code-block:: bash

	processid	sampleid	specimenid	museumid	fieldid	inst	bin_uri	identification	funding_src	kingdom	phylum	class	order	family	subfamily	genus	species	subspecies	identified_by	voucher_type	collectors	collection_date	collection_date_accuracy	life_stage	sex	reproduction	extrainfo	notes	coord	coord_source	coord_accuracy	elev	depth	elev_accuracy	depth_accuracy	country	province	country_iso	region	sector	site	collection_time	habitat	collection_noteassociated_taxa	associated_specimen	species_reference	identification_method	recordset_code_arr	gb_acs	marker_code	nucraw	sequence_run_site	processid_minted_date	sequence_upload_date	identification_rank
	AAASF001-17	CBGSFMX-0101	7804897	None	CBGSFMX-0101	Universidad Autonoma de Nuevo Leon	BOLD:ADP3520	Lutzomyia cruciata	None	AnimaliaArthropoda	Insecta	Diptera	Psychodidae	Phlebotominae	Lutzomyia	Lutzomyia cruciata	None	Jorge J. Rodriguez Rojas	None	Wilbert P	2016-10-28	None	Adult	M	S	None	Slide mounted with Euparal	(19.3786,-88.1892)	None	None	None	None	None	None	Mexico	Quintana Roo	None	Candelaria	None	None	None	None	None	None	None	None	Morphological	['AAASF', 'DS-17IBMWP', 'DS-UNIQUE17']	MK851247	COI-5P	AACATTATATTTTATTTTTGGAGCCTGAGCAGGAATAGTGGGAACATCTTTAAGAATTTTAATTCGAGCAGAATTAGGTCACCCCGGTGCTTTAATTGGTGATGATCAAATTTATAATGTTATTGTTACAGCTCATGCATTTGTAATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAACTGATTAGTTCCTTTAATATTAGGAGCCCCTGATATAGCATTCCCTCGAATAAATAATATAAGATTTTGACTTTTACCCCCCTCTCTTACTCTCCTTCTTACAAGAAGTATAGTTGAAACTGGGGCAGGAACAGGATGAACTGTTTATCCACCTCTTTCAAGAAATATTGCCCATAGAGGAGCTTCTGTTGATTTAGCAATTTTTTCCCTACATTTAGCCGGGATTTCATCTATTCTTGGAGCAGTAAATTTTATTACTACAGTTATTAATATACGATCTGCTGGAATTACATTAGATCGAATACCTTTATTTGTTTGATCTGTAATAATTACTGCGGTACTTCTATTATTATCATTACCTGTTTTAGCAGGTGCAATTACAATACTTCTAACTGATCGTAATCTAAATACTTCTTTTTTTGACCCTGCGGGAGGTGGGGATCCAATTTTATATCAACATTTATTT	Instituto Politecnico Nacional, Centro de Biotecnologia Genomica	30-May-2017	12-Jun-2017	species







.. _lineage_tsv_without_taxid_io:

lineage tsv without taxID 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - phylum
    - class
    - order
    - family
    - subfamily
    - genus
    - species
    - seqIDs

All identical lineages are pooled into a single line, seqIDs are in the last column separated by semicolons

.. code-block:: bash

	phylum	class	order	family	subfamily	genus	species	seqIDs
	Acanthocephala							12418139
	Acanthocephala	Archiacanthocephala	Gigantorhynchida	Gigantorhynchidae		Mediorhynchus		5445424;3143887
	Acanthocephala	Archiacanthocephala	Gigantorhynchida	Gigantorhynchidae		Mediorhynchus	Mediorhynchus gallinarum	15188348;15188349;5445423

.. _lineage_tsv_with_taxid_select_taxa_io:

lineage tsv with taxID (output of :ref:`select_taxa.pl<select_taxa_reference>`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - taxon
    - taxID
    - homonymy
    - number of sequences
    - superkingdom
    - kingdom
    - phylum
    - class
    - order
    - family
    - subfamily
    - genus
    - species

.. code-block:: bash

	taxon	taxID	homonymy	number of sequences	superkingdom	kingdom	phylum	class	order	family	genus	species
	Abylidae	316207	0	33	Eukaryota	Metazoa	Cnidaria	Hydrozoa	Siphonophorae	Abylidae		


.. _lineage_tsv_with_taxid_add_taxids_io:

lineage tsv with taxID (output of :ref:`add_taxids.pl<add_taxids_reference>`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - lowest_taxname
    - lowest_rank
    - lowest_TaxID
    - phylum
    - class
    - order
    - family
    - subfamily
    - genus
    - species
    - seqIDs

.. code-block:: bash

	lowest_taxname	lowest_rank	lowest_TaxID	phylum	class	order	family	subfamily	genus	species	seqIDs
	Acanthocephala	phylum	10232	Acanthocephala							12418139
	Mediorhynchus	genus	60535	Acanthocephala	Archiacanthocephala	Gigantorhynchida	Gigantorhynchidae		Mediorhynchus		3143887;5445424


.. _custom_lineages_tsv_io:

custom lineages tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - phylum
    - class
    - order
    - family
    - subfamily
    - genus
    - species
    - homonymy
    - seqIDs

.. code-block:: bash

	phylum	class	order	family	subfamily	genus	species	homonymy	seqIDs
	Cnidaria	Hydrozoa	Leptothecata	Aglaopheniidae		Aglaophenia		0	OEB_MLR10
	Bryozoa	Gymnolaemata	Cheilostomatida	Margarettidae		Margaretta	Margaretta cereoides	1	OEB_EH13;OEB_EH17;OEB_EH19
	Streptophyta	Magnoliopsida	Gentianales	Apocynaceae		Margaretta	Margaretta cereoides	1	OEB_EH13;OEB_EH17;OEB_EH19


.. _ambiguous_lineages_io:

ambiguous lineages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - match_lineage_proportion
    - ncbi_taxname
    - ncbi_taxlevel
    - ncbi_TaxID
    - phylum
    - class	order
    - family
    - subfamily
    - genus
    - species
    - ncbi_superkingdom
    - ncbi_kingdom
    - ncbi_phylum
    - ncbi_class
    - ncbi_order
    - ncbi_family
    - ncbi_genus
    - ncbi_species
    - ncbi_taxname
    - seqIDs

.. code-block:: bash

	match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbi_TaxID	phylum	class	order	family	subfamily	genus	species	ncbi_superkingdom	ncbi_kingdom	ncbi_phylum	ncbi_class	ncbi_order	ncbi_family	ncbi_genus	ncbi_species	ncbi_taxname	seqIDs
	0.4	Bolbophorus	genus	186184	Platyhelminthes	Trematoda	Diplostomida	Diplostomidae	Bolbophorinae	Bolbophorus		Eukaryota	Metazoa	Platyhelminthes	Trematoda	Strigeidida	Bolbophoridae	Bolbophorus		Bolbophorus	12416284;9942141;15268484;12416286;12416287;12416283;12417832;3490428;12417833;5993483;12416282;12416285;12416280;12416281
	0.33	Sylon hippolytes	species	399056	Arthropoda	Hexanauplia		Clistosaccidae		Sylon	Sylon hippolytes	Eukaryota	Metazoa	Arthropoda	Thecostraca		Sylonidae	Sylon	Sylon hippolytes	Sylon hippolytes	2631808;2631807;2631809;2631789;2631806;2631805



Sequence files
-------------------------------------------------

.. _sequence_tsv_without_taxid_io:

sequence tsv without taxID
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - seqID
    - sequence

.. code-block:: bash

	seqID	sequence
	12418139	AGATATTGGTATATTATATATTTTGTTTGCGTTATGAAGAGGC...
	3143887	GTGATATATATAATGTCATCGGTATGAAGTGGTATTATAGGGGTGAT...


.. _sequence_tsv_with_taxid_io:

sequence tsv with taxID
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - seqID
    - taxID
    - sequence

.. code-block:: bash

	seqID	taxID	sequence
	11611742	10236	GGGATAATATATATTTTGCTTGCATTGTGGAGGG...
	10907577	-9466	TAAGATTTTGAATATTACCTCCATCAATTACATT...
	GU179406_1	2921812	GGACTCCTTGGTACTTCTATAAGATTGCTTCTGT...


.. _custom_sequences_tsv_io:

custom sequences tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tsv file with the following columns:

    - seqID
    - taxon name (any taxonomic level)
    - sequence

.. code-block:: bash

	seqID	taxon_name	sequence
	xxx_10236 Porifera	GGGATAATATATATTTTGCTTGCATTGTGGAGGG...
	xxx_10907577	Margaretta	TAAGATTTTGAATATTACCTCCATCAATTACATT...


.. _rdp_classifier_trainset_fasta_io:

RDP classifier trainset fasta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fasta file in RDP Classifier trainseq format. 
Can be downloaded from  https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/

.. code-block:: bash

	>AJ000684	Mycobacterium heidelbergense str. 2554/91 Type	domain__Bacteria; phylum__Actinobacteria; class__Actinobacteria; order__Mycobacteriales; family__Mycobacteriaceae; genus__Mycobacterium
	gaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctctt
	>EF599163	Vibrio atlanticus str. LMG 24300 Type	domain__Bacteria; phylum__Proteobacteria; class__Gammaproteobacteria; order__Vibrionales; family__Vibrionaceae; genus__Vibrio
	gtttgatcctggctcagattgaacgctggcggcaggcctaacacatgcaa



Database formats
-------------------------------------------------

.. _blast_database_files_io:

BLAST database files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Binary files ready to use by BLAST. 

    - blastdb_name.nhr
    - blastdb_name.nin
    - blastdb_name.nog
    - blastdb_name.nsd
    - blastdb_name.nsi
    - blastdb_name.nsq


.. _full_tsv_io:

full tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sequence tsv and the taxonomy files can be formatted by :ref:`format_db.pl<format_db_reference>` to a full tsv file containing the following columns:

    - seqID
    - taxon
    - taxID
    - taxlevel
    - superkingdom
    - superkingdom_taxID
    - kingdom
    - kingdom_taxID
    - phylum
    - phylum_taxID
    - class
    - class_taxID	order
    - order_taxID
    - family
    - family_taxID
    - genus
    - genus_taxID
    - species
    - species_taxID
    - sequence

.. code-block:: bash

	seqID	taxon	taxID	taxlevel	superkingdom	superkingdom_taxID	kingdom	kingdom_taxID	phylum	phylum_taxID	class	class_taxID	order	order_taxID	family	family_taxID	genus	genus_taxID	species	species_taxID	sequence
	5423724	Aspidoscopulia australia	1001026	8	Eukaryota	2759	Metazoa	33208	Porifera	6040	Hexactinellida	60882	Hexactinosida	98040	Farreidae	98041	Aspidoscopulia	999811	Aspidoscopulia_australia1001026	GGATCTCTATTAGAAGACGACCACACCTATAACGTTGTAGTTACAGCTCACGC...


.. _qiime_io:

QIIME 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _qiime_trainseq_fasta_io:

QIIME trainseq fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fasta file with only seqIDs in the definition line

.. code-block:: bash

	>OEB_CA11
	AGTGGTCTCAGTGCTTTAATTCGCATTGAGTTAAGTCAGCCAGGTGGTTTAATGGGCAATG...
	>OEB_EH10
	AGTGGGTAGAGGGTTAAGAGCTTTGATCCGGGTCGAACTAAGTCAACCTGGAGGTTTACTA...



.. _qiime_taxon_file_io:

QIIME taxon file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

file with the following columns:

    - seqID
    - lineage

The taxonomic levels of the lineage are separated by ;

Negative taxIDs are allowed.
Empty taxlevels are filled out using the name of higher-level taxa.

.. code-block:: bash

	OEB_CA11	k__Metazoa_33208; p__Bryozoa_10205; c__Gymnolaemata_10206; o__Cheilostomatida_10207; f__Adeonidae_558780; g__Reptadeonella_2576536; s__Reptadeonella_violacea_-35055
	OEB_EH46	k__Metazoa_33208; p__Bryozoa_10205; c__; o__; f__; g__; s__


.. _rdp_io:

RDP 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _rdp_trainseq_fasta_io:

RDP trainseq fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fasta file with the definition as follows

.. code-block:: bash

	>seqID cellularOrganisms;superkingdom_taxID;kingdom_taxID;phylum_taxID;class_taxID;order_taxID;family_taxID;genus_taxID;species_taxID


Negative taxIDs are allowed
Empty taxlevels are filled out using the name of higher-level taxa (e.g. Polychaeta_6341_order).

.. code-block:: bash

	>MG655623_1 cellularOrganisms;Eukaryota_2759;Metazoa_33208;Ctenophora_10197;Nuda_1919246;Beroida_37538;Beroidae_37539;Beroe_10199;Beroe_forskalii_140453
	ATTTTAGATAAATGATTAGGTTCTGTTTATCATTACAATATTGCTTCTTTATATTTTTTTTTTTCTATTTCTTTAGGGTTTTGTGCCTTTTTTTATTCTTTTATTATAAGATTGTCTTTAGTTTGGCCTTTTGCATTTCTATCTTCAGGTTCTATCTATTTGCATTACGTTACTT
	>7437763 cellularOrganisms;Eukaryota_2759;Metazoa_33208;Annelida_6340;Polychaeta_6341;Polychaeta_6341_order;Orbiniidae_46603;Orbinia_195262;Orbinia_johnsoni_-91
	CGAACAGAACTAGGCCAACCCGGCTCTCTTCTTGGAAGAGACCAACTATACAATACAATTGTTACCGCTCACGCAGTATTAATAATTTTCTTTCTTGTAATGCCCGTCCTAATTGGAGGATTTGGCAACTGACTTGTCCCTTTAAT


.. _rdp_taxon_file_io:

RDP taxon file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

file with the following columns separated by stars:

    - taxID
    - taxon_name_taxID
    - parent taxID
    - taxonomic rank index ('root',1, 'superkingdom',2, 'kingdom',3, 'phylum',4, 'class',5, 'order',6, 'family',7, 'genus',8, 'species')
    - taxonomic rank 

Negative taxIDs are allowed.
Empty taxlevels are filled out using the name of higher-level taxa.

.. code-block:: bash

	-1*Acanthogyrus_cheni_-1*2493664*8*species
	-10*Amynthas_sexpectatus_-10*195544*8*species
	-100*Meiodrilus_adhaerens_-100*2723626*8*species
	-1000*Runcinia_erythrina_-1000*486328*8*species
	-10000*Psaltoda_claripennis_-10000*1225615*8*species
	-10001*Psaltoda_flavescens_-10001*1225615*8*species
	...
	-35075*Polychaeta_6341_order*6341*5*order


.. _sintax_io:

SINTAX 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sintax_fasta_io:

SINTAX fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fasta file with the definition line as follows

.. code-block:: bash

	>seqID;tax=k:kingdom_taxID,p:phylum_taxID,c:class_taxID,o:order_taxID,f:family_taxID,g:genus_taxID,s:species_taxID


Negative taxIDs are allowed
Empty taxlevels are filled out using the name of higher-level taxa (e.g. Monostilifera_6227_family).

.. code-block:: bash

	>KF935544_1;tax=k:Metazoa_33208,p:Nemertea_6217,c:Enopla_6225,o:Monostilifera_6227,f:Monostilifera_6227_family,g:Vieitezia_1068817,s:Vieitezia_luzmurubeae_1068818
	ATTTTAGATAAATGATTAGGTTCTGTTTATCATTACAATATTGCTTCTTTATATTTTTTTTTTTCTATTTCTTTAGGGTTTTGTGCCTTTTTTTATTCTTTTATTATAAGATTGTCTTTAGTTTG
	>PP587771_1;tax=k:Metazoa_33208,p:Chordata_7711,c:Mammalia_40674,o:Rodentia_9989,f:Sciuridae_55153,g:Sciurus_10001
	CCTCCTCTAGCAGGAAATCTAGCCCATGCAGGAGCCTCAATAGATCTAACTATTTTCTCACTCCACCTGGCAGGTGTTTCCTCCATCTTAGGGGCAATTAATTTTATTACTACTATTATCAATAT
	>BOLD_COI-5P_YBIVV3784-23;tax=k:Metazoa_33208,p:Arthropoda_6656,c:Insecta_50557,o:Diptera_7147,f:Rhagionidae_92609,g:Chrysopilus_124301,s:Chrysopilus_alaskaensis_-9996
	TTTATATTTTATCTTTGGAGCTTGAGCGGGTATAGTAGGAACATCTCTTAGTATATTAATTCGAGCAGAATTAGGCCATCCTGGAGCATTAATTGGTGACGATCAAATTTATAATGTGATTGTAA


.. _vtam_database_files_io:

VTAM database files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BLAST database binary files ready to use by BLAST. 

    - blastdb_name.nhr
    - blastdb_name.nin
    - blastdb_name.nog
    - blastdb_name.nsd
    - blastdb_name.nsi
    - blastdb_name.nsq

Taxonomy file with the following columns:

    - tax_id
    - parent_tax_id
    - rank
    - name_txt
    - old_tax_id (old_tax_id merged to tax_id)
    - taxlevel
    
.. code-block:: bash

	tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel
	1	1	no rank	root		0
	2	131567	superkingdom	Bacteria		1
	6	335928	genus	Azorhizobium		7
	7	6	species	Azorhizobium caulinodans	395	8
	9	32199	species	Buchnera aphidicola	28241	8
	10	1706371	genus	Cellvibrio		7
	11	1707	species	Cellulomonas gilvus		8
	13	203488	genus	Dictyoglomus		7
	14	13	species	Dictyoglomus thermophilum		8


Other
-------------------------------------------------

.. _outdir_io:

outdir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Name of the directory to write output files

.. _out_io:

out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

String for naming output files
