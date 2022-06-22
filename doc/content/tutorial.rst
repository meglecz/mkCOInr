.. _tutorial_tutorial:

Tutorial
============

For the sake of the tutorial, I have created the following file system, and the file names in the examples will contain relative paths.
The scripts directory contains all mkCOInr scripts

.. code-block:: bash

	mkCOInr
	├── custom
	│   └── COI_custom.tsv
	└── scripts
		├── add_taxids.pl
		├── dereplicate.pl
		├── download_bold.pl
		├── download_taxonomy.pl
		├── format_bold.pl
		├── format_custom.pl
		├── format_db.pl
		├── format_ncbi.pl
		├── mkdb.pm
		├── pool_and_dereplicate.pl
		├── select_region.pl
		└── select_taxa.pl


.. _customize_tutorial:

Customize database 
-------------------------------------------------

The creation of the COInr database is explained in the :ref:`Create COInr from BOLD and NCBI section <create_coinr_tutorial>`. You can download this database from `Zenodo <https://doi.org/10.5281/zenodo.6555985>`_ and customize it to your needs.

Download, untar COInr and move the files to the mkCOInr/COInr directory

.. code-block:: bash

	tar-zxvf COInr_2022_05_06.tar.gz
	mkdir -p mkCOInr/COInr
	mv COInr.tsv mkCOInr/COInr
	mv taxonomy.tsv mkCOInr/COInr


This gives the following file structure

.. code-block:: bash

	mkCOInr
	├── COInr
	│   ├── COInr.tsv
	│   └── taxonomy.tsv
	├── custom
	│   └── COI_custom.tsv
	└── scripts
		├── add_taxids.pl
		├── dereplicate.pl
	...



The  :ref:`I/O formats section <io_formats_io>` gives you details about all file formats and examples are provided as well. 

.. _add_custom_sequences_tutorial:

Add custom sequences to a database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _format_custom_tutorial:

Format custom files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`input tsv file <custom_sequences_tsv_io>` (-custom) contains :ref:`seqIDs <seqid_glossary>`, taxon name (can be at any taxonomic level) and sequences.
This script will suggest one or more lineages for each taxon name based on the existing lineages in :ref:`taxonomy.tsv <taxonomy_io>`. It will also consider synonyms.

.. code-block:: bash

	perl format_custom.pl -custom ../custom/COI_custom.tsv -taxonomy ../COInr/taxonomy.tsv -outdir ../custom/format


The output lineage file (custom/format/custom_lineages.tsv) should be checked manually to see if the lineages are coherent.
If homonymy, choose the correct lineage, then delete homonymy column. This revised file will be the input to the add_taxids.pl script.

.. _add_taxids_custom_tutorial:

Add taxIDs to custom sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each lineage in the input file
    - Find an existing taxID at the lowest possible taxonomic level. :ref:`taxIDs <taxid_glossary>` can be either from NCBI, or negative taxID already present in :ref:`taxonomy.tsv <taxonomy_io>`.
    - Add new arbitrary (negative) taxIDs to taxa, that are not yet in taxonomy file 
    - Link each new taxID to existing one as a child and include info to the updated taxonomy file
	Make a :ref:`tsv file with sequences and taxIDs <sequence_tsv_with_taxid_io>`
	Update the :ref:`taxonomy.tsv <taxonomy_io>` file

.. code-block:: bash

	perl add_taxids.pl -lineages ../custom/format/custom_lineages_verified.tsv -sequences ../custom/format/custom_sequences.tsv -taxonomy ../COInr/taxonomy.tsv -outdir ../custom/add_taxids 


.. _dereplicate_custom_tutorial:

Dereplicate custom sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Eliminate sequences that are substring of another sequence of the same :ref:`taxID <taxid_glossary>`. Use :ref:`sequences_with_taxIDs.tsv <sequence_tsv_with_taxid_io>` file (output of the previous script).

.. code-block:: bash

	perl dereplicate.pl -tsv ../custom/add_taxids/sequences_with_taxIDs.tsv -outdir ../custom/dereplicate -out custom_dereplicated_sequences.tsv


The output file is in the same format as the input tsv file.

.. _pool_and_dereplicate_custom_tutorial:

Pool and dereplicate datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use 2 dereplicated :ref:`sequence tsv files <sequence_tsv_with_taxid_io>`

- COInr.tsv  (pool of BOLD and NCBI, downloaded from Zenodo)
- custom_dereplicated_sequences.tsv

Pool the files and dereplicate sequences of the taxIDs that are present in both files

.. code-block:: bash

	perl pool_and_dereplicate.pl -tsv1 ../COInr/COInr.tsv -tsv2 ../custom/dereplicate/custom_dereplicated_sequences.tsv -outdir ../final -out COInr_custom.tsv


The output is the same format as the input tsv file.

Your custom database is composed of two files:

- the dereplicated sequnece file (COInr_custom.tsv)
- the last version of the taxonomy file (custom/add_taxids/taxonomy_updated.tsv)

For simplicity, move the updated taxonomy file to the same folder of the sequence.

.. code-block:: bash

	mv ../custom/add_taxids/taxonomy_updated.tsv ../final/taxonomy_updated.tsv



This database can be further customized, or you can simply format it to be ready for your taxonomic assignment program by the :ref:`format_db.pl <format_db_reference>` script.


.. _select_sequences_custom_tutorial:

Select sequences from existing database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select sequences for a list of taxa with a minimum taxonomic rank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequences can be selected for a list of taxa and/or for a minimum taxonomic level (species/genus/family/order/class/phylum/kingdom/superkingdom/root)

The input file (:ref:`-taxon_list <taxon_list_io>`) contains a list of taxa and eventually their taxIDs. The first line is a heading and will be ignored.

.. code-block:: bash

	perl select_taxa.pl -taxon_list ../final/taxon_list.txt -tsv ../final/COInr_custom.tsv -taxonomy ../final/taxonomy_updated.tsv  -min_taxlevel species  -outdir ../final/selected/ -out COInr_custom_selected.tsv


The main output is a :ref:`sequence tsv file <sequence_tsv_with_taxid_io>` in the same format as the input.
A :ref:`lineage file <lineage_tsv_with_taxID_io>` is also written for all taxa in the taxon_list to check if they are coherent with the target taxon names. 

See more details in the detailed description of the :ref:`select_taxa.pl <select_taxa_reference>` script.

Excluding  sequences of a taxon list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the same script it is also possible to eliminate sequences of taxa instead of selecting them. Set the *negative_list* option to 1 to do that.

.. code-block:: bash

	perl select_taxa.pl -taxon_list ../final/taxon_list.txt -tsv ../final/COInr_custom.tsv -taxonomy ../final/taxonomy_updated.tsv  -min_taxlevel species  -outdir ../final/selected/ -out COInr_custom_reduced.tsv -negative_list 1



.. _select_region_custom_tutorial:

Select region
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequences can be trimmed to a specific region of the COI gene. To define the region, you can either give a fasta file with sequences covering the region of interest, or you can detect them automatically by e-pcr, as it is in this example.

.. code-block:: bash

	perl select_region.pl -tsv ../final/COInr_custom.tsv -outdir ../final/amplicon/ -e_pcr 1 -fw GGNTGAACNGTNTAYCCNCC -rv TAWACTTCDGGRTGNCCRAARAAYCA -trim_error 0.3 -min_amplicon_length 280 -max_amplicon_length 345 -min_overlap 10 -tcov_hsp_perc 0.5 -perc_identity 0.7


See more details in the details description of the :ref:`select_region.pl <select_region_reference>` script.

.. _format_db_custom_tutorial:

Format database 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format the database to one of the following formats
    - blast
    - rdp
    - qiime
    - full

.. code-block:: bash

	perl format_db.pl -tsv ../final/COInr_custom.tsv -taxonomy ../final/taxonomy_updated.tsv -outfmt qiime -outdir ../final/qiime -out COInr_custom_qiime


You should use the rdp_calssifier or qiime's feature-classifier to train the database using the output files of this script if you have used the rdp or qiime options.



The full option, gives a :ref:`tsv file <full_tsv_io>` with seqIDs, ranked lineages, taxIDs for each sequnece, and this is a very easy-to-parse, complete file.

.. code-block:: bash

	perl format_db.pl -tsv ../final/COInr_custom.tsv -taxonomy ../final/taxonomy_updated.tsv -outfmt full -outdir ../final/full -out COInr_custom



For making a BLAST database, the taxonomy file is not necessary and the indexed files in the output folder are ready to use.

.. code-block:: bash

	perl format_db.pl -tsv ../final/COInr_custom.tsv -outfmt blast -outdir ../final/blast -out COInr_custom

See more details in the details description of the :ref:`format_db.pl <format_db_reference>` script.


.. _create_coinr_tutorial:

Create COInr from BOLD and NCBI
-------------------------------------------------
The following steps describe how COInr database (available at `Zenodo <https://doi.org/10.5281/zenodo.6555985>`_ was produced. 

.. _download_ncbi_taxonomy_tutorial:

Download NCBI taxonomy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download NCBI taxonomy dmp file and create :ref:`taxonomy.tsv <taxonomy_io>`.

.. code-block:: bash

	perl download_taxonomy.pl -outdir ../taxonomy


.. _ncbi_sequences_tutorial:

NCBI sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download NCBI sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command will download Coding DNA Sequence (CDS) fasta files of all sequences with COI, CO1, COXI or COX1 in the title lines and complete mitochondrial genomes.
It takes several hours (days) to run this command.

.. code-block:: bash

	nsdpy -r "COI OR COX1 OR CO1 OR COXI OR (complete[Title] AND genome[Title] AND Mitochondrion[Filter])" -T -v --cds


The results are found in the NSDPY_results/yyyy-mm-dd_hh-mm-ss folder.

sequences.fasta contains all CDS sequences. Sequences are correctly oriented but should still be filtered to keep only COI sequences.
TaxIDs.txt contains the sequenceIDs and the TaxIDs.

Move the results of nsdpy to the ../ncbi/nsdpy directory and clean up the directory.

.. code-block:: bash

	mkdir -p ../ncbi
	mv NSDPY_results/yyyy-mm-dd_hh-mm-ss ../ncbi/download
	mv report.tsv ../ncbi/download
	rmdir NSDPY_results


Format NCBI sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    - Select COI sequences and clean them. 
    - Eliminate identical sequences of the same taxID.
    - Clean tax names and taxids.

.. code-block:: bash

	perl format_ncbi.pl -cds ../ncbi/download/sequences.fasta -taxids ../ncbi/download/TaxIDs.txt -taxonomy ../taxonomy/taxonomy.tsv -outdir ../ncbi/format


The major output is a :ref:`sequence tsv file with taxIDs <sequence_tsv_with_taxid_io>`.

Dereplicate NCBI sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Eliminate sequences that are substring of another sequence of the same :ref:`taxID <taxid_glossary>`.

.. code-block:: bash

	perl dereplicate.pl -tsv ../ncbi/format/ncbi_sequences.tsv -outdir ../ncbi/dereplicate -out ncbi_dereplicated_sequences.tsv

The output is the same format as the input tsv file.

.. _bold_sequences_tutorial:

BOLD sequences 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command will download all sequences and lineages for all taxa on the taxon_list from BOLD.
The bold_taxon_list_2022-02-24.txt taxon list file is constructed manually from taxa on `NCBI Taxonomy  <https://www.boldsystems.org/index.php/TaxBrowser_Home>`_. It is available from the data directory of .
Each taxa on the list has less than 500M specimen records on BOLD. 
The taxon_list constructed on 2022-02-24 (bold_taxon_list_2022-02-24.txt) is available with the scripts (data/bold_taxon_list_2022-02-24.txt in `github.com/meglecz/mkCOInr  <https://github.com/meglecz/mkCOInr>`_). This might need to be updated later.

.. code-block:: bash

	perl download_bold.pl -taxon_list ../bold/bold_taxon_list_2022-02-24.txt -outdir ../bold/download -try_download 3


There will be a tsv file for each taxon, where the download was successful. 
The tsv file contains the taxonomic lineage, marker code, sequences and many other information.

NOTE: The download of a long list of taxa takes several days since it is not parallelized. 
You can cut up the input list and run each of them on a separate computers and move the output files to the same folder afterwards.

Format BOLD sequences 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    - Select COI sequences and clean them.
    - Eliminate identical sequences of the same lineage.
    - Clean lineages and make a list with corresponding sequenceIDs.

.. code-block:: bash

	perl format_bold.pl -download_dir ../bold/download/files -outdir ../bold/format


The major output is the following:

-  :ref:`bold_sequences.tsv <sequence_tsv_without_taxid_io>`
-  :ref:`bold_lineages.tsv <lineage_tsv_without_taxid_io>` (all identical lineages are pooled into a same line)


Add taxIDs to BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each lineage this script will:

    - Find an existing :ref:`taxID <taxid_glossary>` at the lowest level possible. TaxIDs can be either from NCBI, or negative taxID already present in :ref:`taxonomy.tsv <taxonomy_io>`.
    - Add new arbitrary (negative) taxIDs to taxa, that are not yet in taxonomy.tsv 
    - Link each new taxID to existing one as a child and include info to the updated taxonomy file
    - Update the input taxonomy file

.. code-block:: bash

	perl add_taxids.pl -lineages ../bold/format/bold_lineages.tsv -sequences ../bold/format/bold_sequences.tsv -taxonomy ../taxonomy/taxonomy.tsv -outdir ../bold/add_taxids


The main output files are the following:

- :ref:`sequences_with_taxIDs.tsv <sequence_tsv_with_taxid_io>`
- :ref:`taxonomy_updated.tsv <taxonomy_io>`


Dereplicate BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Eliminate sequences that are substring of another sequence of the same taxID.

.. code-block:: bash

	perl dereplicate.pl -tsv ../bold/add_taxids/sequences_with_taxIDs.tsv -outdir ../bold/dereplicate -out bold_dereplicated_sequences.tsv


The output is the same format as the input tsv file.

.. _pool_and_dereplicate_tutorial:

Pool and dereplicate datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the dereplicated sequence files from BOLD and NCBI.
Pool the files and dereplicate sequences of a taxID that are present in both files.

.. code-block:: bash

	perl pool_and_dereplicate.pl -tsv1 ../bold/dereplicate/bold_dereplicated_sequences.tsv -tsv2 ../ncbi/dereplicate/ncbi_dereplicated_sequences.tsv -outdir ../COInr -out COInr.tsv


The output is the same format as the input tsv file.

**Move the taxonomy file to the same directory**

.. code-block:: bash

	mv ../bold/add_taxids/taxonomy_updated.tsv ../COInr/taxonomy.tsv

