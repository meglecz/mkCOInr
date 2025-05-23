.. _tutorial_tutorial:

Tutorial
============

After installing mkCOInr you have a file system like this:

.. code-block:: bash

	mkCOInr
	├── data
	│   ├── bold_taxon_list_2022-02-24.txt
	│   ├── example
	│   │   ├── custom_lineages_verified.tsv
	│   │   ├── my_sequences.tsv
	│   │   ├── taxon_list_eukaryota.tsv
	│   │   ├── taxon_list_insecta.tsv
	│   │   └── taxon_list.tsv
	│   └── one_seq_per_order_658.fas
	├── doc
	... (abbreviated)
	└── scripts
		├── add_taxids.pl
		├── dereplicate.pl
		├── download_bold.pl
		├── download_taxonomy.pl
		├── format_bold.pl
		├── format_custom.pl
		├── format_db.pl
		├── format_ncbi.pl
		├── get_subtaxa.pl
		├── mkdb.pm
		├── pool_and_dereplicate.pl
		├── select_region.pl
		└── select_taxa.pl


.. _customize_tutorial:

Customize database 
-------------------------------------------------

In the first part of the tutorial, I will start from the COInr database in each major step to illustrate 
    - :ref:`How to include custom sequences <add_custom_sequences_tutorial>`
    - :ref:`How to select or eliminate sequences of a list of taxa or a minimum resolution <select_sequences_custom_tutorial>`
    - :ref:`How to select a target region <select_region_custom_tutorial>`
    - :ref:`How to format a dataset to different database formats <format_db_custom_tutorial>`
    
These steps can be executed independently. 

The last example shows how to :ref:`create a pipeline <chained_custom_tutorial>` by combining different commands.

The creation of the COInr database is explained in the :ref:`Create COInr from BOLD and NCBI section <create_coinr_tutorial>`. 
You can download this database from `Zenodo <https://doi.org/10.5281/zenodo.6555985>`_ and customize it to your needs.

.. _download_coinr_tutorial:

Download and untar COInr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will need to change the date in the filename, and get the up-to-date link from zenodo for later releases.

.. code-block:: bash

	cd mkCOInr
	wget https://zenodo.org/record/6555985/files/COInr_2022_05_06.tar.gz
	tar -zxvf COInr_2022_05_06.tar.gz
	rm COInr_2022_05_06.tar.gz


For shortening the paths in this tutorial, rename COInr_2022_05_06 directory to COInr.

.. code-block:: bash

	mv COInr_2022_05_06 COInr

This gives the following file structure

.. code-block:: bash

	mkCOInr
	├── COInr
	│   ├── COInr.tsv
	│   └── taxonomy.tsv
	├── data
	│   ├── bold_taxon_list_2022-02-24.txt
	│   ├── example
	│   │   ├── custom_lineages_verified.tsv
	│   │   ├── my_sequences.tsv
	│   │   ├── taxon_list_eukaryota.tsv
	│   │   ├── taxon_list_insecta.tsv
	│   │   └── taxon_list.tsv
	│   └── one_seq_per_order_658.fas
	...(abbreviated)
	└── scripts
		├── add_taxids.pl
		├── dereplicate.pl
		...(abbreviated)


The COInr database is composed of two files
    - :ref:`COInr.tsv <sequence_tsv_with_taxid_io>` contains :ref:`sequenceIDs <seqid_glossary>`, :ref:`taxIDs <taxid_glossary>` and sequences
    - :ref:`taxonomy.tsv <taxonomy_io>` contains all taxIDs and associated information


.. _add_custom_sequences_tutorial:

Add custom sequences to a database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _format_custom_tutorial:

Format custom files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`input tsv file <custom_sequences_tsv_io>` (-custom) contains :ref:`seqIDs <seqid_glossary>`, 
taxon names (can be at any taxonomic level) and sequences (see the example data/example/my_sequences.tsv).
The :ref:`format_custom.pl <format_custom_reference>` script will suggest one or more lineages for each 
taxon name based on the existing lineages in :ref:`taxonomy.tsv <taxonomy_io>`. It will also consider synonyms.

.. code-block:: bash

	perl scripts/format_custom.pl -custom data/example/my_sequences.tsv -taxonomy COInr/taxonomy.tsv -outdir tutorial/custom/1_format


The output lineage file (custom_lineages.tsv) looks like this:

.. code-block:: bash

	phylum	class	order	family	subfamily	genus	species	homonymy	seqIDs
	Mollusca	Bivalvia	Cardiida	Cardiidae		Acanthocardia	Acanthocardia paucicostata	0	Seq113;Seq88
	NA	NA	NA	NA	NA	NA	Ilia nucleus	0	Seq117
	Streptophyta	Magnoliopsida	Ericales	Ericaceae		Leucothoe		1	Seq96
	Arthropoda	Malacostraca	Amphipoda	Leucothoidae		Leucothoe		1	Seq96
	Annelida	Polychaeta	Phyllodocida	Polynoidae				0	Seq65


This output should should be checked manually to see if the lineages are coherent.
If homonymy, choose the correct lineage (e.g. for *Leucothoe* genus), then delete homonymy column. 

If a taxon name is not present in the taxonomy file, the lineage should be completed manually (e.g. *Ilia nucleus* in the example file).

I created a revised version of the lineage file (data/example/custom_lineages_verified.tsv), which will be used in the next step:

.. code-block:: bash

	phylum	class	order	family	subfamily	genus	species	seqIDs
	Mollusca	Bivalvia	Cardiida	Cardiidae		Acanthocardia	Acanthocardia paucicostata	Seq113;Seq88
	Arthropoda	Malacostraca	Decapoda	Leucosiidae		Ilia	Ilia nucleus	Seq117
	Arthropoda	Malacostraca	Amphipoda	Leucothoidae		Leucothoe		Seq96
	Annelida	Polychaeta	Phyllodocida	Polynoidae				Seq65

See details in description section: :ref:`format_custom.pl <format_custom_reference>` script.



.. _add_taxids_custom_tutorial:

Add taxIDs to custom sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`add_taxids.pl <add_taxids_reference>` script will

- For each lineage in the input file
    - Find an existing taxID at the lowest possible taxonomic level. :ref:`taxIDs <taxid_glossary>` can be either from NCBI, or negative taxID already present in :ref:`taxonomy.tsv <taxonomy_io>`.
    - Add new arbitrary (negative) taxIDs to taxa not yet in the taxonomy file 
    - Link each new taxID to an existing one as a child and include info to the updated taxonomy file
- Make a :ref:`tsv file with sequences and taxIDs <sequence_tsv_with_taxid_io>`
- Update the :ref:`taxonomy.tsv <taxonomy_io>` file

.. code-block:: bash

	perl scripts/add_taxids.pl -lineages data/example/custom_lineages_verified.tsv -sequences tutorial/custom/1_format/custom_sequences.tsv -taxonomy COInr/taxonomy.tsv -outdir tutorial/custom/2_add_taxids

See details in description section: :ref:`add_taxids.pl <add_taxids_reference>` script.


.. _dereplicate_custom_tutorial:

Dereplicate custom sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`dereplicate.pl <dereplicate_reference>` script will eliminate sequences that are substrings of another sequence of the same :ref:`taxID <taxid_glossary>`. 
Use :ref:`sequences_with_taxIDs.tsv <sequence_tsv_with_taxid_io>` file (output of the previous script) as the input.

.. code-block:: bash

	perl scripts/dereplicate.pl -tsv tutorial/custom/2_add_taxids/sequences_with_taxIDs.tsv -outdir tutorial/custom/3_dereplicate -out custom_dereplicated_sequences.tsv

The output file is in the same format as the input tsv file.

See details in description section: :ref:`dereplicate.pl <dereplicate_reference>` script.


.. _pool_and_dereplicate_custom_tutorial:

Pool and dereplicate datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use two dereplicated :ref:`sequence tsv files <sequence_tsv_with_taxid_io>`:
    - COInr.tsv  (pool of BOLD and NCBI, downloaded from Zenodo)
    - custom_dereplicated_sequences.tsv (output of the previous script)

:ref:`pool_and_dereplicate.pl <pool_and_dereplicate_reference>` will pool the files and dereplicate sequences 
of the taxIDs that are present in both files.

.. code-block:: bash

	perl scripts/pool_and_dereplicate.pl -tsv1 COInr/COInr.tsv -tsv2 tutorial/custom/3_dereplicate/custom_dereplicated_sequences.tsv -outdir tutorial/custom -out COInr_custom.tsv

The output is the same format as the input tsv file.

See details in description section: :ref:`pool_and_dereplicate.pl <pool_and_dereplicate_reference>` script.


Custom database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your custom database is composed of two files:
    - the dereplicated sequence file (COInr_custom.tsv)
    - the last version of the taxonomy file (taxonomy_updated.tsv)

For simplicity, move the updated taxonomy file to the same folder as the sequence file.

.. code-block:: bash

	mv tutorial/custom/2_add_taxids/taxonomy_updated.tsv tutorial/custom/taxonomy_updated.tsv


This database can be further customized, or you can simply be formated to your taxonomic assignment program by the :ref:`format_db.pl <format_db_reference>` script.




.. _select_sequences_custom_tutorial:

Select sequences from existing database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select sequences for a list of taxa with a minimum taxonomic rank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequences can be selected for a list of taxa and/or for a minimum taxonomic level (species/genus/family/order/class/phylum/kingdom/domain/root)

The input file (:ref:`-taxon_list <taxon_list_io>`) contains a list of taxa and eventually their taxIDs (see example data/example/taxon_list.tsv). 

.. code-block:: bash

	perl scripts/select_taxa.pl -taxon_list data/example/taxon_list.tsv -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv  -min_taxlevel species  -outdir tutorial/select_taxa_0 -out COInr_selected.tsv

The main output is a :ref:`sequence tsv file <sequence_tsv_with_taxid_io>` (COInr_selected.tsv).
A :ref:`lineage file <lineage_tsv_with_taxID_io>` (taxa_with_lineages.tsv) is also written for all taxa in the taxon_list to check if they are coherent with the target taxon names. 

See details in description section: :ref:`select_taxa.pl <select_taxa_reference>` script.




Excluding  sequences of a taxon list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the same script it is also possible to eliminate sequences of taxa instead of selecting them. Set the *negative_list* option to 1 to do that.

.. code-block:: bash

	perl scripts/select_taxa.pl -taxon_list data/example/taxon_list.tsv -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv  -min_taxlevel species -outdir tutorial/select_taxa_1 -out COInr_reduced.tsv -negative_list 1

See details in description section: :ref:`select_taxa.pl <select_taxa_reference>` script.




.. _select_region_custom_tutorial:

Select region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequences can be trimmed to a specific region of the COI gene by the :ref:`select_region.pl <select_region_reference>` script. 
To define the region, you can either give a fasta file with sequences trimmed to the region of interest, 
or you can detect it automatically by e-pcr.


.. _select_region_e_pcr_custom_tutorial:

Select region using the e_pcr option
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primers used in this example are amplifying a Leray fragment (ca. 313 bp of the second half of the barcode region).

.. code-block:: bash

	perl scripts/select_region.pl -tsv COInr/COInr.tsv -outdir tutorial/select_region/ePCR -e_pcr 1 -fw GGNTGAACNGTNTAYCCNCC -rv TAWACTTCDGGRTGNCCRAARAAYCA -trim_error 0.3 -min_amplicon_length 280 -max_amplicon_length 345 -min_overlap 20 -tcov 0.8 -identity 0.7


.. _select_region_bait_fas_custom_tutorial:

Select region using the bait_fas option
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the *e_pcr* option is an easy way to produce some sequences trimmed to the target region, 
and they can be used as a database to align all other sequences to them. 
However, if the parameters of the e_pcr are relaxed, it can produce some false positives. 
An alternative solution is to use a small, taxonomically divers fasta file, with sequences already trimmed to the target region 
(-*bait_fas* option). 
An example of such a file is given in the data directory (data/one_seq_per_order_658.fas). 
It contains one sequence for each taxonomic order among the taxa that have a compete mitochondrial genome available in GenBank. 
Sequences are trimmed to the approximately 658 bp (depending on the taxon) barcode fragment of the COI gene.

.. code-block:: bash

	perl scripts/select_region.pl -tsv COInr/COInr.tsv -outdir tutorial/select_region/bait_fas -e_pcr 0 -bait_fas data/one_seq_per_order_658.fas -tcov 0.8 -identity 0.7


See details in description section: :ref:`select_region.pl <select_region_reference>` script.



.. _format_db_custom_tutorial:

Format database 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format the database to one of the following formats
    - qiime
    - rdp
    - full
    - blast
    - vtam
    - sintax

**qiime**

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -outfmt qiime -outdir COInr/qiime -out COInr_qiime


**rdp**

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -outfmt rdp -outdir COInr/rdp -out COInr_rdp

You should use the rdp_calssifier or qiime's feature-classifier to train the database using the output files of this script if you have used the rdp or qiime options.


**full**

The full option, gives a :ref:`tsv file <full_tsv_io>` with seqIDs, ranked lineages, taxIDs for each sequence, and this is a very easy-to-parse, complete file.

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -outfmt full -outdir COInr/full -out COInr_full

**sintax**

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -outfmt sintax -outdir COInr/sintax -out COInr_sintax


**blast**

For making a BLAST database, the taxonomy file is not necessary and the indexed files in the output folder are ready to use.

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -outfmt blast -outdir COInr/blast -out COInr_blast

**vtam**

The vtam option produces a BLAST database and a taxonomy file adapted to `VTAM <https://github.com/aitgon/vtam>`_ .

.. code-block:: bash

	perl scripts/format_db.pl -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -outfmt vtam -outdir COInr/vtam -out COInr_vtam

See details in description section: :ref:`format_db.pl <format_db_reference>` script.


.. _chained_custom_tutorial:

Chaining steps to make a custom database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the above examples, we have started from the COInr database. However, you can chain the different commands. 

Bellow, I will show you how to create a database with the following characteristics:
    - Eukaryota sequences
    - Excluding insects
    - Enriched with custom sequences
    - Sequences assigned at least to genus level
    - Trimmed to the Leray fragment (ca. 313 nt of the second half of the barcode region) of the COI gene (keep sequences if cover at least 90% of the target region)
    - rdp_classifier format


**Notes**:
    - It is a good idea to start with steps that are relatively quick and reduce the size of the database. 
    - Since, over 70% of the sequences are from Insecta in COInr, we will start by eliminating them. 
    - The custom sequences are all Non-Insect Eukaryotes, so we can add custom sequences to the reduced dataset. Otherwise, we should have started by adding custom sequences. This solution is also fine, but gives large intermediate files.
    - The selection of the target region is the most computationally intensive, and the more diverse the dataset, the less precise it is. So it is preferable to do this at the end of the pipeline.

.. _exclude_insecta_tutorial:

Exclude Insecta and sequences with resolution lower than genus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	perl scripts/select_taxa.pl -taxon_list data/example/taxon_list_insecta.tsv -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv  -min_taxlevel genus -outdir tutorial/chained/1_noInsecta -out COInr_noIns.tsv -negative_list 1


.. _keep_eukaryota_tutorial:

Keep Eukaryota
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	perl scripts/select_taxa.pl -taxon_list data/example/taxon_list_eukaryota.tsv -tsv tutorial/chained/1_noInsecta/COInr_noIns.tsv -taxonomy COInr/taxonomy.tsv -outdir tutorial/chained/2_Eukaryota -out COInr_noIns_Euk.tsv


.. _add_custom_chained_tutorial:

Add custom sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	perl scripts/format_custom.pl -custom data/example/my_sequences.tsv -taxonomy COInr/taxonomy.tsv -outdir tutorial/chained/3_add_custom/1_format

Check and format the custom_lineages.tsv and make custom_lineages_verified.tsv as in :ref:`Add custom sequences to a database <add_custom_sequences_tutorial>` section.

.. code-block:: bash

	perl scripts/add_taxids.pl -lineages data/example/custom_lineages_verified.tsv -sequences tutorial/chained/3_add_custom/1_format/custom_sequences.tsv -taxonomy COInr/taxonomy.tsv -outdir tutorial/chained/3_add_custom/2_add_taxids
	
	perl scripts/dereplicate.pl -tsv tutorial/chained/3_add_custom/2_add_taxids/sequences_with_taxIDs.tsv -outdir tutorial/chained/3_add_custom/3_dereplicate -out custom_dereplicated_sequences.tsv

Add the formatted, dereplicated custom sequences to the sequences in COInr_noIns_Euk.tsv

.. code-block:: bash

	perl scripts/pool_and_dereplicate.pl -tsv1 tutorial/chained/2_Eukaryota/COInr_noIns_Euk.tsv -tsv2 tutorial/chained/3_add_custom/3_dereplicate/custom_dereplicated_sequences.tsv -outdir tutorial/chained/3_add_custom -out COInr_noIns_Euk_custom.tsv
	
	mv tutorial/chained/3_add_custom/2_add_taxids/taxonomy_updated.tsv tutorial/chained/3_add_custom/taxonomy_updated.tsv


.. _keep_genus_tutorial:

Keep only sequences with genus or higher resolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have eliminated sequences with lower than genus resolution from COInr in the first step (-min_taxlevel genus). 
However, among the custom sequences we had a sequence with an unknown genus. 
So let's redo the selection for a minimum taxonomic level. 

Yes, you are right! We could have just avoided adding that sequence to the database in the previous step. 
But if you have many custom sequences, you might just be lazy to check the custom sequences manually, 
and in that case you can use mkCOInr to this for you.

**Attention**: From now on, we have to use the updated taxonomy file (taxonomy_updated.tsv), since some of the taxa of the custom sequences might not be in the original taxonomy.tsv file.

.. code-block:: bash

	perl scripts/select_taxa.pl -tsv tutorial/chained/3_add_custom/COInr_noIns_Euk_custom.tsv -taxonomy tutorial/chained/3_add_custom/taxonomy_updated.tsv -outdir tutorial/chained/4_genus -out COInr_noIns_Euk_custom_genus.tsv


.. _trim_to_leray_tutorial:

Trim to Leray region
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	perl scripts/select_region.pl -tsv tutorial/chained/4_genus/COInr_noIns_Euk_custom_genus.tsv -outdir tutorial/chained/5_select_region -e_pcr 1 -fw GGNTGAACNGTNTAYCCNCC -rv TAWACTTCDGGRTGNCCRAARAAYCA -trim_error 0.3 -min_amplicon_length 280 -max_amplicon_length 345 -min_overlap 20 -tcov 0.9 -identity 0.7



.. _format_rdp_chained_tutorial:

Format for RDP_classifier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	perl scripts/format_db.pl -tsv tutorial/chained/5_select_region/trimmed.tsv -taxonomy tutorial/chained/3_add_custom/taxonomy_updated.tsv -outfmt rdp -outdir tutorial/chained/6_rdp -out COInr_customized





.. _create_coinr_tutorial:

Create COInr from BOLD and NCBI
-------------------------------------------------
The following steps describe how COInr database (available at `Zenodo <https://doi.org/10.5281/zenodo.6555985>`_ ) was produced. 

.. _download_ncbi_taxonomy_tutorial:

Download NCBI taxonomy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download NCBI taxonomy dmp file and create :ref:`taxonomy.tsv <taxonomy_io>`.

.. code-block:: bash

	cd mkCOInr
	perl scripts/download_taxonomy.pl -outdir COInr_new/taxonomy

See details in description section: :ref:`download_taxonomy.pl <download_taxonomy_reference>` script.

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

The sequences.fasta file contains all CDS sequences. Sequences are correctly oriented but should still be filtered to keep only COI sequences.
TaxIDs.txt contains the sequenceIDs and the TaxIDs.

Move the results of nsdpy to the COInr_new/ncbi/download directory and clean up the directory.

.. code-block:: bash

	mkdir -p COInr_new/ncbi
	mv NSDPY_results/yyyy-mm-dd_hh-mm-ss COInr_new/ncbi/download
	mv report.tsv COInr_new/ncbi/download
	rmdir NSDPY_results


Format NCBI sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`format_ncbi.pl <format_ncbi_reference>` script will
    - Select COI sequences and clean them. 
    - Eliminate identical sequences of the same taxID.
    - Clean tax names and taxids.

.. code-block:: bash

	perl scripts/format_ncbi.pl -cds COInr_new/ncbi/download/sequences.fasta -taxids COInr_new/ncbi/download/TaxIDs.txt -taxonomy COInr_new/taxonomy/taxonomy.tsv -outdir COInr_new/ncbi/format

The major output is a :ref:`sequence tsv file with taxIDs <sequence_tsv_with_taxid_io>`.

See details in description section: :ref:`format_ncbi.pl <format_ncbi_reference>` script.

Dereplicate NCBI sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Eliminate sequences that are substring of another sequence of the same :ref:`taxID <taxid_glossary>`.

.. code-block:: bash

	perl scripts/dereplicate.pl -tsv COInr_new/ncbi/format/ncbi_sequences.tsv -outdir COInr_new/ncbi/dereplicate -out ncbi_dereplicated_sequences.tsv

The output is the same format as the input tsv file.

See details in description section: :ref:`dereplicate.pl <dereplicate_reference>` script.

.. _bold_sequences_tutorial:

BOLD sequences 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`download_bold.pl <download_bold_reference>` script is deprecated. The BOLD API used in download_bold.pl do not allow anymore to download large data files.

It is possible, however, to download all public sequences as a data package from 
`https://www.boldsystems.org/index.php/datapackages <https://www.boldsystems.org/index.php/datapackages>`_. 
You need to have a BOLD account for downloading the data package in (tar.gz compressed) format, 
that contains a TSV file with sequences, taxonomic lineages and other metadata. This uncompressed file will be the input of format_bold_package.pl.


Format BOLD sequences 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`format_bold_package.pl <format_bold_package_reference>` script will
    - Select COI sequences and clean them
    - Select sequences with out without BIN_URI according to the delete_noBIN argument
    - Eliminate identical sequences of the same lineage
    - Clean lineages and make a list with corresponding sequenceIDs

.. code-block:: bash

	perl scripts/format_bold_package.pl -bold_data COInr_new/bold/download/BOLD_Public.26-Apr-2024.tsv -outdir COInr_new/bold/format -delete_noBIN 1

The major output is the following:
    - :ref:`bold_sequences.tsv <sequence_tsv_without_taxid_io>`
    - :ref:`bold_lineages.tsv <lineage_tsv_without_taxid_io>` (all identical lineages are pooled into a same line)

See details in description section: :ref:`format_bold_package.pl <format_bold_package_reference>` script.


Add taxIDs to BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each lineage the :ref:`add_taxids.pl <add_taxids_reference>` script will
    - Find an existing :ref:`taxID <taxid_glossary>` at the lowest level possible. TaxIDs can be either from NCBI, or negative taxID already present in :ref:`taxonomy.tsv <taxonomy_io>`.
    - Add new arbitrary (negative) taxIDs to taxa, that are not yet in taxonomy.tsv 
    - Link each new taxID to existing one as a child and include info to the updated taxonomy file
    - Update the input taxonomy file

.. code-block:: bash

	perl scripts/add_taxids.pl -lineages COInr_new/bold/format/bold_lineages.tsv -sequences COInr_new/bold/format/bold_sequences.tsv -taxonomy COInr_new/taxonomy/taxonomy.tsv -outdir COInr_new/bold/add_taxids

The main output files are the following:
    - :ref:`sequences_with_taxIDs.tsv <sequence_tsv_with_taxid_io>`
    - :ref:`taxonomy_updated.tsv <taxonomy_io>`

See details in description section: :ref:`add_taxids.pl <add_taxids_reference>` script.

Dereplicate BOLD sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Eliminate sequences that are substring of another sequence of the same taxID.

.. code-block:: bash

	perl scripts/dereplicate.pl -tsv COInr_new/bold/add_taxids/sequences_with_taxIDs.tsv -outdir COInr_new/bold/dereplicate -out bold_dereplicated_sequences.tsv

The output is the same format as the input tsv file.

See details in description section: :ref:`dereplicate.pl <dereplicate_reference>` script.


.. _pool_and_dereplicate_tutorial:

Pool and dereplicate datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the dereplicated sequence files from BOLD and NCBI.
The :ref:`pool_and_dereplicate.pl <pool_and_dereplicate_reference>` script will pool the files and dereplicate sequences of a taxID that are present in both files.

.. code-block:: bash

	perl scripts/pool_and_dereplicate.pl -tsv1 COInr_new/bold/dereplicate/bold_dereplicated_sequences.tsv -tsv2 COInr_new/ncbi/dereplicate/ncbi_dereplicated_sequences.tsv -outdir COInr_new -out COInr.tsv

The output is the same format as the input tsv file.

See details in description section: :ref:`pool_and_dereplicate.pl <pool_and_dereplicate_reference>` script.

**Move the taxonomy file to the same directory**

.. code-block:: bash

	mv COInr_new/bold/add_taxids/taxonomy_updated.tsv COInr_new/taxonomy.tsv

