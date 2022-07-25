.. _installation_installation:

Installation
=================================================

Special attention was taken to minimize dependencies. 

Third-party programs
-------------------------------------------------

The following tools should be installed. (Scripts using the program are in parantheses)

    - `nsdpy <https://github.com/RaphaelHebert/nsdpy>`_ used for downloading sequences from NCBI (:ref:`Hebert and Meglécz, 2022<Hebert_2022_reflist>`)
   
    - `BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ (:ref:`Altschul et al., 1997<Altschul_1997_reflist>`; :ref:`format_bold.pl <format_bold_reference>`, for the check_orientation option, :ref:`format_db.pl <format_db_reference>` for the blast option option)
   
    - `vsearch <https://github.com/torognes/vsearch>`_ (:ref:`Rognes et al., 2016<Rognes_2016_reflist>`; :ref:`dereplicate.pl <dereplicate_reference>`, :ref:`pool_and_dereplicate.pl <pool_and_dereplicate_reference>`, :ref:`select_region.pl <select_region_reference>`)
   
    - `cutadapt <https://cutadapt.readthedocs.io>`_ (:ref:`Martin, 2011<Martin_2011_reflist>`; :ref:`select_region.pl <select_region_reference>` for the e_pcr option)
   

All third party programs can be easily installed to a conda environment, but it is not essential.

Commands for a quick installation of the conda environment and dependencies

.. code-block:: bash

	conda create --name mkcoinr python=3.9 -y
	conda activate mkcoinr

	python3 -m pip install cutadapt
	conda install -c bioconda blast -y
	conda install -c bioconda vsearch -y
	pip install nsdpy

.. _mkCOInr scripts_installation:

mkCOInr scripts
-------------------------------------------------

The mkCOInr scripts are written in Perl, no installation is necessary to run them apart from the Perl interpreter already installed in all unix systems. 
Just clone the :ref:`mkCOInr git repository <https://github.com/meglecz/mkCOInr>`_ and you are ready.


