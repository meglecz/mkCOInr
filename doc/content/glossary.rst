Glossary
=======================================

.. _ranked_lineage_glossary:

ranked lineage
------------------------------------------------

Taxonomic lineage of a taxon, containing the major taxonomic ranks: domain, kingdom, phylum, class, order, family, genus, species

.. _seqid_glossary:

seqID
------------------------------------------------

Short and unique ID for each sequence. 

They can be NCBI accession numbers, numerical values derived from BOLD, or arbitrary values for custom sequences.  

For arbitrary values, use only alphanumerical characters and underscore (_). 
To avoid potential duplicate seqIDs for custom sequences, I suggest the following format:  xxx_xxx####, where is x a letter, # is a digit.

.. _taxid_glossary:

taxID
------------------------------------------------

Numerical ID specific for each taxon. 

Positive integer values are from NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy). 
Negative values are used for taxa not present in NCBI taxonomy.

Each taxID is linked to a direct parent taxID, a scientific name, a taxonomic rank and a :ref:`taxlevel_index<taxlevel_index_glossary>`

.. _taxlevel_index_glossary:

taxlevel index
------------------------------------------------

Integer associated to a major taxonomic rank.

0 => root, 1=> domain, 2=> kingdom, 3=> phylum, 4=> class, 5=> order, 6=> family, 7=> genus, 8=> species

Levels in between have 0.5 added to the next highest level (e.g. 5.5 for infraorder and for superfamily).


