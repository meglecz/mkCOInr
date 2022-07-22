use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
# tsv file with the following columns:	seqID	taxon_name	sequence

# ALGO
# get lineages using ncbi dmp files (to include synonymes) and if exitst taxonomy file, to include existing arbitrary taxIDs
# 	For all taxa, look for already existing taxid and get the lineage if there is at least one taxid
# 	If more than one taxid creates a lineage for each taxid, and write 1 to the homonomy column
# 	If no taxid for a species, get the taxid for the genus, and get lineage, the species name is kept in the lineage

# check sequnece IDs and write fasta => give a warning for potential pb, but do not exit
#	If some of the sequence IDs are not unique, list the duplicate IDs and exit
#	If sequence IDs format is similar to acessions used in BOLD and NCBI/EMBL/DDBJ => list IDs but continu
#		a fairly safe format is xxx_xxx####, where is x a letter, # is a digit
# clean sequences : gaps deleted, non-TCGA changed to N, external Ns deleted, sequnces with more than $max_n consecutive Ns are deleted


# OUTPUT
# custom_lineages.tsv all identical lineages are pooled into a same line, seqIDs are separated by semi-colon
#	phylum	class	order	family	subfamily	genus	species	seqIDs (with heading)
# custom_sequences.tsv 
#	seqid	sequence (without heading)


# The output file should be ckecked manually to see if the suggested lineages are plausible, and select the correct lineage if there is more than one
# delete the homonymy column for the next step
 

###############################
###############################
my %params = 
(
	'custom' => '', #seqID	taxon	sequence
	'taxonomy' => '', # if sequneces to be included in an already exiting database, where arbitrary taxIDs exist, inlude taxonomy file. Can be empty otherwise.
	'outdir' => '',
	'max_n' => 5, # while clean_seq, eliminate sequences with $max_n or more Ns
	'min_length' => 100, # minimum length of the cleaned sequence; 
	'max_length' => 2000,  # maximum length of the cleaned sequence; 
	'check_seqid_format' => 1 # if 1 check if seq ID resemble to bold and ncbi formats. If yes,  print out warning
);
modify_params_from_tags(\%params, \@ARGV);

my $custom = $params{'custom'};
my $taxonomy = $params{'taxonomy'};
my $outdir = $params{'outdir'};
my $max_n = $params{'max_n'};
my $min_length = $params{'min_length'};
my $max_length = $params{'max_length'};
my $check_seqid_format = $params{'check_seqid_format'};

$outdir = add_slash_to_dir($outdir);
###############################
###############################

#### define filenames
my $t = time;
my $t0 = time;
my $date = get_date();
my $out = $outdir.'custom_lineages.tsv';
my $out_seq = $outdir.'custom_sequences.tsv';
my $log = $outdir.'format_custom.log';
my @ranks = qw(phylum class order family subfamily genus species);
my %stat;

####
unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";


####
#### Read input file 
print "\n####\nRead input file\n";
print LOG "\n####\nRead input file\n";
open(IN, $custom) or die "Cannot open $custom\n";
my $title = <IN>;
my %taxa; #$taxa{taxon}{seqIDs} = ''
my %seq; # $seq{seqID} = seq;
my %non_uniq_seqids; # %non_uniq_seqids{seqid} = '';
while(my $line = <IN>)
{
	++$stat{'1.1 Number of lines in input file: '};
	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split("\t", $line);
	
	$taxa{$line[1]}{$line[0]} = '';
	if(exists $seq{$line[0]})
	{
		$non_uniq_seqids{$line[0]} = '';
	}
	$seq{$line[0]} = $line[2];
}
close IN;

# if some of the sequnece IDs are not unique, print them out and exit
$stat{'1.2 Number of non-uniq seqids: '} = 0;
if(scalar keys %non_uniq_seqids)
{
	$stat{'1.2 Number of non-uniq seqids: '} = scalar keys  %non_uniq_seqids;
	print "\nERROR\nThe following sequence IDs are not unique\nDelete duplicates or modify the sequence IDs\n";
	print join("\n", sort keys  %non_uniq_seqids), "\n";
	print LOG "\nERROR\nThe following sequence IDs are not unique\nDelete duplicates or modify the sequence IDs\n";
	print LOG join("\n", sort keys  %non_uniq_seqids), "\n";
	exit;
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####


####
#### Check format of seqIDs 
print "\n####\nCheck format of seqIDs\n";
print LOG "\n####\nCheck format of seqIDs\n";
if($check_seqid_format)
{
	my %potential_pb = check_seq_id_format(\%seq);
	if(scalar keys %potential_pb)
	{
		$stat{'1.3 Number of numerical seqids: '} = scalar keys  %{$potential_pb{numeric}};
		print "\nWARNING\nNumerical sequence IDs can be confounded with BOLD sequnece IDs\nConsider renaming  the following sequences:\n";
		print join("\n", sort keys  %{$potential_pb{numeric}}), "\n";
		print LOG "\nWARNING\nNumerical sequence IDs can be confounded with BOLD sequnece IDs\nConsider renaming  the following sequences:\n";
		print LOG join("\n", sort keys  %{$potential_pb{numeric}}), "\n";

		$stat{'1.4 Number of ncbi-like seqids: '} = scalar keys  %{$potential_pb{acession}};
		print "\nWARNING\nThe following sequence IDs can be confounded with NCBI/DDBL/EMBL sequnece IDs\nConsider renaming the following sequences:\n";
		print join("\n", sort keys  %{$potential_pb{acession}}), "\n";
		print LOG "\nWARNING\nThe following sequence IDs can be confounded with NCBI/DDBL/EMBL sequnece IDs\nConsider renaming the following sequences:\n";
		print LOG join("\n", sort keys  %{$potential_pb{acession}}), "\n";
	}
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}

####
####


####
#### Check sequences
print "\n####\nCheck sequences\n";
print LOG "\n####\nCheck sequences\n";
$stat{'2.1 Number of sequences that passed the cleaning: '} = 0;
$stat{'2.2 Number of sequences that did not passed the cleaning: '} = 0;
foreach my $sid (sort keys %seq)
{
	my $seq = clean_seq($seq{$sid}, $max_n);
	my $l = length $seq;
	if( $l >= $min_length and $l <= $max_length)
	{
		++$stat{'2.1 Number of sequences that passed the cleaning: '};
		$seq{$sid} = $seq; 
	}
	else
	{
		++$stat{'2.2 Number of sequences that did not passed the cleaning: '};
		delete $seq{$sid}; # delet seqID if seq did not pass the cleaning
	}
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####
####


#### Read taxonomy file
print "\n####\nRead taxonomy\n";
print LOG "\n####\nRead taxonomy\n";

my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
read_taxonomy_to_tax_hash(\%tax, $taxonomy);

my %name_taxids; # $name_taxids{name}{$taxid} = '' # All names, including synonyms and homonyms
my %taxid_names; # $taxid_names{taxid}{names} = '' # All names, including synonyms and homonyms
my %name_taxids_par; #$name_taxids_par{$name}{taxid} = parent taxid # only scientifi names
read_taxonomy_names(\%taxid_names, \%name_taxids, \%name_taxids_par, $taxonomy);
%taxid_names = (); # delete hash. It is produced automatically by the routine, but we do not need it.
%name_taxids_par = (); # delete hash. It is produced automatically by the routine, but we do not need it.

my %taxid_full_lineage_id; # $taxid_full_lineage_id{taxid} = (list of taxids starting from the most distant)
foreach my $taxid (keys %tax)
{
	@{$taxid_full_lineage_id{$taxid}} = get_lineage_list($taxid, \%tax, 0);
}

my %tax_ranked_lineage; # $tax_ranked_lineage{taxid} = (scientific name, species, genus, family, order, class, phylum, kingdom);
foreach my $taxid (keys %tax)
{
	#@lin = ('root', '', '','','','','','',''); # root + 8 taxlevel
	my @lin = make_ranked_lineage($taxid_full_lineage_id{$taxid}, \%tax);
	shift(@lin); # delete root
	@lin = reverse@lin;
	unshift(@lin, $tax{$taxid}[2]);
	@{$tax_ranked_lineage{$taxid}} = @lin;
}
%taxid_full_lineage_id = ();  # delete hash.We dod not need it anymore.

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####
####





####
#### Get lineages
print "\n####\nGet lineages\n";
print LOG "\n####\nGet lineages\n";
#$taxa{taxon}{seqIDs} = ''
my %lin; # %lin{$taxon} = (lineage1, lineage2)
$stat{'3.1 Number of unique taxa in input file: '} = scalar keys %taxa;
$stat{'3.2 Number of taxa without valid sequences : '} = 0;
$stat{'3.3 Number of taxa with lineage: '} = 0;
$stat{'3.4 Number of taxa without lineage: '} = 0;
foreach my $taxon (sort keys %taxa)
{
	my $taxon_orig = $taxon;
	foreach my $sid (keys %{$taxa{$taxon}}) # delete seqids that has been eliminated in the sequence cleaning step
	{
		unless(exists $seq{$sid})
		{
			delete $taxa{$taxon}{$sid};
		}
	}
	
	if(scalar keys %{$taxa{$taxon}} == 0)
	{
		delete $taxa{$taxon};
		++$stat{'3.2 Number of taxa without valid sequences : '};
	}
	else # valid sequences for the taxon => get lineage
	{
		if($taxon =~ /([A-Z][a-z-]+ [a-z-]+) [a-z-]+$/) # subspecies =>  change to species
		{
			$taxon = $1;
		}
		my $homonym = 0;
		# get lineage
		if(exists $name_taxids{$taxon}) # taxid for taxon
		{
			++$stat{'3.3 Number of taxa with lineage: '};
			if(scalar keys %{$name_taxids{$taxon}} > 1)
			{
				$homonym = 1;
			}
			foreach my $taxid (keys %{$name_taxids{$taxon}})
			{
				my $lineage = get_linenage(\%tax_ranked_lineage, \%tax, \@ranks, $taxid); #phylum	class	order	family	subfamily	genus	species
				$lineage .= "\t".$homonym;
				push(@{$lin{$taxon_orig}}, $lineage); # If homonymy there can be more than one lineages
			}
		}
		elsif($taxon =~ /([A-Z-][a-z-]+) [a-z-]+$/ and exists $name_taxids{$1}) # try genus name is no taxid for the taxon that has a latin species name format
		{
			++$stat{'3.3 Number of taxa with lineage: '};
			$taxon = $1;
			if(scalar keys %{$name_taxids{$taxon}} > 1)
			{
				$homonym = 1;
			}
			foreach my $taxid (keys %{$name_taxids{$taxon}})
			{
				my $lineage = get_linenage(\%tax_ranked_lineage, \%tax, \@ranks, $taxid);
				$lineage =~ s/\t[^\t]*$/\t$taxon_orig/; # replace the shortenned name to the original taxon name
				$lineage .= "\t".$homonym;
				push(@{$lin{$taxon_orig}}, $lineage); # If homonymy there can be more than one lineages
			}
		}
		else
		{
			++$stat{'3.4 Number of taxa without lineage: '};
			my $lineage = "NA\t" x 6 .$taxon_orig. "\t0";
			push(@{$lin{$taxon_orig}}, $lineage);
		}
	}
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####



#### Print outfiles
print "\n####\nPrint outfiles\n";
print LOG "\n####\nPrint outfiles\n";
#### lineage file
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "phylum	class	order	family	subfamily	genus	species	homonymy	seqIDs\n";
foreach my $taxon (sort keys %lin)
{
	my @sids = sort keys %{$taxa{$taxon}};
	foreach my $lineage (@{$lin{$taxon}}) # If homonymy there can be more than one lineages
	{
		print OUT $lineage, "\t", join(';', @sids), "\n";
	}
}
close OUT;

#### sequence file
open(OUT, '>', $out_seq) or die "Cannot open $out_seq\n";
print OUT "seqID	sequence\n";
foreach my $sid (sort keys %seq)
{
	print OUT "$sid\t$seq{$sid}\n";
}
close OUT;
#### 

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
print LOG print_stat(\%stat, $t0);
close LOG;

exit;


#############################################

sub check_seq_id_format
{
	my ($seq_hash) = @_;
	
	my %pb;
	foreach my $id (keys %$seq_hash)
	{
		if($id =~ /^[0-9]+$/)
		{
			$pb{numeric}{$id} = '';
		}
		elsif($id =~ /^[A-Z]{1,2}_*[0-9]{5,8}$/)
		{
			$pb{acession}{$id} = '';
		}
	}
	return %pb;
}

##############################################"
sub get_linenage
{
	my ($tax_ranked_lineage, $tax, $ranks, $taxid) = @_;
#tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	
	my @ncbi_lin = @{$$tax_ranked_lineage{$taxid}};
	my @lin = ($ncbi_lin[6], $ncbi_lin[5], $ncbi_lin[4], $ncbi_lin[3], '', $ncbi_lin[2], $ncbi_lin[1]);
	
	for(my $i = 0; $i < scalar @$ranks; ++$i)# put the scientific name to the correct taxlevel
	{
		if($$tax{$taxid}[1] eq $$ranks[$i])
		{
			$lin[$i] = $ncbi_lin[0];
			last;
		}
	}
	return join("\t", @lin);
}

sub print_help
{

##############################################"
print '
usage: perl format_custom.pl [-options] -custom INPUT_TSV_FILE -taxonomy TAXONOMY_TSV -outdir OUTDIR

 ARGUMENTS
   -custom                 Input tsv file: seqID,taxon_name,sequence
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
   -outdir                 Name of the otput directory
 OPTIONS
   -max_n                  Positive integer; Default:5)
                              Eliminates sequences with max_n or more consecutive Ns
   -min_length             Positive integer; Default:100
                              Minimum length of the cleaned sequence;
   -max_length             Positive integer; Default:2000
                              Maximum length of the cleaned sequence;
   -check_seqid_format     [0/1]; Default: 1
                              If 1, print out warnings if seqIDs are similar to bold and ncbi seqID formats
',  "\n";
  exit;
}
