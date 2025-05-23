use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
# -taxon_name OR taxid
# -taxonomy file (tax_id	parent_tax_id	rank	name_txt	old_tax_id	synonyms)

# AIM
# make a list of sub-taxa of the input taxon at the next major taxonomic level

# If taxon is a taxon name, get all taxIDs that correspond to this name (e.g. 1266065 and 50622 for Plecoptera)
# Determine the next lowest major taxonomic rank (phylum, class, order, family, genus, species) for each taxID 
# (e.g. if taxId is an order or suborder or superfamily, the next major tax rank is family)
# List subtaxa of each taxID of this taxonomic rank.

# OUTPUT
# tsv file with the list of sub-taxon names and taxids

###############################
###############################
my %params = (
'taxon' =>  '', # taxon_name OR taxid
'taxonomy' =>  '', 
'outdir' => ''
);

modify_params_from_tags(\%params, \@ARGV);

my $taxon = $params{taxon};
my $taxonomy = $params{taxonomy};
my $outdir = $params{outdir};

$outdir = add_slash_to_dir($outdir);

###############################
###############################


#### define filenames and variables
my $t = time;
my $t0 = time;
my $date = get_date();
unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}
my $log = $outdir.'get_subtaxa.log';
my %stat;
my %taxlevel = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'domain',1, '', 0, 'root', 0);

open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";
###

### Read taxonomy
print "\n####\nRead taxonomy\n";
print LOG "\n####\nRead taxonomy\n";
my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
read_taxonomy_to_tax_hash(\%tax, $taxonomy);

my %name_taxids; # $name_taxids{name}{$taxid} = '' # All names, including synonyms and homonyms
my %taxid_names; # $taxid_names{taxid}{names} = '' # All names, including synonyms and homonyms
my %name_taxids_par; #$name_taxids_par{$name}{taxid} = parent taxid # only scientific names
read_taxonomy_names(\%taxid_names, \%name_taxids, \%name_taxids_par, $taxonomy);
%taxid_names = (); # delete hash. It is produced automatically by the routine, but we do not need it.
%name_taxids_par = (); # delete hash. It is produced automatically by the routine, but we do not need it.
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

### Get sub_taxa
print "\n####\nGet the list of subtaxa at the next major taxonomic rank\n";
print LOG "\n####\nGet the list of subtaxa at the next major taxonomic rank\n";
my @taxids;
if($taxon =~ /^[-0-9]+$/) # taxid is given
{
	$taxids[0] = $taxon;
	unless(exists $tax{$taxids[0]})
	{
		print "$taxids[0] is not in the list of taxIDs in $taxonomy file\n";
		exit;
	}
	$taxon = $tax{$taxon}[2];
}
else
{
	if(exists $name_taxids{$taxon}) 
	{
		foreach my $taxid (keys %{$name_taxids{$taxon}})
		{
			push(@taxids, $taxid);
		}
	}
}

if((scalar @taxids) == 0)
{
	print "$taxon is not in the list of taxa in $taxonomy file\n";
	exit;
}

my $out = $outdir.$taxon.'_subtaxa.tsv'; 
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "taxon	taxID	sub_taxon	taxID\n";
foreach my $taxid (@taxids)
{
	++$stat{'1. Number of taxIDs: '};
	my $taxlevel = int($tax{$taxid}[3] + 1); # get the taxlevel index for the next major taxonomic level
	foreach my $tid (keys %tax)
	{
		if($tax{$tid}[3] == $taxlevel) # correct taxlevel
		{
			my @lin = get_lineage_list($tid, \%tax, 0);
			for(my $i=0; $i < scalar @lin; ++$i)
			{
				if($lin[$i] == $taxid)
				{
					++$stat{'2. Number of subtaxa: '};
					print OUT "$taxon	$taxid	$tax{$tid}[2]	$tid\n";
					last;
				}
			}
		}
	}
}
close OUT;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

print LOG print_stat(\%stat, $t0);
close LOG;

exit;

#############################################

sub get_taxon
{
	my ($line, $name_taxids, $mismatch) = @_;
	
	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split("\t", $line);
	if( scalar @line > 1)
	{
		if(exists $$name_taxids{$line[0]}{$line[1]}) # name matches taxids
		{
			return ($line[0], $line[1]);
		}
		else
		{
			$$mismatch{$line[0]} = $line[1];
			return ($line[0], 0);
		}
	}
	else
	{
		return ($line[0], 0);
	}

}

#############################################
sub get_ranked_lin_from_tax
{
	my ($lin, $tax, $taxlevel_index, $name_index) = @_;
	# $taxlvel_index: index of the element in the list in @{$tax{taxid}}
	#$tax{$taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	# @lin = (list of taxids stating from the smallest)
	
	my @ranked_lin = ('', '', '', '', '', '', '', '');
	
	foreach my $taxid (@$lin)
	{
		my $taxlevel = $$tax{$taxid}[$taxlevel_index];
		if($taxlevel == int($taxlevel) and $taxlevel >= 1) # taxlevel index is an iteger
		{
			$ranked_lin[$taxlevel-1] = $$tax{$taxid}[$name_index];
		}
	}
	
	return @ranked_lin;
	
}
############################################

sub check_lineage_bis
{
	my ($tax, $taxid, $taxids_hash, $anti_taxids_hash) = @_; 
	
	# if taxid is a child of an existing taxid, add to %taxids, add to anti_taxids oterwise
	
	if(exists $$taxids_hash{$taxid})
	{
		return 1;
	}
	elsif(exists $$anti_taxids_hash{$taxid})
	{
		return 0;
	}
	else
	{
		my @lineage = get_lineage_list($taxid, $tax, 0);
		my $bool = 0;
		for(my $i = 0; $i < scalar @lineage; ++$i)
		{
			if(exists $$taxids_hash{$lineage[$i]})
			{
				$$taxids_hash{$taxid} = $$taxids_hash{$lineage[$i]};
				$bool = 1;
				last;
			}
		}

		unless($bool) # none of the taxid in the lineage is in tha %taxids_hash
		{
			$$anti_taxids_hash{$taxid} = '';
		}
		return $bool;
	}
}

#############################################

sub check_lineage
{
	my ($lineage, $taxid, $taxids_hash, $anti_taxids_hash) = @_; 
	
	# if taxid is a child of an existing taxid, add to %taxids, add to anti_taxids oterwise
	
	my $bool = 0;
	for(my $i = 0; $i < scalar @$lineage; ++$i)
	{
		if(exists $$taxids_hash{$$lineage[$i]})
		{
			$$taxids_hash{$taxid} = $$taxids_hash{$$lineage[$i]};
			$bool = 1;
			last;
		}
	}
	
	unless($bool) # none of the taxid in the lineage is in tha %taxids_hash
	{
		$$anti_taxids_hash{$taxid} = '';
	}
	return $bool;
}

#############################################
sub complete_names
{
	my ($ncbi_tax_names, $name_taxids) = @_; # add synonymes to the hash
	#$name_taxids{$name}{taxid} =''

	open(IN, $ncbi_tax_names) or die "Cannot open $ncbi_tax_names\n";
	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split("\t", $line);
		$$name_taxids{$line[1]}{$line[0]} = '';
	}
	close IN;

}

#############################################
sub read_taxonomy_local
{
	my ($taxonomy, $tax, $name_taxids_par) = @_;
	
	#$tax{$taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	#$name_taxids{$name}{taxid} = parent taxid


	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		@{$$tax{$line[0]}} = ($line[1], $line[2], $line[3], $line[5]);
		$$name_taxids_par{$line[3]}{$line[0]} = '';
	}
	close IN;

}

#############################################
sub print_help
{

print '
usage: perl select_taxa.pl [-options] -tsv INPUT_FILE -taxonomy TAXONOMY_TSV -taxon_list TAXON_LIST -outdir OUTDIR -out OUTFILE

 ARGUMENTS
   -taxon                  taxon name or a taxID
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
   -outdir                 Name of the otput directory
',  "\n";

  exit;
}
