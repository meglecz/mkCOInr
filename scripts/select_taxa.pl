use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
# -tsv file (seqID	TaxID	sequence)
# -taxonomy file (tax_id	parent_tax_id	rank	name_txt	old_tax_id	synonyms)
# -taxon_list : one taxon per line or tsv with taxon name in the first column and taxids in the second(optional)
# -min_taxlevel : species/genus/family/order/class/phylum (optional)

# AIM
# select sequences that belong to the taxa in a taxon list, and their rank >= of min_taxlevel
# If taxids are not given in the taxon_list file the script uses all taxids that matches the taxon name
# if taxon_list a lineage file is written for each taxa and the corresponding taxids 
#	-it should be checked manually, if lineages are coherent with the targat taxa
#	-homonymy column indicate if more than one taxid for a taxon
#	-if there are incoherent lineages make a new taxon_list file based on the lineage file including taxon names and taxids and rerun the script with the new taxon_list

# OUTPUT
# -tsv file (seqID	TaxID	sequence)
# lineage and taxid for all input taxa => check manually. if incoherent => make exclude list

###############################
###############################
my %params = (
'tsv' =>  '', # seqID	TaxID	sequence
'taxonomy' =>  '', 
'taxon_list' => '', # optional, taxon name in first column, taxid is available in second
'outdir' =>  '',
'out' => '',
'negative_list' => 0, # if 1 keep all taxa except the ones on the taxon list 
'min_taxlevel' => 'root', #species/genus/family/order/class/phylum/kingdom/superkingdom/root
);

modify_params_from_tags(\%params, \@ARGV);

my $tsv = $params{tsv};
my $taxon_list = $params{taxon_list};
my $negative_list = $params{negative_list};
my $min_taxlevel = $params{min_taxlevel};
my $taxonomy = $params{taxonomy};
my $outdir = $params{outdir};
my $out = $params{out};

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
my $log = $outdir.'select_taxa.log';
my %stat;
$out = $outdir.$out; 
my $out_lin = $outdir.'taxa_with_lineages.tsv';
my %taxlevel = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1, '', 0, 'root', 0);
my $min_taxlevel_index = $taxlevel{$min_taxlevel};

open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";
###

### Read taxonomy
print "\n####\nRead taxonomy\n";
print LOG "\n####\nRead taxonomy\n";

my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
my %merged; #$merged{merged_taxid} = up to date taxid
read_taxonomy_to_tax_hash_include_merged(\%tax, $taxonomy, \%merged);

my %name_taxids; # $name_taxids{name}{$taxid} = '' # All names, including synonyms and homonyms
my %taxid_names; # $taxid_names{taxid}{names} = '' # All names, including synonyms and homonyms
my %name_taxids_par; #$name_taxids_par{$name}{taxid} = parent taxid # only scientifi names
read_taxonomy_names(\%taxid_names, \%name_taxids, \%name_taxids_par, $taxonomy);
%taxid_names = (); # delete hash. It is produced automatically by the routine, but we do not need it.
%name_taxids_par = (); # delete hash. It is produced automatically by the routine, but we do not need it.
my %missing_taxids; # keep list of taxid not in the taxonomy file
print LOG "Runtime: ", time - $t, "s \n";
$t = time;




open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "seqID	taxID	sequence\n";
if($taxon_list) # if select for taxon list
{
	### Read input taxon_list and fill hashes
	my %taxids; #$taxids{Taxids of taxa in taxon list + taxids of descendant taxa} = taxID of the original taxon (the one in the input file); #do not check for homonymy
	my %anti_taxids; #  %anti_taxids{taxid} # taxid that are not in the list and not children of the list
	my %taxon_list; # $taxon_list{taxon name}{Taxids of taxa in taxon list + taxids of descendant taxa} = lineage
	my %seq_count;   # $seq_count{taxid of original taxon} = number of sequences
	my %mismatch; #$mismatch{taxname} =taxid in input file, the input taxon name and taxid do not match

	print "\n####\nRead input taxon_list\n";
	print LOG "\n####\nRead input taxon_list\n";
	open(IN, $taxon_list) or die "Cannot open $taxon_list\n";
	my $title = <IN>;
	open(LIN, '>', $out_lin) or die "Cannot open $out_lin\n";
	print LIN "taxon	taxID	homonymy	number of sequences	superkingdom	kingdom	phylum	class	order	family	genus	species\n";
	
	### get all taxids that correspond to the taxon names in the input file
	while(my $line = <IN>) 
	{
		# get name and taxid, if no taxid given in the input, returns 0, if taxid does not match to name, returns 0 and add name to %mismatch
		my ($taxon, $taxid) = get_taxon($line, \%name_taxids, \%mismatch);
		if($taxon) # avoid empty lines
		{
			
			if($taxid) # taxid is given for the taxon
			{
				my @lin = get_lineage_list($taxid, \%tax, 0);
				@{$taxon_list{$taxon}{$taxid}} = @lin;
				$taxids{$taxid} = $taxid;
				$seq_count{$taxid} = 0;
			}
			else # no taxid for the taxon in the input file => get taxids
			{
				if(exists $name_taxids{$taxon})
				{
					foreach my $taxid (keys %{$name_taxids{$taxon}})
					{
						unless(exists $taxids{$taxid}) # add taxid to %taxids with the list of taxids of its lineage as a value
						{
							my @lin = get_lineage_list($taxid, \%tax, 0);
							@{$taxon_list{$taxon}{$taxid}} = @lin;
							$taxids{$taxid} = $taxid;
							$seq_count{$taxid} = 0;
						}
					}
				}
				else
				{
					print LIN "$taxon	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA\n";
					++$stat{'1.3 Number of input taxa without taxids: '};
				}
			}
		}
	}
	close IN;
	$stat{'1.1 Number of input taxon names: '} = scalar keys %taxon_list;
	$stat{'1.2 Number of taxids for input taxon names: '} = scalar keys %taxids;
	
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	###
	
	### Select sequences from tsv file
	print "\n####\nSelect sequences from tsv file\n";
	print LOG "\n####\nSelect sequences from tsv file\n";
	open(TSV, $tsv) or die "Cannot open $tsv\n";
	$title = <TSV>;
	while(my $line = <TSV>)
	{
		++$stat{'2.1 Number of sequences in the input tsv: '};
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		my $taxid = $line[1];
		if(exists $merged{$taxid}) # update taxid if a merged taxid is used instead of an up to date
		{
			$taxid = $merged{$taxid};
		}
		if($tax{$taxid}[3] >= $min_taxlevel_index) ### taxlevel is at least $min_taxlevel
		{
			# if taxid is a child of an existing taxid, add to %taxids, add to %anti_taxids otherwise
			# returns 1 in taxin is %taxids 0 otherwise
			my $on_the_list = check_lineage_bis(\%tax, $taxid, \%taxids, \%anti_taxids); # if taxid is a child of an existing taxid, add to %taxids, add to %anti_taxids otherwise

			if ($on_the_list and $negative_list ==0) # taxid id on the taxlist or its a child  AND sequences on the list should be printed 
			{
				print OUT $line, "\n"; # print out line if child taxid
				++$stat{'2.2 Number of sequences in the output tsv: '};
			}
			elsif ($on_the_list == 0 and $negative_list) # taxid is NOT on the taxilist nor its a child  AND negative selection
			{
				print OUT $line, "\n"; # print out line if child taxid
				++$stat{'2.2 Number of sequences in the output tsv: '};
			}
			
			if($on_the_list)
			{
				++$seq_count{$taxids{$taxid}};
			}
		}
	}
	close TSV;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####

	
	### Write lineage file for all taxa in the taxon list
	print "\n####\nWrite lineage file for all taxa in the taxon list\n";
	print LOG "\n####\nWrite lineage file for all taxa in the taxon liste\n";
	foreach my $taxon (sort keys %taxon_list)
	{
		my $homonymy = 0;
		if(scalar keys %{$taxon_list{$taxon}} > 1)
		{
			$homonymy = 1;
			++$stat{'1.4 Number of input taxa with more than one taxid: '};
		}
		foreach my $taxid (keys %{$taxon_list{$taxon}})
		{
			my @ranked_lin = get_ranked_lin_from_tax($taxon_list{$taxon}{$taxid}, \%tax, 3, 2); #(\@lin, \%tax, $taxlevel_index, $name_index)
			print LIN "$taxon	$taxid	$homonymy	$seq_count{$taxid}	",join("\t", @ranked_lin),"\n";
		}
	}
	close LIN;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	###
	
	if(scalar keys %mismatch)
	{
		print "WARNING:\nThe following names did not match with the taxIDs given in the taxon_list file.\nWhen possible taxIDs were established from the taxonomy file:\n";
		print join("\n", sort keys %mismatch), "\n";
		print LOG "WARNING:\nThe following names did not match with the taxIDs given in the taxon_list file.\nWhen possible taxIDs were established from the taxonomy file:\n";
		print LOG join("\n", sort keys %mismatch), "\n";
	}
}
else # select for taxlevel, but not for taxon list
{
	### Select sequences annotated to the min_taxlevel
	print "\n####\nSelect sequences annotated to the min_taxlevel\n";
	print LOG "\n####\nSelect sequences annotated to the min_taxlevel\n";
	open(TSV, $tsv) or die "Cannot open $tsv\n";
	my $title = <TSV>;
	while(my $line = <TSV>)
	{
		++$stat{'1.1 Number of sequences in the output tsv: '};
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		my $taxid = $line[1];
		
		if(exists $merged{$taxid}) # old taxid is replaced by up to date one
		{
			$taxid = $merged{$taxid};
		}
		
		if(exists $tax{$taxid})
		{
			if($tax{$taxid}[3] >= $min_taxlevel_index) ### taxlevel is at least $min_taxlevel
			{
				++$stat{'1.2 Number of sequences in the output tsv: '};
				print OUT $line, "\n";
			}
		}
		else
		{
			print "WARNING: $line[1] taxID is not present in taxonomy file\n";
			$missing_taxids{$line[1]} = '';
		}
	}
	close TSV;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;

}
close OUT;
print Dumper(\%missing_taxids);

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
   -taxon_list             Txt file with one taxon per line, OR  a tsv file with taxon_name and taxid in the first two columns
                              Not necessary if sequences are selectes only by their taxonomic resolution
   -tsv                    Input tsv file with seqID,taxID,sequence
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
   -outdir                 Name of the otput directory
   -out                    Name of the output file
 OPTIONS
   -negative_list          [1/0]; Default: 0
                              If 1, keeps all taxa except the ones on the taxon list
   -min_taxlevel           [species/genus/family/order/class/phylum/kingdom/root]; Default: root

',  "\n";
  exit;
}
