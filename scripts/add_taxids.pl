use warnings;
use strict;
use mkdb;
use Data::Dumper;


# INPUT
# lineages (BOLD, custom or pooled) with the following columns. Identical lineages can appear in multiple lines (e.g. if poodd from format_bold and format_custom output)
	# phylum_name	class_name	order_name	family_name	subfamily_name	genus_name	species_name	seqIDs
# sequences tsv file with the following columns (without heading)
	#sequenceID	sequence
# outdir
# taxonomy already existing taxonomy file, can contain negative taxids, so these taxid will be taken into account while giving arbitrary taxids to new taxa 
#	Columns: tax_id	parent_tax_id	rank	name_txt	old_tax_id	rank index	synonyms

# ALGO
# Finds the smallest taxon that matches an alredy existing taxID
	# accept taxID if at least 0.6 the taxa in the upward lineages matches the lineage of the taxID
	# (for species level do not count genus, since it matches necessarily the species name)
	# OR 
	# it is a phylum in BOLD that corresponds to a phylum in taxonomy.tsv
	# if taxID with <= 0.25 matches => go to the next taxlevel
	# if taxID with 0.25 < matches < 0.6 => print to ambiguous file, and it should be checked manually
# 	This search for ncbi taxid, takes into account 
#		- all existing taxiID (ncbi or previous arbitrary) 
#		- existing synonyms
#		- if more than one taxid for name
#				- takes the one with highest proportion of taxa matching the upwards lineage
#				- then the same tax level

# For each lineages that do not have yet taxid, assign arbitrary negative taxIDs to taxa and link it as a child to an existing taxID 

# OUTPUT

#sequences_with_taxIDs.tsv:
#		#seqID	TaxID	sequence

# lineages_with_taxIDs.tsv: One line per lineage
# 	lowest_taxname	lowest_rank	lowest_TaxID	match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbiTaxID	phylum	class	order	family	subfamily	genus	species	seqIDs
# 	lowest_taxname	lowest_rank	lowest_TaxID	match_phylum	class	order	family	subfamily	genus	species	seqIDs
#		(Different lines can have the same taxid if synonyms in ncbi)

# ambiguous_lineages.tsv; this file should be examined and modified manually by experts 
#		match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbi_TaxID	phylum	class	order	family	subfamily	genus	species	ncbi_superkingdom	ncbi_kingdom	ncbi_phylum	ncbi_class	ncbi_order	ncbi_family	ncbi_genus	ncbi_species	ncbi_taxname	seqIDs
# ambiguous_sequneces.tsv => can be used to rerun add_taxids.pl after modifing the ambiguous_lineages.tsv

# taxonomy_updated.tsv updated with the new arbitrary taxIDs

###############################
###############################
my %params = (
'lineages' =>  '', # tsv file with cleaned BOLD and/or custom lineages
'sequences' => '', # optinal. It is just a check if all seqIDs are unique
'outdir' =>  '',
'taxonomy' =>  '' # already existing taxonomy file, can contain negative taxids, tax_id	parent_tax_id	rank	name_txt	old_tax_id	rank index	synonyms
						# if empty, it will be produced from ncbi tax dmp files
);
modify_params_from_tags(\%params, \@ARGV);

my $lineages = $params{lineages};
my $sequences = $params{sequences};
my $outdir = $params{outdir};
my $taxonomy = $params{taxonomy};

$outdir = add_slash_to_dir($outdir);

my %bold_taxlevel = ('species' => 10,'genus' => 9,'subfamily' => 8,'family' => 7,'order' => 6,'class' => 5,'phylum' => 4);
my $sep = "\t";
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
my $log = $outdir.'add_taxids.log';
# all lineages with non-ambiguous ncbi taxids with detais
my $out_details = $outdir.'lineage_comparaison.tsv'; 
# all sequences with ambiguous taxids with details
my $out_boubt = $outdir.'ambiguous_lineages.tsv';
my $out_seq_bp = $outdir.'ambiguous_sequences.tsv';
my $out = $outdir.'lineages_with_taxIDs.tsv'; 
my $out_seq = $outdir.'sequences_with_taxIDs.tsv';  #seqID	TaxID	sequence
my $taxonomy_out = $outdir.'taxonomy_updated.tsv'; # tax_id	parent_tax_id	rank	name_txt	old_tax_id



open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";

my %ncbi_ranked_tax_index = ('ncbi_taxname' => 0, 'species' => 1,'genus' => 2,'family' => 3,'order' => 4,'class' => 5,'phylum' => 6,'kingdom' => 7,'superkingdom' => 8);
my @ncbi_ranked_tax = ('ncbi_taxname', 'ncbi_species','ncbi_genus','ncbi_family','ncbi_order','ncbi_class','ncbi_phylum','ncbi_kingdom','ncbi_superkingdom');
my %stat;
####

my %seq;
####
#### Read sequences to check if all seqids are unique. 
# This is just a check, sequences are needed only to print out ambiguous lineages with their sequences
print "\n####\nReading sequences\n";
print LOG "\n####\nReading sequences\n";

my %non_uniq;
read_seq_file(\%seq, \%non_uniq, $sequences);
$stat{'0. Number of input sequneces : '} = scalar (keys %seq);
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

#### Stop if not all seqIDs are uniq
if(scalar keys %non_uniq)
{
	print "\nERROR\nThe following sequence IDs are not uniq. Rename these sequences:\n";
	print join("\n", sort keys %non_uniq), "\n";
	print LOG "\nERROR\nThe following sequence IDs are not uniq. Rename these sequences:\n";
	print LOG join("\n", sort keys %non_uniq), "\n";
	exit;
}
####
####


####
#### Read taxonomy
print "\n####\nRead taxonomy\n";
print LOG "\n####\nRead taxonomy\n";

my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
read_taxonomy_to_tax_hash(\%tax, $taxonomy);

my @taxids = sort {$a <=> $b} keys %tax;
my $smallest_taxid_hash = $taxids[0]; # used for finding new taxids
my $smallest_taxid_input = $smallest_taxid_hash;

my %name_taxids; # $name_taxids{name}{$taxid} = '' # All names, including synonyms and homonyms
my %taxid_names; # $taxid_names{taxid}{names} = '' # All names, including synonyms and homonyms
my %name_taxids_par; #$name_taxids_par{$name}{taxid} = parent taxid # only scientifi names
read_taxonomy_names(\%taxid_names, \%name_taxids, \%name_taxids_par, $taxonomy);

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

$stat{'1.0 Number of taxids in taxonomy file: '} = scalar keys %tax;
$stat{'1.1 Number of scientific names in taxonomy file: '} = scalar keys %name_taxids_par;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####
####




#### Read lineages and pool identical lineages
print "\n####\nReading lineages\n";
print LOG "\n####\nReading lineages\n";
my %lin; # $lin{lineage} = (list of seqids) # exemple lineage 'Arthropoda	Insecta	Coleoptera	Scarabaeidae	Dynastinae	Oryctes	Oryctes borbonicus	NA'
#@bold_ranks taxlevel list deduced from title line
my @bold_ranks = read_lineage_file(\%lin, $lineages); # pool identical lineages on the same line and their corresponcing seqIDs as well

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####


#### open output files
open(OD, '>', $out_details) or die "Cannot open $out_details\n";
print OD "match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbi_TaxID	", join("\t", reverse @bold_ranks), "\t", join("\t", reverse @ncbi_ranked_tax), "\tseqIDs\n";
open(DOUBT, '>', $out_boubt) or die "Cannot open $out_boubt\n";
print DOUBT "match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbi_TaxID	", join("\t", reverse @bold_ranks), "\t", join("\t", reverse @ncbi_ranked_tax), "\tseqIDs\n";
open(DOUBTSEQ, '>', $out_seq_bp) or die "Cannot open $out_seq_bp\n";
print DOUBTSEQ "seqID	sequence\n";
####

##### Match to existing taxid
	#### Searching for ncbi taxIDs
	print "\n####\nSearching for ncbi taxIDs\n\n";
	print LOG "\n####\nSearching for ncbi taxIDs\n";
	# go thougth all lineages one by one
	foreach my $line (sort keys %lin) #$lin{lineage} = (list of seqids)
	{
		++$stat{'2.0 Number of input lineages to be assigned : '};
		my @line = split("\t", $line);
		my $tmp = pop@line; # eliminate the last dummy element that insured that there will be the correct number of elements in the list, even if that last ones are empty

		# make @lineage_bold from BOLD fields starting from the lowest rank (species	genus	subfamily	family	order	class	phylum)
		my @lineage_bold = reverse @line; # from species to phylum

		my $taxid = 0;
		my $match = 0;
		# start from smallest taxon and check is it has a ncbi tax => classed as ok or doubious (print out) go to next level otherwise
		# at this step assign a ncbi taxid (or an already existing negative taxid) to each sequence at the smallest level possible, ignore lower(higher resuolution) levels
		for(my $i = 0; $i < scalar @lineage_bold; ++$i)# get ncbi taxid for the smallest taxon present in ncbi tax
		{
			# 0 if no ncbi taxid at all
			# if more than one taxid => choose the one with highest portion of the BOLD lineage (above $i) that matches ncbi lineage => then the one that has a the same rank
			# $match: the proportion of taxnames matching between ncbi lineage and @lineage_bold (above $i)
			($taxid, $match) = get_taxid_new($i, \@lineage_bold, \%name_taxids, \%taxid_full_lineage_id, \@bold_ranks, \%tax, \%taxid_names);
			
			if($taxid) # valid taxid
			{
				if($match >= 0.6 or ($i == (scalar @lineage_bold-1) and $tax{$taxid}[1] eq $bold_ranks[$i])) # good match or highest taxlevel with matching rank => print and quit the loop
				{
					print OD  $match, "\t",  $tax_ranked_lineage{$taxid}[0], "\t", $bold_ranks[$i], "\t", $taxid, "\t", join("\t", reverse @lineage_bold), "\t", join("\t", reverse @{$tax_ranked_lineage{$taxid}}), "\t", join(";", @{$lin{$line}}), "\n";
					++$stat{'2.1 Number of lineages with unambiguous ncbi taxid: '};
					last; # stop if there if convincing taxid
				}
				elsif($match <= 0.25) # lineage is soo different that it is very likely to be a homonym to a incorrect taxa => go to the next taxlevel
				{
					$taxid = 0;
					$match = 0;
					next;
				}
				else # Print out if the match id between 0.25 and 0.6 not convincing enough but not good enough either.
						# This needs human checking to avoid going further up in the lineage and assign a new taxid lower taxa that has already
				{
					print DOUBT $match,  "\t",  $tax_ranked_lineage{$taxid}[0], "\t", $bold_ranks[$i], "\t", $taxid, "\t", join("\t", reverse @lineage_bold), "\t", join("\t", reverse @{$tax_ranked_lineage{$taxid}}), "\t",  join(';', @{$lin{$line}}), "\n";
					++$stat{'2.2 Number of lineages with ambiguous ncbi taxid: '}; 
					foreach my $sid ( @{$lin{$line}}) # print out sequnces with no precise taxid
					{
						print DOUBTSEQ "$sid	$seq{$sid}\n";
					}
					last; # stop if boubtful taxid, to avoid adding non-valid taxon name latter to the taxid system 
				}
			}
		}
		unless($taxid) # BOLD lineage cannot be matched to an NCBI lineage
		{
			++$stat{'2.3 Number of lineages with no ncbi taxid: '};
			print DOUBT "\t" x 4, join("\t", reverse @lineage_bold),  "\t" x 10, join(';', @{$lin{$line}}), "\n"; # print lineage and seq ids
			foreach my $sid ( @{$lin{$line}}) # print out sequnces with no precise taxid
			{
				print DOUBTSEQ "$sid	$seq{$sid}\n";
			}
		}
	}
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	close IN;
	close OD;
	close DOUBTSEQ;
	####
#### end Match to existing taxid


my %tax_par;
my %tax_rank;
my %tax_name;
#### add_new taxID
	#### Get non-ncbi taxids
	print "\n####\nAssigning nonNCBI taxiIDs\n";
	print LOG "\n####\nAssigning nonNCBI taxiIDs\n";

	my $taxid_lineages = $out_details; 
	open(IN, $taxid_lineages) or die "Cannot open $taxid_lineages\n";
	my $title = <IN>;
	$title =~ s/\s*$//;
	my @title = split($sep, $title);

	open(OUT, '>', $out) or die "Cannot open $out\n";
	print OUT "lowest_taxname	lowest_rank	lowest_TaxID	phylum	class	order	family	subfamily	genus	species	seqIDs\n";
	open(SEQ, '>', $out_seq) or die "Cannot open $out_seq\n";#seqID	TaxID	sequence
	print SEQ "seqID	taxID	sequence\n";


	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split($sep, $line);
		#match_lineage_proportion	ncbi_taxname	ncbi_taxlevel	ncbi_TaxID	phylum	class	order	family	subfamily	genus	species	ncbi_superkingdom	ncbi_kingdom	ncbi_phylum	ncbi_class	ncbi_order	ncbi_family	ncbi_genus	ncbi_species	ncbi_taxname	seqIDs

		# my %bold_taxlevel = ('species' => 10,'genus' => 9,'subfamily' => 8,'family' => 7,'order' => 6,'class' => 5,'phylum' => 4);
		my $taxid_parent = $line[3]; # use it for lower levels
		my $ncbi_taxrank = $line[2];
		my $smallest_taxon_index = $bold_taxlevel{$ncbi_taxrank}; # the smallest taxlevel with ncbi taxid
		my $taxrank =  $title[$smallest_taxon_index];
		my $taxid = $line[3];

		for(my $i = ($smallest_taxon_index +1); $i <= $bold_taxlevel{species} ; ++$i) # go through taxon names starting from highest taxlevel without ncbi taxid
		{
			if($line[$i]) # if field is not empty, get a negative taxid for a taxon name 
							#check if the name with the same parent is already in tax hashes => if yes, take the taxid, otherwise establish a new)
			{
				$taxrank =  $title[$i];
				# check if name exists already with the same parent taxid. If not create a new taxid, and complete taxonomy hashes
				# $smallest_taxid_hash the smallest taxids in the %tax hash
				($taxid, $smallest_taxid_hash) = get_new_taxid_simple($line[$i], $taxid_parent, $taxrank, \%name_taxids_par, \%tax, $smallest_taxid_hash);
				
				$taxid_parent = $taxid; # reintialize for the next level
				$smallest_taxon_index = $i;
			}
		}
		my $seq_id_list = pop @line;
		my $lowest_taxname = $line[$smallest_taxon_index];
		@line = splice(@line, 4, 7);
		print OUT "$lowest_taxname	$taxrank	$taxid	", join("\t", @line), "\t$seq_id_list\n";
		my @seqids = split(';', $seq_id_list);
		foreach my $sid (@seqids)
		{
			print SEQ "$sid	$taxid	$seq{$sid}\n";
			++$stat{'3.3 Number of sequences with taxID : '};
		}
	}
	close IN;
	close OUT;
	close SEQ;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####

	####
	print "\n####\nWriting updated taxonomy file\n";
	print LOG "\n####\nWriting updated taxonomy file\n";
	$stat{'3.0 Number of taxids in output taxonomy file: '} = scalar keys %tax;
	$stat{'3.1 Number of taxnames in output taxonomy file: '} = scalar keys %name_taxids_par;
	$stat{'3.2 Number of new taxids added to the taxonomy file: '} = $smallest_taxid_input - $smallest_taxid_hash;
	update_taxonomy($taxonomy, $taxonomy_out, \%tax);
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####
#### $add_non_ncbi_taxid option

system 'rm '.$out_details;

print LOG print_stat(\%stat, $t0);
close LOG;
exit;


####################################################
####################################################


################################################
sub read_lineage_file
{
	my ($hash, $file) = @_;
	
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>; # phylum	class	order	family	subfamily	genus	species	seqIDs
	
	$title =~ s/\s*$//; 
	$title =~ s/"//g;
	my @bold_ranks = split("\t", $title);
	@bold_ranks = reverse @bold_ranks;
	my $temp = shift(@bold_ranks); # species	genus	subfamily	family	order	class	phylum
	
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		my $seqids = pop @line;
		my $l = join("\t", @line);
		$l .= "\tNA"; # add ; so if the last element of the list is empty, it will be still taken into account
		my @seqids = split(';', $seqids);
		foreach my $sid (@seqids)
		{
			push(@{$$hash{$l}}, $sid);
		}
	}
	close IN;
	return @bold_ranks;
}


################################################

sub read_seq_file
{
	my ($hash, $non_uniq, $file) = @_;
	
	open(IN, $file) or die "Cannot open $file\n";
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		if(exists $$hash{$line[0]})
		{
			$$non_uniq{$line[0]} = '';
		}
		$$hash{$line[0]} = $line[1];
	}
	close IN;
}


#############################################

sub update_taxonomy
{
	my ($taxonomy, $taxonomy_out, $tax) = @_;

## in taxonomy file, taxid can be duplicated if there are more than one merged taxid
# e.g.
#2834507	139133	species	Sistotremastrum vigilans	2698470
#2834507	139133	species	Sistotremastrum vigilans	2704103
# only the last line is kept in the %tax hash, which is not a problem for defining new taxids, 
# but it is a problem when printing out updated taxonomy file => print taxonomy to taxonomy_out and add new taxids afterwards

	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	open(OUT, '>', $taxonomy_out) or die "Cannot open $taxonomy_out\n";
#	print OUT "tax_id	parent_tax_id	rank	name_txt	old_tax_id	tax_rank_index	synonyms\n";
	
	while(my $line = <IN>) # print out original lines, 
	{
		print OUT $line;
		my @line = split("\t", $line);
		if(exists $$tax{$line[0]}) 
		{
			delete $$tax{$line[0]}; # delete taxid present in input file
		}
	}
	close IN;
	
	foreach my $taxid (reverse sort {$a <=> $b} keys %$tax ) # print out new taxids
	{
		print OUT "$taxid	$$tax{$taxid}[0]	$$tax{$taxid}[1]	$$tax{$taxid}[2]		$$tax{$taxid}[3]	\n"; # tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel	synonyms
	}
	close OUT;
}
#############################################
sub get_new_taxid_simple
{
	my ($name, $parent_taxid, $rank, $name_taxids_par, $tax, $smallest_taxid) = @_;
#	$tax{$taxid} = (parent_tax_id	rank	name_txt	taxrank index)
# tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxrank index
#	$name_taxids_par{$name}{taxid} = parent taxid
	
	my %taxlevel_ind = ('species',8,'genus',7,'family',6,'subfamily',6.5,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1); # numerical score for each major tax level

	
	my $taxid = 0;
	if(exists $$name_taxids_par{$name}) # name is already in tax hashes
	{
		foreach my $tid (keys %{$$name_taxids_par{$name}}) # parents match
		{
			if($$name_taxids_par{$name}{$tid} == $parent_taxid)
			{
				$taxid = $tid;
				last;
			}
		}
	}
	
	unless($taxid) # new taxid should be established => one less then the smallest existing taxid
	{
		--$smallest_taxid;
		$taxid = $smallest_taxid;

		$$name_taxids_par{$name}{$taxid} = $parent_taxid;
		@{$$tax{$taxid}} = ($parent_taxid, $rank, $name, $taxlevel_ind{$rank});
	}
	
	return ($taxid, $smallest_taxid);
}


#############################################

sub get_taxid_new
{
	my ($i, $lin, $name_taxids, $taxid_full_lineage_id, $bold_ranks, $tax, $taxid_names) = @_;
#($i, \@lineage_bold, \%name_taxids, \%taxid_full_lineage_id, \@bold_ranks, \%tax, \%taxid_names);

	if(exists $$name_taxids{$$lin[$i]}) # ncbi taxid exists
	{
		my %taxid_temp; # $taxid_temp{$taxid} = ($match, $nb_tax)
		foreach my $taxid ( keys %{$$name_taxids{$$lin[$i]}}) # get taxid with the highest proportion of matches in the bold lineages above $i, or if it is <0.5, the taxon where ncbi and bold ranks are matching with at least 1,2 matching taxon name in the lineage
		{
			# return 
				# $match: the proportion of taxname matching between ncbi linegae and @lin (bold above $i) ($match_p)
				# $nb_tax: number of taxa in bold lineage above $i (ignore empty)
			# Include in ncbi lineage all synonymes
			@{$taxid_temp{$taxid}} = compare_lineages_new($i, $lin, $$taxid_full_lineage_id{$taxid}, $taxid_names);
		}
		my ($taxid, $match) = get_best_match_new(\%taxid_temp, $$bold_ranks[$i], $tax); # take the taxid with the highest $match_p > same taxlevel
		return ($taxid, $match);
	}
	else # no ncbi taxid
	{
		return (0, 0);
	}

}

####################################################

sub get_best_match_new
{
	my ($taxids, $bold_rank, $tax) = @_; # take the taxid with the highest $match_p => same taxlevel
 # $taxids{$taxid} = ($match_p, $nb_tax)
 #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
 
	my $best_taxid = 0;
	my $best_match = -1;
	my $best_nb_tax = -1;
	foreach my $taxid (keys %$taxids)
	{
		my $match = $$taxids{$taxid}[0];
		my $nb_tax = $$taxids{$taxid}[1];
		if($match > $best_match)
		{
			$best_taxid = $taxid;
			$best_nb_tax = $nb_tax;
			$best_match = $match;
		}
		elsif( $match == $best_match )
		{
			if($$tax{$taxid}[1] eq $bold_rank)
			{
				$best_taxid = $taxid;
				$best_nb_tax = $nb_tax;
				$best_match = $match;
			}
		}
	}
 
	return ($best_taxid, $best_match);

}


####################################################
sub compare_lineages_new
{
	my ($i, $lin, $ncbi_lin, $taxid_names) = @_;

	my $match = 0; #number of taxname matching between ncbi linegae and @lin (bold above $i)
	my $nb_tax = 0; # number of tax in bold lienage above $i (ignore empty)

	my %ncbi_hash; # $ncbi_hash{$taxname} = '';
	# make a hash with all ncbi taxon names in the lineage (including synonymes)
	foreach my $id (@$ncbi_lin) # go though ncbi lineage ids
	{
		foreach my $taxname (keys %{$taxid_names{$id}}) # take all synonymes of a taxon names in the lineage
		{
			$ncbi_hash{$taxname} = '';
		}
	} 

# count matching taxa above the actual taxon
	if($i == 0) # if species level, do not compare the genera, since it suppposed to be OK
	{
		++$i;
	}
	for(my $j = $i+1; $j < scalar @$lin; ++$j)
	{
		if($$lin[$j]) # taxon is not empty
		{
			++$nb_tax;
			if (exists $ncbi_hash{$$lin[$j]} )
			{
				++$match;
			}
		}
	}
	
	if($nb_tax)
	{
		my $p = $match/$nb_tax;
		return ($p, $nb_tax); #proportion of the taxon names in the lienage above the actual taxon, that matches the ncbi lineage
	}
	else
	{
		return (0,0);
	}
}

############################################

sub print_help
{

print '
usage: perl add_taxids.pl -lineages INPUT_LINEAGE_TSV -sequences INPUT_SEQ_TSV -outdir OUTDIR -taxonomy TAXONOMY_TSV

 ARGUMENTS
   -lineages               Input tsv file with lineages: phylum,class,order,family,subfamily,genus,species,seqIDs
   -sequences              Input tsv file with sequences: seqID,sequence
   -outdir                 Name of the otput directory
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
', "\n";
  exit;
}


