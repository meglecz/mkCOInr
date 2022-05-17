use warnings;

############################################

sub add_slash_to_dir
{
 my ($dir) = @_;

 unless($dir eq '')
 {
	 $dir =~ s/\\/\//g;
	 unless($dir =~ /\/$/)
	 {
			$dir .= '/';
	 }
 }
 return $dir;
}

#############################################

sub check_name
{
	my ($name, $rank) = @_;


#		/^[A-Z][a-z-]+\s[a-z-]+\s[a-z-]+$/ # subspecies level
#		/^[A-Z][a-z-]+\s[a-z-]+$/ # species level
#		/^[A-Z][a-z-]+$/ All other levels	

	my $bool = 0;
	if($rank eq 'subspecies')
	{
		if($name =~ /^[A-Z][a-z-]+\s[a-z-]+\s[a-z-]+$/)
		{
			 $bool = 1;
		}
	}
	elsif ($rank eq 'species')
	{
		if($name =~ /^[A-Z][a-z-]+\s[a-z-]+$/)
		{
			$bool = 1;
		}
	}
	else
	{
		if($name =~ /^[A-Z][a-z-]+$/)
		{
			$bool = 1;
		}
	}
	return $bool;
}


###################################

sub clean_seq
{
	my ($seq, $max_n) = @_;
	
	$seq = uc $seq; 
	$seq =~ s/-//g;
	$seq =~ s/[^ATCG]/N/g;
	$seq =~ s/^N+//;
	$seq =~ s/N+$//;
	if($seq =~ /N{$max_n,}/) # eliminate sequences with $max_n or more Ns
	{
		$seq = '';
	}
	return $seq;
}

#################################################
sub get_date
{
	my @date = localtime;

	my $y = $date[5] + 1900;
	my $m = $date[4] + 1;
	if((length $m) == 1)
	{
		$m = '0'.$m;
	}
	if((length $date[3]) == 1)
	{
		$date[3] = '0'.$date[3];
	}
	if((length $date[2]) == 1)
	{
		$date[2] = '0'.$date[2];
	}
	if((length $date[1]) == 1)
	{
		$date[1] = '0'.$date[1];
	}
	if((length $date[0]) == 1)
	{
		$date[0] = '0'.$date[0];
	}
	my $d = $y.'-'.$m.'-'.$date[3].'-'.$date[2].'-'.$date[1].'-'.$date[0];
	return $d;
}

#################################################
sub get_file_list_from_folder
{
 my ($folder, $file_motif) = @_;
 
  unless ( opendir(FOLDER, $folder) )
  {
      print "Cannot access to folder $folder\n";
      exit;
  }

my @filenames = grep ( !/^\.\.?$/, readdir(FOLDER) );

#print "@filenames\n";
closedir(FOLDER);
my @files = ();

foreach my $file (sort @filenames)
{
	if ($file =~ /$file_motif/)
	{
		push(@files, $file);
	}
}
@filenames = ();
#print "@files\n";
return @files;
}

###################################
sub get_lineage_list
{
	my ($taxid, $tax, $tax_par_ind) = @_;
	# $tax{taxid} = (tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel)
	#$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
	#$tax_par_ind: index of taxpar_id in the list
	
	my @lin;
	while($taxid != 1)
	{
		push(@lin, $taxid);
		my $tpar = $$tax{$taxid}[$tax_par_ind];
		unless($tpar)
		{
			print "ERROR: $taxid taxID do not have a parent in tax hash\n";
			last;
		}
		$taxid = $tpar;
	}
	return @lin;
}


##############################################
sub get_lowest_correct_latin_name
{
	my ($taxid, $tax_rank, $tax_name, $tax_parent) = @_;
	
#	my %tax_name; #%tax_name{taxid} = scientific name
#	my %tax_rank; # %tax_rank{taxid} = rank
#	my %tax_par; # %tax_par{taxid} = taxid parent

	my $name = $$tax_name{$taxid};
	my $rank = $$tax_rank{$taxid};
	
	while(check_name($name, $rank) == 0) # not a correct latin name format for the rank
	{
		$taxid = $$tax_parent{$taxid};
		$name = $$tax_name{$taxid};
		$rank = $$tax_rank{$taxid};
		if($taxid == 1)
		{
			last;
		}
	}
	return $taxid;
}

##############################################

sub get_taxlevel
{
	my ($taxid, $tax_ranked_lineage, $tax_rank) = @_;
	#my %tax_ranked_lineage =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
	my $taxlevel = 0;
	
	my %hash = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1); # numerical score for each majot tax level
	my %hash2 = (0=>9, 1=>8, 2 => 7,  3=> 6, 4=> 5, 5=> 4,  6=>3 , 7=> 2, => 8, => 1,  9=> 0); # transform position in @{$ranked_lin{taxid}} => to numerical score

	if(exists $hash{$$tax_rank{$taxid}}) # get numerical score for major taxlevels
	{
		$taxlevel =  $hash{$$tax_rank{$taxid}};
	}
	else # taxlevel between major taxlevels => get an itermediate score (e.g. 7.5 for taxa between genus and species, 8.5 for bellow species)
	{
		for(my $i = 1; $i < 9; ++$i)
		{
			if($$tax_ranked_lineage{$taxid}[$i])
			{
				$taxlevel = $hash2{$i} + 0.5;
				last;
			}
		}
	}
	return $taxlevel;
}

##############################################
sub make_ranked_lineage
{
	my ($lineage_list, $tax) = @_;
	#tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	
	my @lin = ('root', '', '','','','','','',''); # root + 8 taxlevel
	foreach my $taxid (@$lineage_list)
	{
		my $tli = $$tax{$taxid}[3];
		if($tli == int($tli)) # major tax level
		{
			$lin[$tli] = $$tax{$taxid}[2];
		}
	}
	return @lin;
}

##############################################

sub make_taxonomy_with_rank_levels_synonyms
{
	my ($ncbitax_dir, $outdir, $motif) = @_;
	
	# output is a taxonomy file with the following columns:
	# tax_id	parent_tax_id	rank	name_txt	old_tax_id(taxIDs merged to another)	taxlevel	list of all names corresponding to the taxid
	# one extra line per merged taxids (old_tax_id)
	
	my $ncbi_tax_rankedlin = $ncbitax_dir.'rankedlineage.dmp';
	my $ncbi_tax_nodes = $ncbitax_dir.'nodes.dmp';
	my $merged_dump = $ncbitax_dir.'merged.dmp';
	my $names_dump = $ncbitax_dir.'names.dmp';
	
	my %tax_names; # $tax_names{taxid}{all names} = '';
	my %tax_par; # $tax_par{taxid} = taxid parent
	my %tax_rank; # $tax_rank{taxid} = rank
	my %tax_ranked_lineage; # $ranked_lineage{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
							# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	my %merged; # $merged{taxid} = (list of old taxids)

	# $tax_names{taxid}{all names} = '';
	read_new_names_dmp_to_hash_simple($names_dump, \%tax_names);
	# $tax_ranked_lineage{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
	# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	# $tax_ranked_lineage{taxid}[0] scientific name of the taxon
	read_ncbitax_rankedlin($ncbi_tax_rankedlin, \%tax_ranked_lineage);
	# $tax_par{taxid} = taxid parent
	# $tax_rank{taxid} = rank
	read_new_nodes_dmp_to_hash($ncbi_tax_nodes, \%tax_par, \%tax_rank);

	read_new_merged_dmp_to_hash($merged_dump, \%merged);

	my $taxonomy = $outdir.'taxonomy.tsv';
	if($motif)# add motif to filename
	{
		$taxonomy =~ s/\.tsv/$motif.tsv/;
	}

	open(OUT, '>', $taxonomy) or die "Cannot open $taxonomy\n";
	print OUT "tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel	synonyms\n";
	
	foreach my $taxid (sort {$a <=> $b} keys %tax_rank)
	{
		my $taxlevel = get_taxlevel($taxid, \%tax_ranked_lineage, \%tax_rank);
		my $sci_name = $tax_ranked_lineage{$taxid}[0];
		
		delete $tax_names{$taxid}{$sci_name};  # eliminate the scientific name from the synonyms

		if(exists $merged{$taxid})
		{
			foreach my $old_taxid (@{$merged{$taxid}})
			{
				print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$sci_name	$old_taxid	$taxlevel	",join(';', keys %{$tax_names{$taxid}}),"\n";
			}
		}
		else
		{
			print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$sci_name		$taxlevel	",join(';', keys %{$tax_names{$taxid}}),"\n";
		}
	}

	close OUT;
	
	return $taxonomy
}
##############################################

sub make_taxonomy_with_rank_levels
{
	my ($ncbitax_dir, $outdir, $motif) = @_;
	
	# output is a taxonomy file with the following columns:
	# tax_id	parent_tax_id	rank	name_txt	old_tax_id(taxIDs merged to another)	taxlevel
	# one extra line per merged taxids (old_tax_id)
	
	my $ncbi_tax_rankedlin = $ncbitax_dir.'rankedlineage.dmp';
	my $ncbi_tax_nodes = $ncbitax_dir.'nodes.dmp';
	my $merged_dump = $ncbitax_dir.'merged.dmp';
	
	my %tax_par; # $tax_par{taxid} = taxid parent
	my %tax_rank; # $tax_rank{taxid} = rank
	my %tax_ranked_lineage; # $ranked_lineage{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
							# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	my %merged; # $merged{taxid} = (list of old taxids)

	# $tax_ranked_lineage{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
	# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	# $tax_ranked_lineage{taxid}[0] scientific name of the taxon
	read_ncbitax_rankedlin($ncbi_tax_rankedlin, \%tax_ranked_lineage);
	# $tax_par{taxid} = taxid parent
	# $tax_rank{taxid} = rank
	read_new_nodes_dmp_to_hash($ncbi_tax_nodes, \%tax_par, \%tax_rank);

	read_new_merged_dmp_to_hash($merged_dump, \%merged);

	my $taxonomy = $outdir.'taxonomy.tsv';
	if($motif)# add motif to filename
	{
		$taxonomy =~ s/\.tsv/$motif.tsv/;
	}

	open(OUT, '>', $taxonomy) or die "Cannot open $taxonomy\n";
	print OUT "tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel\n";
	foreach my $taxid (sort {$a <=> $b} keys %tax_rank)
	{
		my $taxlevel = get_taxlevel($taxid, \%tax_ranked_lineage, \%tax_rank);
		if(exists $merged{$taxid})
		{
			foreach my $old_taxid (@{$merged{$taxid}})
			{
				print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$tax_ranked_lineage{$taxid}[0]	$old_taxid	$taxlevel\n";
			}
		}
		else
		{
			print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$tax_ranked_lineage{$taxid}[0]		$taxlevel\n";
		}
	}
	close OUT;
	return $taxonomy
}
#######################################################
sub modify_params_from_tags_orig
{
	my ($param, $inp) = @_;

	my @bad_tags = ();
	for(my $i = 0; $i<scalar@$inp; $i=$i+2)
	{
		$$inp[$i] =~ s/^-*//;
		if(exists $$param{$$inp[$i]})
		{
			$$param{$$inp[$i]} = $$inp[$i+1];
		}
		else
		{
			push(@bad_tags, $$inp[$i]);
		}
	}
	if(scalar @bad_tags > 0)
	{
		print "The following tags are not accepted: @bad_tags\n";
#		print_usage();
		exit;
	}
}

#######################################################
sub modify_params_from_tags
{
	my ($param, $inp) = @_;

	my @bad_tags = ();
	my $version = 0;
	my $help = 0;
	for(my $i = 0; $i<scalar@$inp; $i=$i+2)
	{
		$$inp[$i] =~ s/^-*//;
		if($$inp[$i] =~ /version/i)
		{
			$version = 1;
		}
		elsif($$inp[$i] eq 'h' or $$inp[$i] =~ /help/i)
		{
			$help = 1;
		}
		else
		{
			if(exists $$param{$$inp[$i]})
			{
				$$param{$$inp[$i]} = $$inp[$i+1];
			}
			else
			{
				push(@bad_tags, $$inp[$i]);
			}
		}
	}
	
	if($version)
	{
		print_version();
		exit;
	}
	
	if($help)
	{
		print_help();
		exit;
	}
	if(scalar @bad_tags > 0)
	{
		print "The following tags are not accepted: @bad_tags\n";
		print_help();
		exit;
	}
}

#######################################################

sub print_version
{
	print "####################\nmkCOInr-0.1.0\n";
	print "2022-05-15\n####################\n";
}
#######################################################
sub print_params_hash_to_log
{
	my ($params, $date) = @_;

	my @print = ('####', $date, 'PARAMETERS:');

	foreach my $param (sort keys %$params)
	{
		push(@print, "$param: $$params{$param}");
	}
	push(@print, '####');
	push(@print, '');

	return @print;
}

##############################################
sub print_stat
{
	my ($stat, $t) = @_;
	
	my @lines = ("\n####\n");
	foreach my $k (sort keys %$stat)
	{
		push(@lines, $k.$$stat{$k});
	}
	$t = time - $t;
	push(@lines, "Total runtime: ".$t." s");
	
	return join("\n", @lines);

}

############################################
sub read_fasta_to_hash
{
	my ($file) = @_;
	

	open(IN, $file) or die "Cannot open $file\n";
	my $id = '';
	my %hash;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		if($line =~ /^>([^\s]+)/)
		{
			$id = $1;
		}
		else
		{
			$hash{$id} .= uc $line;
		}
	}
	close IN;
	return %hash;
}

#########################################################
sub read_ncbitax_rankedlin
{
my ($file, $ranked_lin) = @_;
# $ranked_lin{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon

	unless(open(IN, $file))
	{
		print "Cannot open $file\n";
	}

	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split("\t", $line);
		my $taxid = shift@line;
		@{$$ranked_lin{$taxid}} = @line;
	}
	close IN;
}

#############################################
sub read_new_merged_dmp_to_hash
{
my ($file, $merged) = @_;
	# $merged{valid_taxid} = (old taxid list)


	open(IN, $file) or die "Cannot open $file\n";


	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split('\t', $line);;
		push(@{$$merged{$line[1]}}, $line[0]);
	}
	close IN;
}

#############################################
sub read_new_names_dmp_to_hash
{
my ($file, $name_taxids, $taxid_names) = @_;

	# $name_taxids{name} = (list of taxids) # homonymy
	# $taxid_names{taxid} = (list of all synonymes)

unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/\s$//;
	$line =~  s/\t\|//g;
	my @line = split("\t", $line);

	push(@{$$name_taxids{$line[1]}}, $line[0]);
	push(@{$$taxid_names{$line[0]}}, $line[1]);
}

close IN;
}

#############################################
sub read_new_names_dmp_to_hash_simple
{
my ($file, $taxid_names) = @_;

	# $taxid_names{taxid}{all names} = ''

unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/\s$//;
	$line =~  s/\t\|//g;
	my @line = split("\t", $line);

	$$taxid_names{$line[0]}{$line[1]} = '';
}

close IN;
}

#############################################
sub read_new_nodes_dmp_to_hash
{
my ($file, $par, $rank) = @_;

	# $tax_par{taxid} = taxid parent
	# $tax_par{taxid} = rank
	
unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/\s$//;
	$line =~  s/\t\|//g;
	my @line = split('\t', $line);
	
	$$par{$line[0]} = $line[1];
	$$rank{$line[0]} = $line[2];
}
close IN;
}

#############################################
	
sub read_new_taxidlineage_dmp_to_hash
{
	my ($file, $taxid_full_lineage_id) = @_;
# $taxid_full_lineage_id{taxid} = (list of taxids starting from the most distant)

	unless(open(IN, $file))
	{
		print "Cannot open $file\n";
	}

	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split("\t", $line);

		if(scalar @line > 1)
		{
			my @ids = split(' ', $line[1]);
			@{$$taxid_full_lineage_id{$line[0]}} = @ids;
		}
		else
		{
			@{$$taxid_full_lineage_id{$line[0]}} = ();
		}
	}
}

############################################
sub read_taxonomy_to_tax_hash
{
	my ($tax, $taxonomy) = @_;
	#$tax_{taxid} = (parent_tax_id	rank	name_txt	taxlevel_index)
	
	#tax_id	parent_tax_id	rank	name_txt	old_tax_id	(taxlevel_index)	Synonyms
	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		@{$$tax{$line[0]}} = ($line[1], $line[2], $line[3], $line[5]);
	}
	close IN;
}

############################################
sub read_taxonomy_and_merged
{
	my ($tax_name, $tax_rank, $tax_parent, $merged, $taxonomy) = @_;
	#my %tax_name; #%tax_name{taxid} = scientific name
	#my %tax_rank; # %tax_rank{taxid} = rank
	#my %tax_parent; # %tax_parent{taxid} = taxid parent
	#my  %merged{merged taxid} = new taxid
	
	# merged taxid are read into a separate hash, so later they can be replaced by the new taxid
	# other hashes do not contain the merged taxids
	
	#tax_id	parent_tax_id	rank	name_txt	old_tax_id	(taxlevel_index)	Synonyms
	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		$$tax_name{$line[0]} = $line[3];
		$$tax_rank{$line[0]} = $line[2];
		$$tax_parent{$line[0]} = $line[1];
		if(scalar @line > 4 and $line[4]) # old taxid exists
		{
			$$merged{$line[4]} = $line[0];
		}
	}
	close IN;
}

############################################

sub read_taxonomy_names
{
	my ($taxid_names, $name_taxids, $name_taxids_par, $taxonomy) = @_;
#	my %name_taxids; # $name_taxids{name}{$taxid} = '' # homonymy All names, including synonyms and homonyms
#	my %taxid_names; # $taxid_names{taxid}{names} = ''
#	#$name_taxids_par{$name}{taxid} = parent taxid # only scientifi names

	#tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel_index	synonyms
	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split("\t", $line);

		$$name_taxids{$line[3]}{$line[0]} = ''; #add taxid and scientific name
		$$taxid_names{$line[0]}{$line[3]} = '';
		$$name_taxids_par{$line[3]}{$line[0]} = $line[1]; # only scientific names
		if(scalar @line > 6) # there are synonyms
		{
			my @syns = split(';', $line[6]);
			foreach my $synonym (@syns)
			{
				$$name_taxids{$synonym}{$line[0]} = '';
				$$taxid_names{$line[0]}{$synonym} = '';
			}
		}
	}
	
	delete $$name_taxids{''}; # due to error in ncbi dmp, sometimes thara is not scientific name
	close IN;

}

############################################

sub reverse_complement
{
    my ($temp) = @_;
    my $revcomp = reverse $temp;
    $revcomp =~ tr/ATCGRYMKBHDVatcgrymkbhdv\[\]/TAGCYRKMVDHBtagcyrkmvdhb\]\[/;
    return $revcomp;
}

############################################################
sub sum
{
my $sum = 0;

foreach my $l (@_)
{
	$sum += $l;
}
return $sum;
}

1;
