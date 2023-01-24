use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT:
# dereplicated sequence file with taxIDs
#	seqID	taxID	sequence

# taxonomy.tsv (the most recent version for all ncbi and arbitrary taxids in it)
#	tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel


# AIM:
# Make databses in different formats
# Options:
#	blast
#	rdp #taxid*name*parent_taxid*taxlevel_index*taxlevel
#		#>KJ592630 cellularOrganisms;Eukaryota;undef_Eukaryota;undef_undef_Eukaryota;Florideophyceae;Hapalidiales;Hapalidiaceae;Mesophyllum;Mesophyllum_lichenoides
#	qiime
#	vtam
#	full

###############################
###############################
my %params = (
'tsv' => '', # seqID	taxID	sequence
'taxonomy' => '', # tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel
'outfmt' => '', # blast, rdp, qiime, full, vtam
'outdir' => '',
'out' => '',
'blast_path' => ''
);
modify_params_from_tags(\%params, \@ARGV);

my $tsv = $params{'tsv'};
my $taxonomy = $params{'taxonomy'};
my $outfmt = $params{'outfmt'};
my $outdir = $params{'outdir'};
my $out = $params{'out'};
my $blast_path = $params{'blast_path'};

$outdir = add_slash_to_dir($outdir);
$blast_path = add_slash_to_dir($blast_path);
###############################
###############################

####
#### define filenames
my $t = time;
my $t0 = time;
my $date = get_date();
my $tmpdir = $outdir.'temp/';
my $log = $outdir.'format_db.log';
my %stat;

unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";
####



if($outfmt eq 'blast' or $outfmt eq 'vtam')
{
	my $tmpdir = $outdir.'temp_'.$t.'/';
	unless(-e $tmpdir)
	{
		system 'mkdir -p '.$tmpdir;
	}
	####
	#### Make blast database 
	print "\n####\nMake blast database\n";
	print LOG "\n####\nMake blast database\n";
	
	my $fas = $tmpdir.'seq.fas';
	my $taxids = $tmpdir.'taxids.tsv'; # seqID	taxID
	open(FAS, '>', $fas) or die "Cannot open $fas\n";
	open(TID, '>', $taxids) or die "Cannot open $taxids\n";
	open(IN, $tsv) or die "Cannot open $tsv\n";
	my $title = <IN>;
	my %taxids; 
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		print FAS ">$line[0]\n$line[2]\n";
		print TID "$line[0]\t$line[1]\n";
		++$stat{'Number of sequences: '};
		$taxids{$line[1]} = '';
	}
	close IN;
	close FAS;
	close TID;
	$stat{'Number of unique taxids: '} = scalar keys %taxids;
	
	my $makeblastdb = $blast_path.'makeblastdb -dbtype nucl -in '.$fas.' -parse_seqids -taxid_map '.$taxids.' -out '.$outdir.$out;
	system $makeblastdb;
	
	my $taxonomy_out = $outdir.$out.'_taxonomy.tsv';
	if($outfmt eq 'vtam')
	{
		open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
		open(OUT, '>', $taxonomy_out) or die "Cannot open $taxonomy_out\n";
		while(my $line = <IN>)
		{
			$line =~ s/\s*$//;
			my @line = split("\t", $line);
			@line = splice(@line, 0, 6);
			print OUT join("\t", @line), "\n";
		}
		close OUT;
		close IN;
	}
	

	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	
	## clean up
	if(-e $tmpdir)
	{
		system 'rm -rf '.$tmpdir;
	}
}
elsif($outfmt eq 'rdp' or $outfmt eq 'qiime' or $outfmt eq 'full')
{
	### Read taxonomy
	print "\n####\nRead taxonomy\n";
	print LOG "\n####\nRead taxonomy\n";
	my %tax; #$tax{$taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	my %new_names; # $new_names{name} = taxid # new names are the undef_xxx names for the missing taxlevels
	read_taxonomy_local($taxonomy, \%tax);
	#####
	# temporary patch, since in taxonomy.tsv the taxID 0 is attributed to Acanthogyrus_cheni
		my $patch_taxid = 0;
		if(exists $tax{0} and $tax{0}[3] != 0) # taxid 0 is attributed something else then root
		{
			my $smallest_taxid = get_smallest_taxid(\%tax);
			$patch_taxid = $smallest_taxid -1;
			$tax{$patch_taxid} = $tax{0};
			delete $tax{0};
			foreach my $taxid (keys %tax)
			{
				if($tax{$taxid}[0] == 0) # replace the eventual taxparent 0 by  the new
				{
					$tax{$taxid}[0] = $patch_taxid;
				}
			}
#			print Dumper($tax{$patch_taxid});
		}
	####
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	
	if($out)
	{
#		$out .= '_';
	}


	#### Make $outfmt files
	print "\n####\nMake $outfmt files\n";
	print LOG "\n####\nMake $outfmt files\n";
	
	### get filenames, open files
	my $fas;
	my $taxids;
	my $full;
	if($outfmt eq 'full')
	{
		$full = $outdir.$out.'.tsv';
		open(FAS, '>', $full) or die "Cannot open $full\n";
		print FAS "seqID	taxon	taxID	taxlevel	superkingdom	superkingdom_taxID	kingdom	kingdom_taxID	phylum	phylum_taxID	class	class_taxID	order	order_taxID	family	family_taxID	genus	genus_taxID	species	species_taxID	sequence\n";
	}
	else
	{
		#### get filenames and open files
		# rdp >seqID cellularOrganisms;Eukaryota_2759;Metazoa_33208;Bryozoa_10205;Gymnolaemata_10206;Cheilostomatida_10207;Adeonidae_558780;Reptadeonella_2576536;Reptadeonella_violacea_-35055
		#	>seqID cellularOrganisms;Eukaryota_2759;Metazoa_33208;Arthropoda_6656;Malacostraca_6681;Decapoda_6683;Paguridae_6745;Paguridae_6745_genus;Paguridae_6745_genus_species
		# qiime >seqID
		$fas = $outdir.$out.'_trainseq.fasta';
		open(FAS, '>', $fas) or die "Cannot open $fas\n";
		##rdp taxid*name*parent_taxid*taxlevel_index*taxlevel
		##qiime1  339039	Bacteria;Proteobacteria;Alphaproteobacteria;Rhodospirillales;unclassified_Rhodospirillales
		##qiime2  229854	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__	
		$taxids = $outdir.$out.'_taxon.txt'; 
		open(TID, '>', $taxids) or die "Cannot open $taxids\n";
	}


	open(IN, $tsv) or die "Cannot open $tsv\n";
	my $title = <IN>;
	my %taxid_ranked_lin; # taxid_ranked_lin{taxid} = ranked lineage # list of taxids of the selected sequences
	my %tax_short; # $tax_short{taxid} = (name	parent_taxid	taxlevel_index	taxlevel) # taxids are the ones that appear in the ranked lineages in the selected sequences
	@{$tax_short{0}} = ('cellularOrganisms',0,0,'cellularOrganisms'); # initialize with root

	my $smallest_taxid = get_smallest_taxid(\%tax); # help to get arbitrary taxids for undef taxa in the ranked lineges; these taxids are not integrated to the %tax

	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		my $taxid = $line[1];
		
		##### patch to correct the fact that in taxonomy, the taxIµD 0 is attributed to non-root
		if($taxid == 0)
		{
			$taxid = $patch_taxid;
		}
		####

		unless(exists $taxid_ranked_lin{$taxid}) # taxid is new => get ranked lineage
		{
			my @lin = get_lineage_list($taxid, \%tax, 0); # list of taxids in the lineage
			my $ranked_lin = '';
			# identify major taxlevel in lineages, fill empty taxlevels with undef_higher_level_taxon, get taxid for them
			# make a hash with taxids used in the ranked lienages
			($ranked_lin, $smallest_taxid) = get_rdp_ranked_lin_from_tax(\@lin, \%tax, \%new_names, \%tax_short, $smallest_taxid, $outfmt);
			$taxid_ranked_lin{$taxid} = $ranked_lin;
			
		}
		if($outfmt eq 'rdp')
		{
			#>seqID cellularOrganisms;Eukaryota_2759;Metazoa_33208;Bryozoa_10205;Gymnolaemata_10206;Cheilostomatida_10207;Adeonidae_558780;Reptadeonella_2576536;Reptadeonella_violacea_-35055
			print FAS ">$line[0] $taxid_ranked_lin{$taxid}\n$line[2]\n"; 
		}
		elsif($outfmt eq 'qiime')
		{
			#>seqID 
			print FAS ">$line[0]\n$line[2]\n"; 
			##qiime2  229854	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
			print TID "$line[0]	$taxid_ranked_lin{$taxid}\n";
		}
		elsif($outfmt eq 'full')
		{
			print FAS "$line[0]	$tax{$taxid}[2]	$taxid	$tax{$taxid}[3]	$taxid_ranked_lin{$taxid}	$line[2]\n"; 
		}
		++$stat{'1. Number of sequences: '};
	}
	$stat{'2. Number of taxa: '} = scalar keys %taxid_ranked_lin;
	close IN;
	close FAS;
	
	#### Make taxonomy file
	##rdp taxid*name*parent_taxid*taxlevel_index*taxlevel
	if($outfmt eq 'rdp')
	{
		foreach my $taxid (sort keys %tax_short)
		{
			print TID $taxid, '*', join('*', @{$tax_short{$taxid}}), "\n";
			++$stat{'3. Number of taxids in taxonomy file: '};
		}
		close TID;
	}

	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####
}
else
{
	print "Un expected value for outmft. Select one of the followinf options: blast, rdp, qiime, full\n";
	exit;
}
####


	
print LOG print_stat(\%stat, $t0);
close LOG;


exit;


#############################################
sub get_rdp_ranked_lin_from_tax
{
	my ($lin, $tax, $new_names, $tax_short, $smallest_taxid, $outfmt) = @_;

	#$tax{$taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
	# @lin = (list of taxids stating from the smallest)
	# $tax_short{taxid} = (name	parent_taxid	taxlevel_index	taxlevel) # taxids are the one that apprear in the rankes lineages in the selected sequences

	
	my %ranks = (0, 'root',1, 'superkingdom',2, 'kingdom',3, 'phylum',4, 'class',5, 'order',6, 'family',7, 'genus',8, 'species');


	my @ranked_lin_taxid = (0, '', '', '', '', '', '', '', '');
	my @ranked_lin_name = ('cellularOrganisms', '', '', '', '', '', '', '', '');
	

	my %lin; # $lin{taxlevel_index} = taxid , some non-major level disappear, but we do not need them anyway
	foreach my $taxid (@$lin)
	{
		$lin{$$tax{$taxid}[3]} = $taxid;
	}
	
	# select the taxid of the major taxlevel
	if(0)
	{
		foreach my $taxid (@$lin)
		{
			my $taxlevel = $$tax{$taxid}[3];
			if($taxlevel == int($taxlevel) and $taxlevel >= 1) # taxlevel index is an integer
			{
				$ranked_lin_taxid[$taxlevel] = $taxid;
				$ranked_lin_name[$taxlevel] = $$tax{$taxid}[2].'_'.$taxid;
				$ranked_lin_name[$taxlevel] =~ s/ *\(.*//;
				$ranked_lin_name[$taxlevel] =~ s/ /_/g;
				@{$$tax_short{$taxid}} = ($ranked_lin_name[$taxlevel], $ranked_lin_taxid[$taxlevel-1], $taxlevel, $ranks{$taxlevel}); # complete hash for the rdp taxonomy file
			}
		}
	}
	
	
	# complete the taxelevels WO names and add new names to %tax and %new_names
	for(my $i = 1; $i < scalar @ranked_lin_taxid; ++$i)
	{
		if(exists $lin{$i}) # there is a taxid in the @lin with the correct taxlevel
		{
			my $taxid = $lin{$i};
			$ranked_lin_taxid[$i] = $taxid;
			if($outfmt eq 'full')
			{
				$ranked_lin_name[$i] = $$tax{$taxid}[2].';'.$taxid;
			}
			else
			{
				$ranked_lin_name[$i] = $$tax{$taxid}[2].'_'.$taxid;
			}
			$ranked_lin_name[$i] =~ s/ *\(.*//;
			$ranked_lin_name[$i] =~ s/ /_/g;
			@{$$tax_short{$taxid}} = ($ranked_lin_name[$i], $ranked_lin_taxid[$i-1], $i, $ranks{$i}); # complete hash for the rdp taxonomy file
		} 
		else
		{
			my $name = $ranked_lin_name[$i-1].'_'.$ranks{$i}; # add to previous_taxname the taxlevel
			$name =~ s/ /_/g;
			my $tid = 0;
			if(exists $$new_names{$name})
			{
				$tid = $$new_names{$name};
			}
			else
			{
				--$smallest_taxid;
				$tid = $smallest_taxid;
				$$new_names{$name} = $tid;
			}

			$ranked_lin_taxid[$i] = $tid;
			$ranked_lin_name[$i] = $name;
			@{$$tax_short{$tid}} = ($name, $ranked_lin_taxid[$i-1], $i, $ranks{$i}); # complete hash for the rdp taxonomy file
		}
	}
	
	my $ranked_lin;
	
	### make $ranked_lin if the appropiate format
	if($outfmt eq 'rdp')
	{
		$ranked_lin = join (';', @ranked_lin_name);
	}
	elsif ($outfmt eq 'qiime') # qiime or full #k__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__
	{
		my %temp_tl = (0 => '_kingdom', 1 => '_phylum', 2 => '_class',3 => '_order',4 => '_family',5 => '_genus',6 => '_species');
		my %temp_t = (0 => 'k__', 1 => 'p__', 2 => 'c__',3 => 'o__',4 => 'f__',5 => 'g__',6 => 's__');
		
		splice(@ranked_lin_name, 0, 2); # Delete the firts 2 levels above kingdom
		
		# add letters reffereing to the taxlevel
		for(my $i = 0; $i< scalar @ranked_lin_name; ++$i)
		{
			$ranked_lin_name[$i] = $temp_t{$i}.$ranked_lin_name[$i];
		}
		# delete last unvalid taxa
		for(my $i = 6; $i >0; --$i)
		{
			if($ranked_lin_name[$i] =~ /$temp_tl{$i}/) # last taxlevel is made from the name of a higher level taxa (e.g. Parnassius_species)
			{
				$ranked_lin_name[$i] = $temp_t{$i};
			}
			else # stop if a valid name has been found
			{
				last;
			}
		}
		
		$ranked_lin = join ('; ', @ranked_lin_name);
	}
	else # full
	{
		my %temp_tl = (0 => 'superkingdom', 1 => 'kingdom', 2 => 'phylum', 3 => 'class',4 => 'order',5 => 'family',6 => 'genus',7 => 'species');
		
		splice(@ranked_lin_name, 0, 1); # Delete the highest level (cellular organism)
		for(my $i = 7; $i >0; --$i)
		{
			if($ranked_lin_name[$i] =~ /$temp_tl{$i}$/) # last taxlevel is made from the name of a higher level taxa (e.g. Parnassius_species)
			{
				$ranked_lin_name[$i] = ';';
			}
			else # stop if a valid name has been found
			{
				last;
			}
		}
		$ranked_lin = join ("\t", @ranked_lin_name);
		$ranked_lin =~ s/;/\t/g; # separate taxids and taxon names
	}
	return ($ranked_lin,  $smallest_taxid);
}

#############################################

sub get_smallest_taxid
{
	my ($tax) = @_;
	
	my @l = sort {$a <=> $b} keys %$tax;
	return $l[0];
}

#############################################

sub get_largest_taxid
{
	my ($tax) = @_;
	
	my @l = sort {$a <=> $b} keys %$tax;
	return $l[-1];
}
#############################################
sub read_taxonomy_local
{
	my ($taxonomy, $tax) = @_;
	
	#$tax{$taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)

	open(IN, $taxonomy) or die "Cannot open $taxonomy\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//;
		my @line = split("\t", $line);
		@{$$tax{$line[0]}} = ($line[1], $line[2], $line[3], $line[5]);
	}
	close IN;

}

#############################################
sub print_help
{

print '
usage: perl format_db.pl [-options] -tsv INPUT_FILE -taxonomy TAXONOMY_TSV -outdir OUTDIR -out STRING -outfmt OUTPUT_FORMAT

 ARGUMENTS
   -tsv                    Input tsv file with seqID,taxID,sequence
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
   -outdir                 Name of the otput directory
   -out                    String to name the output files
   -outfmt                 [rdp/blast/qiime/full/vtam]; choose the format of the database
 OPTIONS
   -blast_path             Path to the blast executables; Not necessarry if it is in the PATH
',  "\n";
  exit;
}


