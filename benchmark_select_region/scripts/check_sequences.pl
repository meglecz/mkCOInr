use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use Data::Dumper;

# input 
# tsv file with complete COI sequences from whole mitogenomes (out format_ncbi.pl)

# ALGO
# Complete lines with lineage
# delete sequences bellow and above lenght threshold

# Output
# tsv file with ranked lineage of all selected sequences

my %params = (
'tsv' => '/home/meglecz/mkCOInr/no_git/benchmark_select_region/positive_dataset/mitogenomes_COI/1_format_ncbi/ncbi_sequences.tsv', # seqid	taxid	seq
'outdir' => '/home/meglecz/mkCOInr/no_git/benchmark_select_region/positive_dataset/mitogenomes_COI/2_check_sequences/', # seqid	taxid	seq
'taxonomy' => '/home/meglecz/makeCOIdb/mkCOInr_2022-05-10/COInr_zenodo/taxonomy.tsv', # needed only if there is no ranked ineage file make_lineage =1
'taxrank' => 'phylum',
'min_length' => 1100,
'max_length' => 2000,
'seq_n' => 9999999999999999999999999999999 # get random sequneces if more sequneces in iput than seq_n
);
modify_params_from_tags(\%params, \@ARGV);

my $tsv = $params{tsv};
my $outdir = $params{outdir};
my $taxonomy = $params{taxonomy};
my $taxrank = $params{taxrank};
my $min_length = $params{min_length};
my $max_length = $params{max_length};
my $seq_n = $params{seq_n};

$outdir = add_slash_to_dir($outdir);

my %taxlevel = ('root' => 0, 'superkingdom' => 1, 'kingdom' => 2, 'phylum' => 3, 'class' => 4, 'order' => 5, 'family' => 6, 'genus' => 7, 'species' => 8);
my $taxrank_ind = $taxlevel{$taxrank};
my $out_random = $outdir.'random_sample.tsv';
my $out_ok = $outdir.'sequences_checked.tsv';
my $out_len = $outdir.'sequences_length_pb.tsv';

unless ( -e $outdir)
{
	system 'mkdir -p '.$outdir;
}


my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
read_taxonomy_to_tax_hash(\%tax, $taxonomy);

open(IN, $tsv) or die "Cannot open $tsv\n";
open(OUT, '>', $out_ok) or die "Cannot open $out_ok\n";
print OUT "seqID	taxID	sequence	root	superkingdom	kingdom	phylum	class	order	family	genus	species	seqlen\n";
open(LEN, '>', $out_len) or die "Cannot open $out_len\n";
print LEN "seqID	taxID	sequence	root	superkingdom	kingdom	phylum	class	order	family	genus	species	seqlen\n";
my $t = <IN>;
my %l;

my $sn = 0; # number of sequences in checked_file
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	my @line = split("\t", $line);
	my $seqid = $line[0];
	my $taxid = $line[1];
	my $seq = $line[2];
	my $l = length $seq;
	
	my @full_lineage =  get_lineage_list($taxid, \%tax, 0);
	my @ranked_lin = make_ranked_lineage(\@full_lineage, \%tax);
	
	++$l{length $seq};
	
	if(($l < $min_length) or ($l > $max_length)) # length out of range
	{
		print LEN "$line	", join("\t", @ranked_lin), "	$l\n";
	}
	else # clean sequence
	{
		print OUT "$line	", join("\t", @ranked_lin), "	$l\n";
		++$sn;
	}
}
close IN;
close OUT;
close LEN;

if($sn > $seq_n)
{
	open(IN, $out_ok) or die "Cannot open $out_ok\n";
	my @data = <IN>;
	close IN;
	my $t = shift @data;
	my $nlines = scalar @data;
	
	# get seq_n ranom numbers bellow nlines
	my %pos;
	while (scalar keys %pos < $seq_n)
	{
		my $ln = random_number_1($nlines);
		$pos{$ln} = '';
	}
	
	# print out randomly selected lines
	open(OUT, '>', $out_random) or die "Cannot open $out_random\n";
	print OUT $t;
	for(my $i = 0; $i < $nlines; ++$i)
	{
		if(exists $pos{$i})
		{
			print OUT $data[$i];
		}
	}
	close OUT
}

foreach my $l (sort {$a <=> $b} keys %l)
{
#	print "$l	$l{$l}\n";
}


exit;

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


##############################################################
sub random_number_1
{
#If number begins with 0s, zeros are deleted; random number is smaller than number
my ($number) =@_;
	my $digit = length $number;
	my $ok = 1;
	while ($ok == 1)
	{
		my $random_numb = '';
		for (my $i = 0; $i < $digit; ++$i)
		{
			my $l = int rand 10;
#			print "$l ";
		$random_numb .= $l;
 	#end for $i
		}
		if ($random_numb < ($number))
		{
			$random_numb =~ s/^0+//;
			if ($random_numb eq '')
			{
				$random_numb = 0;
			}
			$ok = 0;
			return $random_numb;	
		}
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

