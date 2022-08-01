use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use Data::Dumper;

# input 
# tsv file with complete COI sequences from whole mitogenomes

# ALGO
# take one random sequence from each taxon of a given taxon level 

# Output
# fasta file with selected sequences
# tsv file with ranked lineage of all selected sequence and the total number of sequeces for the taxa for which it was selected

my %params = (
'coi_lin' => '/home/meglecz/mkCOInr/no_git/benchmark_select_region/positive_dataset/mitogenomes_COI/2_check_sequences/coi_sequences_checked.tsv', 
'taxrank' => 'phylum',
'n_bait' => 10, # make $n_bait bait files 
'n_per_taxon' => 1, # number of sequneces chosen randomly from each taxon
'outdir' => '/home/meglecz/mkCOInr/no_git/benchmark_select_region/bait_phylum/'
);
modify_params_from_tags(\%params, \@ARGV);

my $coi_lin = $params{coi_lin};
my $taxrank = $params{taxrank};
my $n_bait = $params{n_bait};
my $n_per_taxon = $params{n_per_taxon};
my $outdir = $params{outdir};

$outdir = add_slash_to_dir($outdir);

my %taxlevel = ('root' => 3, 'superkingdom' => 4, 'kingdom' => 5, 'phylum' => 6, 'class' => 7, 'order' => 8, 'family' => 9, 'genus' => 10, 'species' => 11);
my $taxrank_ind = $taxlevel{$taxrank};

unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

my %taxa; # @taxa{taxon at a  taxrank_ind taxlevel}{seqid} = @line
open(IN, $coi_lin) or die "Cannot open $coi_lin\n";
my $t = <IN>;
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	$line =~ s/ /_/g;
	my @line = split("\t", $line);
	@{$taxa{$line[$taxrank_ind]}{$line[0]}} = @line;
}
close IN;

for(my $i = 0; $i < $n_bait; ++$i)
{
	my $out = $outdir.$taxrank.'_'.$i.'.fas';

	open(OUT, '>', $out) or die "Cannot open $out\n";
	foreach my $taxon (sort keys %taxa)
	{
		my @seqids = sort keys %{$taxa{$taxon}};
		my $n = scalar @seqids;
		my %pos; 
		if($n > $n_per_taxon) # if more sequnece for a taxon than n_per_taxon
		{
			# make list of random numbers 
			while(scalar keys %pos < $n_per_taxon)
			{
				my $rn = random_number_1($n);
				$pos{$rn} = '';
			}
		}
		else # take all seq of the taxon
		{
			for(my $i = 0; $i < $n; ++$i)
			{
				$pos{$i} = '';
			}
		}
		# write selected sequences
		foreach my $p (sort {$a <=> $b} keys %pos)
		{
			my $seqid = $seqids[$p];
			my @info = @{$taxa{$taxon}{$seqid}};
			print OUT ">$seqid $info[1];$info[5];$info[6];$info[7];$info[8];$info[9];$info[10];$info[11]\n$info[2]\n";
#			print OUT join("\t", @{$taxa{$taxon}{$seqid}}), "\n";
		}

	}
	close OUT;
}

exit;

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










