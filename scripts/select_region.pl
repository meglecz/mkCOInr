use warnings;
use strict;
use mkdb;
use Data::Dumper;


my %params = (
'tsv' => '', # sequence tsv with taxids
'outdir' => '',
'e_pcr' => 1, # 0/1 if 1 identify target sequences of the target region by e-pcr followed by vsearch --usearch_global 
'fw' => '', # GGNTGAACNGTNTAYCCNCC
'rv' => '', # TAWACTTCDGGRTGNCCRAARAAYCA
'trim_error' => 0.3,
'min_amplicon_length' => 100,
'max_amplicon_length' => 2000,
'min_overlap' => 10,
'target_region_fas' => '', # if $e_pcr == 0 a fasta file with diverse sequences limited to the target region shoul be given as an input
'tcov_hsp_perc' => 0.5, # Coverage in vsearch --usearch_global
'perc_identity' => 0.7, # Min identity in vsearch --usearch_global
'cutadapt_path' => '',
'vsearch_path' => ''
);

modify_params_from_tags(\%params, \@ARGV);

my $tsv = $params{tsv}; # sequence tsv with taxids
my $outdir = $params{outdir};
my $e_pcr = $params{e_pcr}; # 0/1 if 1 identify target sequences of the target region by e-pcr followed by vsearch --usearch_global 
my $fw = $params{fw};
my $rv = $params{rv};
my $trim_error = $params{trim_error};
my $min_amplicon_length = $params{min_amplicon_length};
my $max_amplicon_length = $params{max_amplicon_length};
my $min_overlap = $params{min_overlap};
my $target_region_fas = $params{target_region_fas}; # if there is a fasta file with diverse sequences limited to the target region, the e-pcr can be avoided
my $tcov_hsp_perc = $params{tcov_hsp_perc}; # Coverage in vsearch --usearch_global
my $perc_identity = $params{perc_identity}; # Min identity in vsearch --usearch_global
my $cutadapt_path = $params{cutadapt_path};
my $vsearch_path = $params{vsearch_path};

my $delete_tmp =0;

#### define filenames and variables
$outdir = add_slash_to_dir($outdir);
my $t = time;
my $t0 = time;
my $date = get_date();
my %stat;
my $tmpdir = $outdir.'tmp_'.$t.'/';
unless (-e $tmpdir)
{
	system 'mkdir -p '.$tmpdir;
}
my $input_fas = $tmpdir.'input.fasta';
my $log = $outdir.'select_region.log';
open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";

#my $out_trimmed =  $outdir.'cutadapt_trimmed.fas';

# selected regions
# region 26-350 TNTCNACNAAYCAYAARGAYATTGG (jgLCO1490)- GGNTGAACNGTNTAYCCNCC (III-C-R-revN) $expected_amplicon_length = 324;
# region 371-683 IIICBR GGNTGAACNGTNTAYCCNCC (III-C-R-revN) - TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) $expected_amplicon_length = 313;
# region 710-1458  IIICBR TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) - GGGCAGCCRTGRATTCAYTC (H8121) $expected_amplicon_length = 749;

# tested regions 
# region 47-350 GAYATTGGWACHTTATAYTTHATHTTHGG (Zfdg2)- GGNTGAACNGTNTAYCCNCC (III-C-R-revN) $expected_amplicon_length = 304;
# region 710-1153 (1155)  IIICBR TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) - GCWGTMTTTGCTATTATAGCAGG (F2640) $expected_amplicon_length = 446;
# region 710-1149  IIICBR TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) - ATAGGRGCWGTATTTGCYATTAT (C1 J2636) $expected_amplicon_length = 442;
# region 710-1523  IIICBR TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) - TYCATTGCACTAATCTGCCATATTAG (R3037) $expected_amplicon_length = 814;
# region 710-883  IIICBR TAWACTTCDGGRTGNCCRAARAAYCA (HBR2d) - GTRGCNGAYGTRAARTATGCTCG (COI907aH2) $expected_amplicon_length = 174;
# region 795-883  IIICBR ATTCCTATGTAGCCGAATGGTTCTTT (AWCR6) - GTRGCNGAYGTRAARTATGCTCG (COI907aH2) $expected_amplicon_length = 89;
# region 906-1458  IIICBR GTRGCNGAYGTRAARTATGCTCG (COI907aH2) - GGGCAGCCRTGRATTCAYTC (H8121) $expected_amplicon_length = 553;

# https://www.researchgate.net/figure/Universal-Primer-Sequences-for-COI_tbl1_271027042
# https://www.researchgate.net/publication/271027042_Molecular_Species_Identification_of_Forensically_Important_Flies_in_Korea/link/56307db608ae76226de04fd8/download


#### Make fasta file from input tsv
my %untrimmed; # %untrimmed when trimmig sequences they are deleted from this hash
my %trimmed; # $trimmed{seqid} = trimmed seq; %trimmed will be gradualy filled
if(1)
{
	print "\n####\nMake fasta file from input tsv\n";
	print LOG "\n####\nMake fasta file from input tsv\n";
	make_fasta($tsv, $input_fas, \%untrimmed);
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}
####


my $cluster_in = $target_region_fas;
my $usearch_in = $input_fas;
my $trimmed_cutadapt = $outdir.'cutadapt_trimmed.fas';
my $untrimmed = $tmpdir.'untrimmed_cutadapt.fas';

if($e_pcr) # make an e_pcr on the input sequences to get sequences limited to the target region
{
	#### e-PCR with cutadapt
	print "\n####\ne-PCR with cutadapt\n";
	print LOG "\n####\ne-PCR with cutadapt\n";
	$rv = reverse_complement($rv);
	my $cmd = $cutadapt_path.'cutadapt --trimmed-only --report=minimal --cores=0 -O '.$min_overlap.' -e '.$trim_error.' --no-indels --minimum-length '.$min_amplicon_length.' --maximum-length '.$max_amplicon_length.' -g "'.$fw.'...'.$rv.'" --output '.$trimmed_cutadapt.' '.$input_fas;
	print $cmd, "\n";
	system $cmd;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####


	#### Make a fasta with sequences untrimmed by cutadapt
	print "\n####\nMake a fasta with sequences untrimmed by cutadapt\n";
	print LOG "\n####\nMake a fasta with sequences untrimmed by cutadapt\n";
	my $n = scalar keys %untrimmed;
	%trimmed = read_fasta_to_hash($trimmed_cutadapt);
	
	open(OUT, '>', $untrimmed) or die  "Cannot open $untrimmed\n";
	foreach my $sid (keys %untrimmed)
	{
		if(exists $trimmed{$sid})
		{
			delete $untrimmed{$sid}; # delete trimmed sequence from %untrimmed
		}
		else
		{
			print OUT ">$sid\n$untrimmed{$sid}\n";
		}
	}
	close OUT;
	$stat{'1. Number of sequences in the input fasta file: '} = $n;
	$stat{'2. Number of sequences in the cutadapt trimmed file: '} =  scalar keys %trimmed;
	$stat{'3. Number of sequences in the cutadapt untrimmed file: '} = scalar keys %untrimmed;

	#### redifine input file names for cluster and usearch
	$usearch_in = $untrimmed;  # use this file as an input to vsearch
	$cluster_in = $trimmed_cutadapt; # use this file as an input to cluster
	
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}
####


my $centroids = $outdir.'target_centroids.fas';
if(1)
{
	#### Cluster trimmed sequneces or input target_region_fas
	print "\n####\nCluster trimmed sequneces or input target_region_fas\n";
	print LOG "\n####\nCluster trimmed sequneces or input target_region_fas\n";
	my $vsearch = $vsearch_path.'vsearch --cluster_fast '.$cluster_in.' --centroids '.$centroids.' --id 0.90';
	system $vsearch;
	# get sequences count in $centroids
	my $cmd = 'grep ">" '.$centroids.' | wc -l';
	my $temp = `$cmd`;
	chomp $temp;
	$stat{'4. Number of centroids: '} = $temp;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####
}

#### vsearch --usearch_global 
my $usearch_out = $tmpdir.'usearch_out.tsv';
my $outfmt = 'query+target+id+alnlen+qilo+qihi+tilo+tihi+tcov+qstrand+aln';
if(1)
{
	print "\n####\nvsearch --usearch_global sequences against centroids\n";
	print LOG "\n####\nvsearch --usearch_global sequences against centroids\n";
	my $cmd = 'vsearch --usearch_global '.$usearch_in.' --db '.$centroids.' --userout '.$usearch_out.' --userfields '.$outfmt.' --id '.$perc_identity.' --strand both --target_cov '.$tcov_hsp_perc;
	system $cmd;
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	####
}

#### trim usearch_global hits
print "\n####\nTrim usearch_global hits\n";
print LOG "\n####\nTrim usearch_global hits\n";
open(IN, $usearch_out) or die "Cannot open $usearch_out\n";
my $c = 0;
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	my @line = split('\t', $line);
	if(exists $untrimmed{$line[0]})
	{
		my $sequence = $untrimmed{$line[0]};
		if($line[9] eq '-')
		{
			$sequence = reverse_complement($sequence);
		}
		my $l = $line[5] - $line[4] +1;
		$sequence = substr($sequence, $line[4]-1, $l);
		$trimmed{$line[0]} = $sequence;
		delete $untrimmed{$line[0]};
		++$c;
	}
	else
	{
		print LOG "$line[0] at least twice in the vsearch output file\n";
	}
}
close IN;
$stat{'5. Number of sequences trimmed based on usearch_global: '} = $c;
$stat{'6. Number of sequences untrimmed after usearch_global: '} = scalar keys %untrimmed;

my $vsearch_trimmed = $outdir.'trimmed.tsv';
my $vsearch_untrimmed = $outdir.'untrimmed.tsv';
open(IN, $tsv) or die "Cannot open $tsv\n";
my $title = <IN>;
open(T, '>', $vsearch_trimmed ) or die "Cannot open $vsearch_trimmed \n";
print T "seqID	taxID	sequence\n";
open(UT, '>', $vsearch_untrimmed ) or die "Cannot open $vsearch_untrimmed \n";
print UT "seqID	taxID	sequence\n";
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	my @line = split("\t", $line);
	if(exists $trimmed{$line[0]})
	{
		$line[2] = $trimmed{$line[0]};
		print T join("\t", @line), "\n";
		++$stat{'7. Total number of trimmed sequences: '};
	}
	elsif(exists $untrimmed{$line[0]})
	{
		print UT $line, "\n";
		++$stat{'8. Total number of untrimmed sequences: '};
	}
	else
	{
		print "WARNING: $line[0] is neither in trimmed nor in untrimmed hash\n";
	}
}
close T;
close UT;
close IN;

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
#####

if($delete_tmp)
{
	system 'rm -f '.$tmpdir.'*';
	system 'rmdir '.$tmpdir;
}

print LOG print_stat(\%stat, $t0);
close LOG;
exit;


###############################################################

sub make_fasta
{
	my ($tsv, $fas, $seq) = @_;
	
	open(FAS, '>', $fas) or die "Cannot open $fas\n";
	open(IN, $tsv) or die "Cannot open $tsv\n";
	my $title = <IN>;
	my %double;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		if(scalar @line < 3)
		{
			print "ERROR: Input tsv file must have at least 3 columns, with sequences in the third\n ";
			exit;
		}
		if(exists $$seq{$line[0]})
		{
			$double{$line[0]} = '';
		}
		$$seq{$line[0]} = $line[2];
		print FAS ">$line[0]\n$line[2]\n";
	}
	close FAS;
	close IN;
	
	if (scalar keys %double)
	{
		print "ERROR: The following seqIDs are not unique:\n";
		print join("\n", keys %double), "\n";
		print LOG "ERROR: The following seqIDs are not unique:\n";
		print LOG join("\n", keys %double), "\n";
	}
}

############################################

sub make_reg_exp
{
	my ($seq) = @_;
	
	my %code = ('A' => 'A','C' => 'C','G' => 'G','T' => 'T','M' => '[AC]','R' => '[AG]','W' => '[AT]','S' => '[CG]','Y' => '[CT]','K' => '[GT]','V' => '[ACG]','H' => '[ACT]','D' => '[AGT]','B' => '[CGT]','N' => '[GACT]');
		
	$seq = uc $seq;
	my @l = split('', $seq);
	for(my $i = 0; $i < scalar @l; ++$i)
	{
		$l[$i] = $code{$l[$i]};
	}
	my $reg_exp = join('', @l);
	return $reg_exp;
}

#########################################################################

sub print_params
{
	my ($file, $beg, $end) = @_;
	
	open(IN, $file) or die "Cannot open $file\n";
	my $bool = 0;
	while(my $line = <IN>)
	{
		if($line =~ $beg)
		{
			$bool = 1;
		}
		elsif($line =~ $end)
		{
			print LOG $line, "\n\n";
			close IN;
			last;
		}
		if($bool)
		{
			print LOG $line;
		}
	}
	close IN;

}

############################################

sub print_help
{

print '
usage: perl select_region.pl -tsv INPUT_TSV -outdir OUTDIR

 ARGUMENTS
   -tsv                    Input tsv file with seqID,taxID,sequence
   -outdir                 Name of the otput directory
 OPTIONS/PARMATERS
   -target_region_fas      Optional; Can be produced by E-pcr included in the script 
                              Fasta file with sequences already trimmed to the target region
  E-pcr related parameters
   -e_pcr                  [0/1]; If 1, identify the target region of the sequences by e-pcr
   -fw                     Optional; Forward primer that amplifies the target region
   -rv                     Optional; Reverse primer that amplifies the target region
   -trim_error             Real [0-1]; Default : 0.3
                              The proportion of mismatches allowed between the primer and the sequence during the e_pcr
   -min_amplicon_length    Integer; Default: 100
                              The minimum length of the amplicon after primer trimming
   -max_amplicon_length    Integer; Default: 2000
                              The maximum length of the amplicon after primer trimming
   -min_overlap            Integer; Default: 10
                              The minimum overlap between primer and the sequence during e-pcr
   -cutadapt_path          Path to the cutadapt executables; Not necessarry if it is in the PATH
   
  usearch_global related parameters
   -tcov_hsp_perc          Real [0-1]; Default : 0.5
                              Minimum coverage of the target sequence in usearch_global hits
   -perc_identity          Real [0-1]; Default : 0.7
                              Minimum percentage of identity between the sequence and the target in usearch_global hits
   -vsearch_path           Path to the vsearch executables; Not necessarry if it is in the PATH

', "\n";
  exit;
}
