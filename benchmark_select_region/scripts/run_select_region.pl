use warnings;
use strict;
use mkdb;
use Data::Dumper;

# input 
# name of the dir with bait files 
# all parametres necessary for select_region

# ALGO
# run select_region for different bait files and clean up output

# Output
#

my %params = (
'dir' => '',
'tsv' => '', # sequence tsv with taxids
'outdir' => '',
'e_pcr' => 1, # 0/1 if 1 identify target sequences of the target region by e-pcr followed by vsearch --usearch_global 
'fw' => '', # GGNTGAACNGTNTAYCCNCC
'rv' => '', # TAWACTTCDGGRTGNCCRAARAAYCA
'trim_error' => 0.3,
'min_amplicon_length' => 100,
'max_amplicon_length' => 2000,
'min_overlap' => 10,
'bait_fas' => '', # if $e_pcr == 0 a fasta file with diverse sequences limited to the target region shoul be given as an input
'tcov' => 0.5, # Coverage in vsearch --usearch_global
'identity' => 0.7, # Min identity in vsearch --usearch_global
'cutadapt_path' => '',
'vsearch_path' => ''
);

modify_params_from_tags(\%params, \@ARGV);

my $dir = $params{dir}; # dir with bait files
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

$outdir = add_slash_to_dir($outdir);
$dir = add_slash_to_dir($dir);


unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

my @files = get_file_list_from_folder($dir, '.fas');


foreach my $bait (@files)
{
	my $out_tmp = $bait;
	$out_tmp =~ s/\..+//;
	my $tmp_name = $out_tmp;
	$out_tmp = $outdir.$out_tmp.'/';
	my $bait = $dir.$bait;
	my $cmd = 'perl /home/meglecz/mkCOInr/scripts/select_region.pl -tsv '.$tsv.' -outdir '.$out_tmp.' -target_region_fas '.$bait.' -e_pcr '.$e_pcr.' -fw '.$fw.' -rv '.$fw.' -trim_error '.$trim_error.' -min_amplicon_length '.$min_amplicon_length.' -max_amplicon_length '.$max_amplicon_length.' -min_overlap '.$min_overlap.' -tcov_hsp_perc '.$tcov_hsp_perc.' -perc_identity '.$perc_identity;
	system $cmd;
	
	system 'mv '.$out_tmp.'trimmed.tsv '.$outdir.$tmp_name.'_trimmed.tsv';
	system 'mv '.$out_tmp.'select_region.log '.$outdir.$tmp_name.'_select_region.log';
	system 'rm -rf '.$out_tmp;
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










