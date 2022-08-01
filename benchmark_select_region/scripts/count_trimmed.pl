use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use Data::Dumper;

my %params = (
'dir' => '',
'motif' => '',
'e_pcr' => 1
);
modify_params_from_tags(\%params, \@ARGV);

my $dir = $params{dir};
my $motif = $params{motif};
my $e_pcr = $params{e_pcr};
$dir = add_slash_to_dir($dir);
my $out = $dir.'count_trimmed_';
if($e_pcr)
{
	$out .= 'epcr.tsv';
}
else
{
	$out .= 'bait.tsv';
}

my @folders = get_file_list_from_folder($dir, $motif);

open(OUT, '>', $out) or die "Cannot open $out\n";
if($e_pcr)
{
	print OUT "neg/pos/	dataset	option	trim_error	min_amplicon_length	identity	trimmed	untrimmed	trimmed_prop	untrimmed_prop\n";
}
else
{
	print OUT "neg/pos/	dataset	bait_taxlevel	N_seq_per_taxon	identity	replicate	trimmed	untrimmed	trimmed_prop	untrimmed_prop\n";
}
foreach my $folder (@folders)
{
	my $name = $folder;
	my @name = split('_', $name);
	$folder = $dir.$folder.'/';
	
	my $cmd = 'wc -l '.$folder.'trimmed.tsv';
	my $trimmed = `$cmd`;
	$trimmed =~ s/ .*\s*$//;
	--$trimmed;
	
	
	$cmd = 'wc -l '.$folder.'untrimmed.tsv';
	my $untrimmed = `$cmd`;
	$untrimmed =~ s/ .*\s*$//;
	--$untrimmed;
	
	my $all = $trimmed + $untrimmed;
	print OUT join("\t", @name), "	$trimmed	$untrimmed	", $trimmed / $all, "\t", $untrimmed / $all, "\n";
}
close OUT;

exit;



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
