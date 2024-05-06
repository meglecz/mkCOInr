use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# select metadata of BOLD sequences that appear in the final file of COInr.tsv

# INPUT
# 	BOLD public datapackage in CSV format
# 	COIn.tsv (seq_id	sequence)

# AIM
# 	select metadata of BOLD sequences that appear in the final file of COInr

#OUTPUT
# BOLD_metadata: tsv file with all valid sequences downloaded from BOLD and reatined in COInr, and their metadata for authorship tracability.

###############################
###############################

my %params = (
	'input_metadata' => '',
	'coinr' => '',
	'out' => ''
);
modify_params_from_tags(\%params, \@ARGV);

my $input_metadata = $params{'input_metadata'};
my $coinr = $params{'coinr'};
my $out = $params{'out'};


####
# Read BOLD sequnec IDs to hash
####
open(IN, $coinr) or die "Cannot open $coinr\n";
my %coi; # $coi{BOLD ids} = ''
while(my $line = <IN>)
{
	my @line = split("\t", $line);
	if($line[0] =~ /^BOLD_/)
	{
		$coi{$line[0]} = '';
	}
}
close IN;


####
# Select lines from BOLD data package that are in COInr and print them to $out
####
open(IN, $input_metadata) or die "Cannot open $input_metadata\n";
my $title = <IN>;
my $id_ind = get_indices('processid', $title); # get the index of the column with processid
my $seq_ind = get_indices('nucraw', $title); # get the index of the column sequences; delete sequnece, since it takes place and it is already in COInr
my $marker_ind = get_indices('marker_code', $title); # get the index of the column with marker code; specimenid is the same for different sequneces from different markers, chosse the one that corresponds to COI

$title =~ s/\tnucraw//;
$title = 'COInrID'."\t".$title;
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT $title;
while(my $line = <IN>)
{
	my @line = split("\t", $line);
	my $processid = $line[$id_ind];
	my $marker = $line[$marker_ind];
	my $id = 'BOLD_'.$marker.'_'.$processid;
	if(exists $coi{$id})
	{
		splice(@line, $seq_ind, 1); # delete sequnece to save spave. It is in COInr.tsv
		$line = join("\t", @line);
		print OUT $id, "\t", $line;
	}
}
close IN;
close OUT;


###################################

sub get_indices
{
	my ($name, $line) = @_;
	# find indice in line that corresponds to name

	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split("\t", $line);
	for(my $i = 0; $i < scalar @line; ++$i)
	{
		if($name eq $line[$i])
		{
			return $i;
		}
	}
	return '';
}

###########################################

sub print_help
{

print '
usage: perl format_bold_package.pl -input_metadata BOLD_PACKAGE_TSV_FILE -out OUTFILE -coinr COINR_TSV

 ARGUMENTS
   -input_metadata         Name of input data package downloaded from BOLD https://www.boldsystems.org/index.php/datapackages
   -coinr                  COInr.tsv TSV file with sequenceID and sequences
   -out                    Name of the otput file
',  "\n";
  exit;
}
