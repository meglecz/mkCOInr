use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
# RDP raw training data downloaded from https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo18_rawtrainingdata.zip

# This contains
#	trainset18_db_taxid.txt
#		taxid*name*parent_id*taxlevel_index*taxlevel (taxlevel_index tells place of the taxon in the lineage, it does not correspobds to a fix taxrank)
#		7*Acidimicrobiaceae*5*5*family
#		8*Acidimicrobium*7*6*genus
#
#	rdp trainset18_xxxxx_speciesrank.fa
#		>Mycobacterium heidelbergense str. 2554/91 Type	domain__Bacteria; phylum__Actinobacteria; class__Actinobacteria; order__Mycobacteriales; family__Mycobacteriaceae; genus__Mycobacterium

# AIM
#	make a taxonomy.tsv file using the rdp_taxids (does not correspond to NCBI taxids)
#	make a sequence tsv file

# OUTPUT
# -rdp_sequences.tsv file (seqID	TaxID	sequence)
# -rdp_taxonomy.tsv taxonomy file in mkCOInr format (taxid	parent_tax_id	rank	name_txt	old_tax_id	taxlevel)

###############################
###############################
my %params = (
'rdp_tax' =>  '/home/emese/DB/RDP_classifier/RDPClassifier_16S_trainsetNo18_rawtrainingdata/trainset18_db_taxid.txt', # taxid*name*parent_id*taxlevel_index*taxlevel
'fasta' => '/home/emese/DB/RDP_classifier/RDPClassifier_16S_trainsetNo18_rawtrainingdata/trainset18_062020_speciesrank.fa', #Mycobacterium heidelbergense str. 2554/91 Type	domain__Bacteria; phylum__Actinobacteria; class__Actinobacteria; order__Mycobacteriales; family__Mycobacteriaceae; genus__Mycobacterium
'outdir' =>  '/home/emese/DB/RDP_classifier/vtam_16S/', 
);

modify_params_from_tags(\%params, \@ARGV);

my $rdp_tax = $params{rdp_tax};
my $fasta = $params{fasta};
my $outdir = $params{outdir};

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
my $log = $outdir.'format_rdp.log';
my %stat;

my $taxonomy = $outdir.'rdp_taxonomy.tsv'; # taxonomy file: tax_id  parent_tax_id   rank    name_txt        old_tax_id      taxlevel        synonyms
my $sequence_tsv = $outdir.'rdp_sequences.tsv'; 

my %taxlevel_ind = ('rootrank', 0, 'domain', 1, 'phylum', 3, 'class', 4, 'subclass', 4.5, 'order', 5, 'suborder', 5.5, 'family', 6, 'genus', 7);

#####
#####
# Read taxonomy file to %tax
#####
#####
my %tax; # $tax{taxid} = (parent_tax_id   rank    name_txt        old_tax_id      taxlevel)
open(IN, $rdp_tax) or die "Cannot open $rdp_tax\n";
#7*Acidimicrobiaceae*5*5*family
#8*Acidimicrobium*7*6*genus
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	$line =~ s/\*/;/g;
	my @line = split(';', $line);
	if($line[4] eq 'rootrank') # root must have 1 as taxid and taxpar
	{
		@{$tax{1}} = (1, 'no rank', 'root', '', 0);
		next;
	}
	else
	{
		my $taxid = $line[0] +1; ### Since bacteria has 1 as taxid in input, add 1 to each taxid
		my $taxid_par = $line[2] +1;
#		my $taxid = $line[0]; 
#		my $taxid_par = $line[2];
		@{$tax{$taxid}} = ($taxid_par, $line[4], $line[1], '', $taxlevel_ind{$line[4]}); 
	}
}
close IN;
#print Dumper(\%tax);


#####
#####
# Write taxonomy file, make %taxaname
#####
#####
open(OUT, '>', $taxonomy) or die "Cannot open $taxonomy\n";
print OUT "tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel\n";
my %taxname; # taxaname{name}{taxrank} = taxid
foreach my $taxid (sort {$a <=> $b} keys %tax)
{
	print OUT "$taxid	", join("\t", @{$tax{$taxid}}), "\n";
	my $name = $tax{$taxid}[2];
	my $rank = $tax{$taxid}[1];
	if(exists $taxname{$name}{$rank})
	{
		print "$name exists more than once in taxonomy file\n";
	}
	else
	{
		$taxname{$name}{$rank} = $taxid;
	}
}
close OUT;
#print Dumper(\%taxname);




#####
#####
# Find taxid for each acc
#####
#####
#>AJ000684	Mycobacterium heidelbergense str. 2554/91 Type	domain__Bacteria; phylum__Actinobacteria; class__Actinobacteria; order__Mycobacteriales; family__Mycobacteriaceae; genus__Mycobacterium
#gaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtaatctgccctgcacatcgggataagcctgggaaactgggtctaataccgaataggacctcgaggcgcatgccttgtggtggaaagcttttgcggtgtgggatgggcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggncttcgggttgtaaacctctttcagcagggacgaagcgcaagtgacggtacctgcagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaaaccgggggcttaaccctcggcgtgcgggcgatacgggcagactggagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgccggcagagatgtcggttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcgggtaatgccggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcaaggttaagcgaatccttttaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaacctttgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgt
open(IN, $fasta) or die "Cannot open $fasta\n";
my %seq; # $seq{acc} = seq
my %taxids; # $seq{acc} = taxid
my $acc;
my $bool = 0;
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	if($line =~ />([^\s]+)\t.+\t(domain.+)$/)
	{
		$acc = $1;
		my $lin = $2;
		if($acc =~ /^[0-9\.]+$/) # numerical seqID, this causes pb with older BLAST versions
		{
			$acc = 'X'.$acc; # add a letter before numerical seqID
		}

		my @lin = split('; ', $lin);
		my $taxon = pop @lin;
		my $taxid = 0;
		if($taxon =~ s/genus__//)
		{
			if(exists $taxname{$taxon}{genus})
			{
				$taxids{$acc} = $taxname{$taxon}{genus};
				$bool=1;
			}
			else
			{
				print "WARNINIG: No taxid for $acc	$taxon\n";
				$bool=0;
			}
		}
		else
		{
				print "WARNINIG: Smallest taxon is not at the genus level: $acc	$taxon\n";
				$bool=0;
		}
	}
	elsif($line =~ />/)
	{
		print "WARNINIG: Fasta title line is not a correct format: $line\n";
		$bool=0;
	}
	elsif($bool)
	{
		if(exists $seq{$acc})
		{
			print "ERROR: $acc exists already as aseqID\n";
			exit;
		}
		$seq{$acc} = uc $line;
	}
}
close IN;


#####
#####
# Print sequneces tsv with taxids
#####
#####
open(OUT, '>', $sequence_tsv) or die "Cannot open $sequence_tsv\n";
print OUT "seqID	taxID	sequence\n";
foreach my $acc (sort keys %taxids)
{
	print OUT "$acc	$taxids{$acc}	$seq{$acc}\n";
}
close OUT;

exit;


#############################################
sub print_help
{

print '
usage: perl format_rdp.pl [-options] -rdp_tax INPUT_FILE -fasta INPUT_FASTA -outdir OUTDIR

 ARGUMENTS
   -rdp_tax                  Taxon file in RDP taxonomy format
   -fasta                    Fasta file in RDP fasta format
   -outdir                   Name of the otput directory
',  "\n";

  exit;
}
