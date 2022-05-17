use warnings;
use strict;
use mkdb;
use Data::Dumper;

# INPUT: 
#	cds fasta file (created by nsdpy)
#	TaxIDs.txt (created by nsdpy)
#	taxonomy.txt (0_download_taxonomy.pl)
#		Tab separated columns: tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel

# ALGO:
# 	Check if gene name and protein name correspondrs to COI
# 	Eliminate genes if they have intron
#	Can have more than 1 COI gene in the same mithochondrion
#	Accept only if valid taxid, replace old non-valid taxids by up to date taxids
#	Eliminate sequences from environmental and metagenomic samples
#	If $check_name is activated take the taxID of the smallest taxon with a valid latin name, otherwise keep the original taxid 
#	Sequences are already in a correct orientation in the input file, since that are comming from CDS files

#	clean sequences
#		Upper case
#		replace non ATCG by N
#		Delete gaps and external Ns
#		Delete sequence if more than max_n conscutive Ns
#	Keep sequneces if length is between min_length and max_length


# OUTPUT: 
#	ncbi_sequences_date.txt (seqID	taxID	sequence)
#	logfile

###############################
###############################
my %params = 
(
	'cds' => '', #output of nsdpy -r "request" -T --cds
	'taxids' => '', #output of nsdpy -r "request" -T --cds; # MT044302.1  273913
	'taxonomy' => '', # created by 0_download_taxonomy.pl
	'check_name' => 1, # if 1, keep only taxa vith valid latin name format
	'max_n' => 5, # while clean_seq, eliminate sequences with $max_n or more Ns
	'min_length' => 100, # minimum length of the cleaned sequence; 100 for COI, 1000 from COI from mitogenome
	'max_length' => 2000,  # maximum length of the cleaned sequence; 2000 for both COI and mitogenome
	'outdir' => ''
);
modify_params_from_tags(\%params, \@ARGV);

my @regexp_gene = ('^\"*cox*-*_*[1i]\"*$'),
my @regexp_protein = ('cytochrome', '^cox*-*_*[1i]$'), # case insensitive

my $cds = $params{'cds'};
my $tax = $params{'taxids'};
my $taxonomy = $params{'taxonomy'};
my $check_name = $params{'check_name'};
my $max_n = $params{'max_n'};
my $min_length = $params{'min_length'};
my $max_length = $params{'max_length'};
my $outdir = $params{'outdir'};

$outdir = add_slash_to_dir($outdir);
###############################
###############################

#### define filenames
my $t = time;
my $t0 = time;
my $date = get_date();
my $log = $outdir.'format_ncbi.log';
my $out = $outdir.'ncbi_sequences.tsv';

unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}
#### 

open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";

#### Read taxonomy file
print LOG "\n####\nRead taxonomy file\n";
print "\n####\nRead taxonomy file\n";
my %tax_name; #%tax_name{taxid} = scientific name
my %tax_rank; # %tax_rank{taxid} = rank
my %tax_parent; # %tax_parent{taxid} = taxid parent
my %merged; # %merged{merged taxid} = new taxid
read_taxonomy_and_merged(\%tax_name, \%tax_rank, \%tax_parent, \%merged, $taxonomy);
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
#### 


#### Read TaxID file
# replace merged taxids by new one
print LOG "\n####\nRead TaxID file\n";
print "\n####\nRead TaxID file\n";
open(IN, $tax) or die "Cannot open $tax\n";
my %acc_tax; # $acc_tax{accession} = taxid
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split(' ', $line);
	if(exists $tax_rank{$line[1]}) # valid taxid
	{
		$acc_tax{$line[0]} = $line[1];
	}
	elsif(exists $merged{$line[1]}) # old taxid => replace it with new
	{
		$acc_tax{$line[0]} = $merged{$line[1]};
	}
	else
	{
#		print "****$line[0]****$line[1]****\n";
	}
}
close IN;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
#### 


#### Read CDS file
print LOG "\n####\nRead CDS file\n";
print "\n####\nRead CDS file\n";
# >lcl|MT044302.1_cds_QTW90623.1_1 [gene=ND1] [protein=NADH dehydrogenase subunit 1] [transl_except=(pos:955..956,aa:TERM)] [protein_id=QTW90623.1] [location=2746..3701] [gbkey=CDS]
# >lcl|MT326627.2_cds_QXF78160.1_3 [gene=COX1] [protein=cytochrome c oxidase subunit I] [protein_id=QXF78160.1] [location=5517..7073] [gbkey=CDS]
open(IN, $cds) or die "Cannot open $cds\n";
my %seq; # $seq{$seq_id} = seq
my $seq_bool = 0;
my $acc = 0;
my %copies; # $copies{acc} = number of copies;
my $seqid = ''; # accession_copy_n
my %stat;
#open(GENE, '>', 'gene.txt') or die "Cannot open gene.txt\n";
#open(PROT, '>', 'prot.txt') or die "Cannot open prot.txt\n";
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	if($line =~ /^>/)# new sequence
	{
		++$stat{'01. Total number of input sequences: '};
		$seq_bool = 0;
		# get info; this is a standard order for fields
		if($line =~ /^>([^\s]+).*\[gene\=([^\]]+)\].*\[protein\=([^\]]+).*\[location\=([^\]]+)/) 
		{
			my $id = $1;
			my $gene = $2;
			my $protein = $3;
			my $location = $4;
			### not a correct gene name
			unless(check_regexp($gene, \@regexp_gene, 0)) 
			{
#				print GENE $gene, "\n";
				next;
			}
			++$stat{'02. Total number of sequences with regular expression in gene name: '};
			### not a correct protein name
			unless(check_regexp($protein, \@regexp_protein, 0)) 
			{
#				print PROT $protein, "\n";
				next;
			}
			++$stat{'03. Total number of sequences with regular expression in protein name: '};
			
			### intron in gene
			if($location =~ /join/) 
			{
				next;
			}
			++$stat{'04. Total number of sequences without intron: '};
			
			$acc = $id;
			$acc =~ s/lcl\|//;
			$acc =~ s/_cds_.*//;
			### no Taxid in TaxIDs.txt and taxonomy
			unless(exists $acc_tax{$acc}) 
			{
				next;
			}
			++$stat{'05. Total number of sequences with valid TaxID: '};
			
			
			### check environemental sequences 
			my $taxid = $acc_tax{$acc};
			if($tax_name{$taxid} =~ /environmental/ or $tax_name{$taxid} =~ /metagenom/)
			{
				next;
			}
			++$stat{'06. Total number of sequences with valid TaxID exluding environmental and metagenomic sequences: '};
			
			
			### get_lowest_correct_latin_name
			if($check_name) 
			{
				$taxid = get_lowest_correct_latin_name($taxid, \%tax_rank, \%tax_name, \%tax_parent);
				if($taxid == 1)# no correct taxon name
				{
					next;
				}
				++$stat{'07. Total number of sequences with correct latin name: '};
			}
			
			$acc_tax{$acc} = $taxid; # redefine taxid, for the taxid of the smallest correct latin name in %acc_hash
			
			### redefine seqid to accession_copynumber
			++$copies{$acc}; # get the copy number of the gene
			$seqid = $acc;
			$seqid =~ s/\..*$//; # delete version from accession
			$seqid = $seqid.'_'.$copies{$acc}; # add copy number to acession
			$acc_tax{$seqid} = $taxid; # add seqid to hash
			$seq_bool = 1;
		}
	}
	elsif($seq_bool) # get sequence 
	{
		$seq{$seqid} .= $line;
	}
}
close IN;
#close GENE;
#close PROT;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

#### Clean sequences and eliminate exact matches within taxID
print LOG "\n####\nClean sequences and eliminate exact matches\n";
print "\n####\nClean sequences and eliminate exact matches\n";
# $acc_tax{accession} = taxid
# $seq{$seq_id} = seq
my %tax_seq; #$tax_seq{taxid}{seq} = seqid
foreach my $id (sort keys %seq)
{
	my $seq = clean_seq($seq{$id}, $max_n);
	my $taxid = $acc_tax{$id};
	if($seq) # Sequence passed the cleaning
	{
		++$stat{'08. Total number of sequences passed sequence cleaning: '};
		my $len = length $seq{$id};
		if($len >= $min_length and $len <= $max_length)
		{
			++$stat{'09. Total number of sequences with length in the correct range: '};
			$tax_seq{$taxid}{$seq} = $id;
		}
	}
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####


#### Print output file
print LOG "\n####\nPrint output file\n";
print "\n####\nPrint output file\n";
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "seqID	taxID	sequence\n";
foreach my $taxid (sort {$a<=>$b} keys %tax_seq)
{
	foreach my $seq (keys %{$tax_seq{$taxid}})
	{
		print OUT "$tax_seq{$taxid}{$seq}	$taxid	$seq\n";
		++$stat{'10. Total number of sequences after eliminating exact matches: '};
	}
}
$stat{'11. Total number of taxids after eliminating exact matches: '} = scalar keys %tax_seq;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

print LOG print_stat(\%stat, $t0);
exit;

#####################################################


sub check_regexp
{
	my ($str, $regexp_list, $case_sensitive) = @_;
	
	foreach my $r (@$regexp_list)
	{
		if($case_sensitive)
		{
			if($str =~ /$r/)
			{
				return 1;
			}
		}
		else
		{
			if($str =~ /$r/i)
			{
				return 1;
			}
		}
	}
	return 0;
}

###########################################

sub print_help
{

print '
usage: perl format_ncbi.pl [-options] -cds CDS_FASTA -taxids TAXIDS -taxonomy TAXONOMY_TSV -outdir OUTDIR

 ARGUMENTS
   -cds                    Input CDS fasta file
   -taxids                 Input tsv file: seqID,taxID
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms
   -outdir                 Name of the otput directory
 OPTIONS
   -check_name             [0/1]; Default 1
                               If one keeps only taxa with valid Latin name format
   -max_n                  Positive integer; Default:5)
                              Eliminates sequences with max_n or more consecutive Ns
   -min_length             Positive integer; Default:100
                              Minimum length of the cleaned sequence;
   -max_length             Positive integer; Default:2000
',  "\n";
  exit;
}


