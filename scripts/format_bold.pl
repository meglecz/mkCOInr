use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
#	download_dir: Name of the directory wiht tsv files downloaded from bold (output of 0_download_bold). There should not be oter tsv files in the directory


# AIM

# Clean downloaded files and pool information to lineage and sequence files
#		Eliminate partial lines (mostly errors in the database)
#		Select sequences with a given marker list
#		Clean sequences: correct sequence Ids, gaps deleted, non-TCGA changed to N, external Ns deleted, sequences with more than $max_n consecutive Ns are deleted
#		Clean lineages: Keep only names matching a correct latin name format (^[A-Z][a-z-]+\s[a-z-]+$/ # species level, /^[A-Z][a-z-]+$/ All other levels)
#		Pool identical lineages into one line with the list  of valid sequence Ids in the last field
#		Eliminate lines with environmental and metagenomic samples
#		Keep sequences with length in a $min_length and max_length range 

# Orient sequence
#	Count codon STOPs (TAA, TAG) in each reading frame
#		0 if codon stop in all frames OR stand + and - among the frames WO codon stop => AMBIGUOUS
		# +1 or -1 otherwise
#	make a small "reference" db form randomly sampled oriented sequences
#	blast ambiguous sequences to check orientation
#	eliminate sequences sequences WO hit
#	NOTE: Almost all sequences are correctly oriented in BOLD. Trying to orient sequneces require  quite a lot of runtime,  and the method is not fully tested. If you intend to use the database for BLAST based assignment, the corrrect orientation is not important, so you can skip this step


#OUTPUT
# BOLD_sequences: tsv file with all valid sequences downloaded from BOLD 
#	seqid	sequence
# BOLD_lineages: all identical lineages are pooled into a same line
#	phylum	class	order	family	subfamily	genus	species	seqIDs
# BOLD_partial_lines => partial lines 


###############################
###############################

my %params = (
	'download_dir' => '',
	'outdir' => '',
	'marker_list' => 'COI-5P COI-3P',
	'check_name' => 1,
	'max_n' => 5, # while clean_seq, eliminate sequences with $max_n or more consecutive Ns
	'min_length' => 100, # minimum length of the cleaned sequence; 
	'max_length' => 2000,  # maximum length of the cleaned sequence; 
	'check_orientation' => 0,
	'blast_path' => ''
);
modify_params_from_tags(\%params, \@ARGV);

my $download_dir = $params{'download_dir'};
my $outdir = $params{'outdir'};
my $marker_list = $params{'marker_list'};
my $check_name = $params{'check_name'};
my $max_n = $params{'max_n'};
my $min_length = $params{'min_length'};
my $max_length = $params{'max_length'};
my $check_orientation = $params{'check_orientation'};
my $blast_path = $params{'blast_path'};

$outdir = add_slash_to_dir($outdir);
$download_dir = add_slash_to_dir($download_dir);

###############################
###############################

#### define filenames and variables
unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

my $t = time;
my $t0 = time;
my $date = get_date();
my $bold_lineages = $outdir.'bold_lineages.tsv'; 
my $seq_file = $outdir.'bold_sequences.tsv'; 
my $partial_lines = $outdir.'bold_partial_lines.txt'; # No clear download error, but some lines are too short => delete from the dataset, and documented here
my $ambiguous_orientation = $outdir.'bold_ambiguous_orientation.fas'; # No clear download error, but some lines are too short => delete from the dataset, and documented here
my $log = $outdir.'format_bold.log';
open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";

my %markers = make_marker_hash($marker_list);
my %stat;
my @taxlevels = qw(phylum_name class_name order_name family_name subfamily_name genus_name species_name);
my @fields = @taxlevels;
push(@fields, 'sequenceID');
push(@fields, 'markercode');
push(@fields, 'nucleotides');
####


#### Read tsv filenames
my @files = get_file_list_from_folder($download_dir, '\.tsv');
####


#### clean sequences, make lineages, eliminate identical sequences of the same lineage
print "\n####\nClean sequences, make lineages, eliminate identical sequences of the same lineage\n";
print LOG "\n####\nClean sequences, make lineages, eliminate identical sequences of the same lineage\n";
open(SH, '>', $partial_lines) or die "Cannot open $partial_lines\n";
my %lineage; # $lineage{join(';', @lineage)} = (seqids)
my %sids; # $sids{seqid} = ''; # to make sure that the the seq id is unique 
my %lineage_seq; # $lineage_seq{lineage}{seq} = seqid # eliminate identical sequences 
my %full_marker_list; #%full_marker_list{marker} = count
foreach my $file (@files)
{
	my $taxon = $file;
	$file = $download_dir.$file;
	$taxon =~ s/\.tsv//;
	
	++$stat{'01. Number of input files: '};

	my @data = read_data($file);
	#### check if file is empty
	$stat{'02. Number of input files without sequneces: '} = 0;
	unless(scalar @data) # empty file
	{
		++$stat{'02. Number of taxa without sequneces: '};
		print LOG "$taxon\tEmpty file\n";
		next; # go to next file
	}

	# $fields_hash{col name} = col ind
	my %fields_hash = get_field_indices(\@fields, \@data);
	# occasionally some lines are too short; error in input data => print these lines to BOLD_partial lines.txt
	delete_partial_lines(\@data, "\t", $file); 
	
	#### clean lines, and pool data
	foreach my $line (@data)
	{
		$line =~ s/"//g;
		my @line = split("\t", $line);
		my $marker = $line[$fields_hash{markercode}];
		++$full_marker_list{$marker};
		if(exists $markers{$marker}) # marker is OK
		{
			++$stat{'03. Number of lines with correct marker: '};
			my $seqid = $line[$fields_hash{sequenceID}];
			if($seqid =~ /^[0-9]+$/) # SequenceID is a positive integer
			{
				if(exists $sids{$seqid}) # if more than one line with the same SeqID (probably error in downloading) takes only the first line
				{
					next;
				}
				++$stat{'04. Number of lines with correct sequence id: '};

				### clean and check sequence
				my $seq = clean_seq($line[$fields_hash{nucleotides}], $max_n);
				if(length $seq <= $min_length or length $seq >= $max_length)
				{
					$seq = '';
				}
				unless($seq) # bad quality sequence
				{
					next;
				}
				++$stat{'05. Number of lines with correct sequence: '};
				### order and clean lineage
				my @new_lineage; # (phylum_name class_name order_name family_name subfamily_name genus_name species_name)
				foreach my $taxon_level (@taxlevels)
				{
					my $taxon = $line[$fields_hash{$taxon_level}]; 
					if($taxon =~ /environmental/ or $taxon =~ /metagenome/) # one of the taxon in lineage contains envinomental or metagenome
					{
						@new_lineage = ();
						last;
					}
					push(@new_lineage, $taxon);
				}
				if(scalar @new_lineage) # Not an environmental sample
				{
					++$stat{'06. Number of lines without environmental or metagenomic sequences: '};
					if($check_name)
					{
						clean_tax_names(\@new_lineage); # keep only taxon names with format of valid names
													# works only if the last element is a species. If not code should be adjusted
					}
					my $lin = join("\t", @new_lineage);
					if($lin =~ /[^\t]/)# at least one valid taxon name in lineage
					{
						++$stat{'07. Number of lines with at least one correct taxon name: '};
						$sids{$seqid} = ''; # keep seqids in a hash
						# keep sequence in hash
						$lineage_seq{$lin}{$seq} = $seqid; # if identical sequneces for the same linegae, keep just one
					}
				}
			}# seqid
		}# marker
	}# foreach line
}
close SH;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

# print full marker list in bold files
if(0)
{
	foreach my $marker (sort keys %full_marker_list)
	{
			print "$marker	$full_marker_list{$marker}\n";
	}
}

#### Count sequences after eliminating exact matches
print "\n####\nCount sequences after eliminating exact matches\n";
print LOG "\n####\nCount sequences after eliminating exact matches\n";
foreach my $lin (sort keys %lineage_seq)
{
	foreach my $seq (sort keys %{$lineage_seq{$lin}})
	{
		++$stat{'08. Number of sequences after eliminating exact matches within taxa : '};
	}
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;


#### Orient sequences
if($check_orientation)
{
	print "\n####\nOrient sequences\n";
	print LOG "\n####\nOrient sequences\n";
	# make temp dir
	my $tmpdir = $outdir.'tmp/';
	unless(-e $tmpdir)
	{
		system 'mkdir -p '.$tmpdir;
	}

	print "\n##\nCheck orientation with codons stops\n";
	print LOG "\n##\nCheck orientation with codons stops\n";
	#$ambiguous_orientation{$lin}{$seq} with sequences with ambiguous orientation: # if codon stop in all frames OR stand + and - among the frames WO codon stop
	# delete amb sequences from %lineage_seq
	# rev comp sequences in %lineage_seq if necessary
	my %revcomp_codon;
	my %ambiguous_orientation = orient_with_codon_stop(\%lineage_seq, \%stat, \%revcomp_codon);
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;

	print "\n##\nMake a small reference file from correctly oriented sequences\n";
	print LOG "\n##\nMake a small reference file from correctly oriented sequences\n";
	# make a reference file by taking the longest sequence of randomly sampled lineages in the %lineage_seq
	my $ref_file = subsample_seq(\%lineage_seq, 10000, $tmpdir);
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;

	print "\n##\nOrient ambiguous sequences by BLAST\n";
	print LOG "\n##\nOrient ambiguous sequences by BLAST\n";
	my $e = 1e-20,
	my $task =  'megablast';
	my %revcomp_blast;
	orient_with_blast(\%lineage_seq, $ref_file, \%ambiguous_orientation, $tmpdir, $blast_path, $e, $task, \%stat, \%revcomp_blast, $ambiguous_orientation);
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
	
#	my $revcomp_codon = $outdir.'reverse_complemented_based_on_codon_stop_orientation.fas';
#	my $revcomp_blast = $outdir.'reverse_complemented_based_on_blast_orientation.fas';
#	print_seq($revcomp_codon, \%revcomp_codon);
#	print_seq($revcomp_blast, \%revcomp_blast);
	
	system 'rm '.$tmpdir.'*';
	system 'rmdir '.$tmpdir;
}


#### print sequences
print "\n####\nPrint sequences\n";
print LOG "\n####\nPrint sequences\n";
open(SEQ, '>', $seq_file) or die "Cannot open $seq_file\n";
print SEQ "seqID	sequence\n";
foreach my $lin (sort keys %lineage_seq)
{
	foreach my $seq (sort keys %{$lineage_seq{$lin}})
	{
		print SEQ $lineage_seq{$lin}{$seq}, "\t", $seq, "\n";
		++$stat{'16. Number of sequences in output: '};
	}
}
close SEQ;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

#### print lineages
print "\n####\nPrint lineages\n";
print LOG "\n####\nPrint lineages\n";
open(OUT, '>', $bold_lineages) or die "Cannot open $bold_lineages\n";
my $tls= join("\t", @taxlevels);
$tls =~ s/_name//g;
print OUT "$tls\tseqIDs\n";
#$lineage_seq{$lin}{$seq}
foreach my $lin (sort keys %lineage_seq)
{
	my @sids;
	foreach my $seq (keys %{$lineage_seq{$lin}})
	{
		push(@sids, $lineage_seq{$lin}{$seq});
	}
	print OUT $lin, "\t", join(';', @sids), "\n";
}
$stat{'17. Number of unique lineages in output: '} = scalar (keys %lineage_seq);
close OUT;

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
print LOG print_stat(\%stat, $t0);
close LOG;

exit;

###############################

sub print_seq
{
	my ($out, $hash) = @_;
	#$hash{$lin}{$seq}  = seqid
	
	open(OUT, ">", $out) or die "Cannnot open $out\n";
	foreach my $l (sort keys %$hash)
	{
		foreach my $s (sort keys %{$$hash{$l}})
		{
			print OUT ">$$hash{$l}{$s}\n$s\n";
		}
	}
	close OUT;
}

###############################

sub make_marker_hash
{
	my ($marker_list) = @_;
	my %markers;
	
	my @l = split(' ', $marker_list);
	
	foreach my $m (@l)
	{
		$markers{$m} = '';
	}
	return %markers;
}
###############################
sub orient_with_codon_stop
{
	my ($lineage_seq, $stat, $revcomp) = @_;
	#$lineage_seq{$lin}{$seq}
	#$revcomp{$lin}{$seq} # just for tracking, these sequences are also in the lineage_seq
	
	my %stop = ('TAA', '', 'TAG', '');
#	my %stop_2 = ('TAA', '', 'TAG', '', 'AGA', '', 'AGG', '');#Vertebrata
#	my %stop_16 = ('TAA', '', 'TGA', '');# Chlorophyceae, Spizellomyces punctatus
#	my %stop_22 = ('TCA', '', 'TAA', '', 'TGA', ''); # Scenedesmus obliquus
#	my %stop_23 = ('TTA', '', 'TAA', '', 'TAG', '', 'TGA', ''); # Thraustochytrium aureum
#	my %stop_33 = ('TAG', ''); # Cephalodiscidae (Hemichordata)
	
	my %amb; #$amb{$lin}{$seq} with ambigous sequences
	# Make a hash with sequences with ambiguous orientation: # 0 if codon stop in all frames OR stand + and - among the frames WO codon stop
	
	foreach my $lin (keys %$lineage_seq)
	{
		foreach my $seq (keys %{$$lineage_seq{$lin}})
		{
			my $strand = check_codon_stop(\%stop, $seq); # 0 if codon stop in all frames OR stand + and - among the frames WO codon stop
														# +1 if only stand + with 0 codon stop, -1 if only stand - with 0 codon stop
			if($strand == -1) # reverse complement sequences and correct info in the hash
			{
				my $rc_seq = reverse_complement($seq);
				$$lineage_seq{$lin}{$rc_seq} = $$lineage_seq{$lin}{$seq};
				$$revcomp{$lin}{$rc_seq} = $$lineage_seq{$lin}{$rc_seq};
				delete $$lineage_seq{$lin}{$seq};
				++$$stat{'10. Number of sequences in reverse orientation according to codon stop search: '};
			}
			elsif($strand == 0) # make  fasta file with ambiguos sequences and delete from %lineage_seq
			{
				$amb{$lin}{$seq} = $$lineage_seq{$lin}{$seq};
				delete $$lineage_seq{$lin}{$seq};
				++$$stat{'11. Number of sequences with ambiguous orientation according to codon stop search: '};
			}
			else
			{
				++$$stat{'09. Number of sequences correctly oriented according to codon stop search: '};
			}
		}
		unless(scalar keys %{$$lineage_seq{$lin}}) # no sequences left in lineage
		{
			delete $$lineage_seq{$lin};
		}
	}
	return %amb;
}
###############################

sub orient_with_blast
{
	my ($lineage_seq, $ref_file, $amb, $tmpdir, $blast_path, $e, $task, $stat, $revcomp, $ambiguous_orientation) = @_;

	#$amb{$lin}{$seq}
	#$lineage_seq{$lin}{$seq}
	#$revcomp{$lin}{$seq}  # just for tracking, these sequences are also in the lineage_seq

	# make fasta file with sequences with ambiguous orientation 
	my $ambigue = $tmpdir.'ambiguous_orientation.fas';
	open(OUT, '>', $ambigue) or die "Cannot open $ambigue\n";
	my %sid; #sid{seq id} = (seq, $lineage);
	foreach my $lin (keys %$amb)
	{
		foreach my $seq (keys %{$$amb{$lin}})
		{
			print OUT ">$$amb{$lin}{$seq}\n$seq\n";
			@{$sid{$$amb{$lin}{$seq}}} = ($seq, $lin);
		}
	}
	close OUT;
	my $outfile = $tmpdir.'megablast.tsv';
	blast($ambigue, $ref_file, $blast_path, $outfile, $e, $task);
	
	open(IN, $outfile) or die "Cannot open $outfile\n";
	# read info from BLAST file to hash
	#my $outfmt = '6 qseqid sseqid pident qstart qend sstart send qcovs sstrand';
	my %hit; # $hit{seqid} = (qseqid sseqid pident qstart qend sstart send qcovs sstrand)
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		unless(exists $hit{$line[0]}) # if multiple hits, take only the first (best Evalue)
		{
			my $seq = $sid{$line[0]}[0];
			my $lin = $sid{$line[0]}[1];
			if($line[8] eq 'minus')
			{
				$seq = reverse_complement($seq);
				$$revcomp{$lin}{$seq} = $line[0];
				++$$stat{'13. Number of sequences in reverse orientation according to BLAST: '};
			}
			else
			{
				++$$stat{'12. Number of sequences in correct orientation according to BLAST: '};
			}
			$$lineage_seq{$lin}{$seq} = $line[0];
			delete $$amb{$lin}{$seq};
		}
		$hit{$line[0]}= '';
	}
	close IN;
	
	my %missing_lineage;
	open(OUT, '>', $ambiguous_orientation) or die "Cannot open $ambiguous_orientation\n";
	foreach my $lin (keys %$amb)
	{
		foreach my $seq (keys %{$$amb{$lin}})
		{
			print OUT ">$$amb{$lin}{$seq} $lin\n$seq\n";
			++$$stat{'14. Number of sequences with ambiguous orientation according to BLAST: '};
		}
		unless(exists $$lineage_seq{$lin}) # lineage is present in ambiguous, but not in oriented lineages
		{
			$missing_lineage{$lin} = '';
		}
	}
	close OUT;
	print LOG "Lineages without oriented sequences:\n";
	$$stat{'15. Number of lineages lost in orientation: '} = 0;
	foreach my $lin (keys %missing_lineage)
	{
		print LOG "$lin\n";
		++$$stat{'15. Number of lineages lost in orientation: '};
	}
}

#####################################################

sub blast
{
	my ($fasta, $db, $blast_path, $outfile, $e, $task) = @_;

	my $outfmt = '6 qseqid sseqid pident qstart qend sstart send qcovs sstrand';
	my $dust = 'yes';
	my $num_threads =8;
	my $max_target_seqs = 1;
	
	mkblastdb($db, $blast_path);
	my $blast = $blast_path.'blastn -task '.$task.' -db "'.$db.'" -query "'.$fasta.'" -evalue '.$e.' -out "'.$outfile.'" -outfmt "'.$outfmt.'" -dust '.$dust.' -num_threads '.$num_threads.' -max_target_seqs '.$max_target_seqs;
	system $blast;

}
#####################################################
sub mkblastdb
{
	my ($db, $blast_path) = @_;

	my $mkdb = $blast_path.'makeblastdb -in '.$db.' -dbtype nucl'; 
	system $mkdb;
}

###############################

sub subsample_seq
{
	my ($lineage_seq, $seqn, $tmpdir) = @_;
	# make a reference file by taking the longest sequence of randomly sampled lineages in the %lineage_seq
	
	my @lineages = keys %$lineage_seq;
	my $lineage_n = scalar @lineages;
	# check if enough lineages in the hash for subsampling WO replacement
	
	if($seqn > $lineage_n)
	{
		print "Not enough linegaes in linegae hash. The number of lineages is reduced to ", int($lineage_n/10), "\n";
		$seqn = int($lineage_n/10);
	}
	
	my %select_lin;
	for(my $i = 0; $i <$seqn; ++$i)
	{
		my $r = int(rand($lineage_n));
		while(exists $select_lin{$lineages[$r]})
		{
			$r = int(rand($lineage_n));
		}
		$select_lin{$lineages[$r]} = '';
	}

	my $ref = $tmpdir.'ref.fas';
	open(OUT, '>', $ref) or die "Cannot open $ref\n";
	foreach my $lin (keys %select_lin)
	{
		my $max_length = 0;
		my $select_seq = '';
		foreach my $seq (keys %{$$lineage_seq{$lin}}) # take the longest sequences from each selected lineage
		{
			if((length $seq) > $max_length)
			{
				$select_seq = $seq;
			}
		}
		print OUT ">$$lineage_seq{$lin}{$select_seq}\n$select_seq\n";
	}
	close OUT;
	return $ref;
}

###############################

sub check_codon_stop
{
	my ($stop, $seq) = @_;
	
	my %count_stop; #$count_stop{number of codon stop} = (frames)
	
	for(my $i = 0; $i < 3; ++$i) # count the codont stop in each frame +
	{
		my $c = 0;
		for(my $j = $i; $j < (length $seq) -2; $j=$j+3) # frame
		{
			my $codon = substr($seq, $j, 3);
			if(exists $$stop{$codon})
			{
				++$c;
			}
		}
		push(@{$count_stop{$c}}, $i+1);
	}
	
	$seq = reverse_complement($seq);
	for(my $i = 0; $i < 3; ++$i) # count the codont stop in each frame -
	{
		my $c = 0;
		for(my $j = $i; $j < (length $seq) -2; $j=$j+3) # frame
		{
			my $codon = substr($seq, $j, 3);
			if(exists $$stop{$codon})
			{
				++$c;
			}
		}
		push(@{$count_stop{$c}}, ($i+1)*-10); # frame -10, -20,-30
	}
	
	my @stop_n = sort {$a <=> $b} keys %count_stop;
	my $min_stop_n = $stop_n[0];
	my $somme_min_frame = sum(@{$count_stop{$min_stop_n}}); # sum the frames (1,2,3,-10,-20,-30) that has the lowest number of codon stops

	if($min_stop_n == 0) # there is at least one frame WO codon stop
	{
		if ($somme_min_frame > 0) # only frame + among the best
		{
			return 1;
		}
		elsif($somme_min_frame/10 == int($somme_min_frame/10)) # only frame - among the best
		{
			return -1;
		}
		else # there is a frame + and a frame - among the best frames
		{
			return 0;
		}
	}
	else # all frames have codon stops
	{
		return 0;
	}
}
################################################

sub clean_tax_names
{
	my ($lineage) = @_; 
# keep only tax names that has a correct latin name format
#  =~ /^[A-Z][a-z-]+\s*[a-z-]*$/

	unless($$lineage[-1] =~ /^[A-Z][a-z-]+\s[a-z-]+$/) # species level
	{
		$$lineage[-1] = '';
	}
	for(my $i = 0; $i < (scalar @$lineage -1); ++$i)
	{
		unless($$lineage[$i] =~ /^[A-Z][a-z-]+$/)
		{
			$$lineage[$i] = '';
		}
	}
}

######################################################


sub delete_partial_lines
{	my ($data, $sep, $filename) = @_;

	my @tmp;
	my @partial;

	for(my $i = 0; $i<scalar @$data; ++$i)
	{
		my $e = $$data[$i];
		$e =~ s/"//g;
		my @l = split($sep, $e);
		if (scalar @l < 80) # there are lines with less than 80 columns; Taxa do not exists in bold, file ends with xml lke lines 
		{
			########### patch to automatically repare a few lines in Archnida, where the line has been accidentally cut into two with a \n within a field
			#ASDEM271-21	CCDB-38232-G09	12795193		CCDB-38232-G09	Research Collection of M. Alex Smith			20	Arthropoda	63	Arachnida	298	Trombidiformes	506346	Demodicidae			506347	Demodex					Alex Smith							ZOO2700 students													Europe	YJMV|European- Polish|European- Polish|F|20|Nose|no|yes||	43.533	-80.226			334				Canada	Ontario			University of Guelph	5006947	http://www.boldsystems.org/pics/ASDEM/YJMV+1614267614.jpg	Dorsal		Alex Smith	2098	CreativeCommons - Attribution Non-Commercial
			#Share-Alike	University of Guelph	Sheri Hincks	14243487	COI-5P		---------------CATCAAATCTGCAATTTAACAAATGAATTAACTATGCGATGAATATTTTCCACTAACCACAAAGACATTGGAACAATATATTTCCTACTAAGAATATGATCTGGCATAATAGGACTAAGTTTAAGCATAATCATTCGTATAGAATTAAGACAACCTGGCAAAATTTTCCAATCAGACCACACCTATAATTGCATTGTCACATCCCACGCAATCGTAATAATCTTCTTTATAGTGATACCTATACTAATCGGAGGATTCGGCAACTGAATAACACCCATTATACTAATAACACAAGATATAGCATTCCCACGAATGAACAACCTAAGATTTTGACTACTACCCCCATCTATTAACCTAGCAATAATAGCATCACTAACAGACAAAGGAAGAGGAACTGGCTGAACATTCTACCCACCCCTTTCTCTAGCCCCATACCACCCTGGACACTCCGTAGATTTAATAATTTTTTCCTTACATGTAGCAGGAATTTCTTCAATTATCGGATCAATTAACTTCCTAGCCACCATTATAAATATAAAACACAAATCAATAGCCCCAGATCGCGTACCCCTATTCATTTGAGCAGTAGCCACAACTACAATCCTACTCCTCCTATCTTTACCCGTACTAGCAGGTGGGATTACCATAATCCTAACAGACCGAAACTTCAGAAGCTCTTTCTTCGACCCAGCAGGAGGAGGTGACCCCATTCTATTCCAACACCTATTT--								
			if($e =~ /Non-Commercial\s*$/)
			{
				$e =~ s/\s*$//;
				$e .= $$data[$i+1];
				++$i;
				@l = split($sep, $e);
				if (scalar @l < 80)
				{
					push(@partial, $e);
				}
				else
				{
					push(@tmp, $e);
				}
			}
			##############
			else
			{
				push(@partial, $e);
			}
		}
		else
		{
			push(@tmp, $e);
		}
	}
	
	if(scalar @partial)
	{
		print SH "\n$filename\n";
		print SH join("", @partial);
	}
	@$data = @tmp;
}


######################################################

sub read_data
{
	my ($file) = @_;
	
	#### read data from one file
	open(IN, $file) or die "Cannot open $file\n";
	my @data = <IN>;
	close IN;
	return @data;
}


###################################

sub get_field_indices
{
	my ($name_list, $data) = @_;

	my %name_hash;
	foreach my $name (@$name_list)
	{
		$name_hash{$name} = '';
	}
	
	my %fields_hash;
	my $line = $$data[0];
	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split("\t", $line);
	for(my $i = 0; $i < scalar @line; ++$i)
	{
		if(exists $name_hash{$line[$i]})
		{
			$fields_hash{$line[$i]} = $i;
		}
	}
	return %fields_hash;
}

###########################################

sub print_help
{

print '
usage: perl format_bold.pl [-options] -download_dir DOWNLOAD_DIR -outdir OUTDIR

 ARGUMENTS
   -download_dir           Name of the directory containing the tsv files downloaded from BOLD 
   -outdir                 Name of the otput directory
 OPTIONS
   -marker_list            Default:  "COI-5P COI-3P"; List of markers to be selected
   -check_name             [0/1]; Default 1
                               If one keeps only taxa with valid Latin name format
   -max_n                  Positive integer; Default:5)
                              Eliminates sequences with max_n or more consecutive Ns
   -min_length             Positive integer; Default:100
                              Minimum length of the cleaned sequence;
   -max_length             Positive integer; Default:2000
                              Maximum length of the cleaned sequence;
   -check_orientation      [0/1]; Default 0
                              If 1, checks the orientation of the sequences (Experimental)
   -blast_path             Path to the blast executables; Not necessarry if it is in the PATH
',  "\n";
  exit;
}

