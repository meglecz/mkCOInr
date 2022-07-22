use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;


# INPUT
# 2 tsv files with the following columns:
#	seqID	taxID	sequence
# Both input files have already been dereplicated separatelly


# ALGO
# take the list of taxids common to the 2 input files
# dereplicate only sequneces of these taxids
# 	foreach taxid, take all sequneces => is they are a substring of another sequence => delete
# 	If more than 10000 sequences with a taxid => cluster first sequences with 100 %id, than use the substring serch for each cluster

# OUTPUT
# tsv file with the following columns:
#	seqID	taxID	sequence


###############################
###############################
my %params = (
'tsv1' => '', # seqID	taxID	sequence
'tsv2' => '', # seqID	taxID	sequence
'outdir' => '',
'out' => '',
'vsearch_path' => ''
);

modify_params_from_tags(\%params, \@ARGV);

my $tsv1 = $params{tsv1};
my $tsv2 = $params{tsv2};
my $out = $params{out};
my $outdir = $params{outdir};
my $vsearch_path = $params{vsearch_path};

$outdir = add_slash_to_dir($outdir);
$vsearch_path = add_slash_to_dir($vsearch_path);
###############################
###############################

#### define filenames and variables
my $t = time;
my $t0 = time;
my $date = get_date();
my $log = $outdir.'pool_and_dereplicate.log'; 
$out = $outdir.$out;
my %stat;

my $tmpdir = $outdir.'tmp_'.$date.'/';
unless(-e $tmpdir)
{
	my $cmd = 'mkdir -p '.$tmpdir;
	system $cmd;
}
open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";
#### 


####
print "\n####\nReading sequences\n";
print LOG "\n####\nReading sequences\n";

my %seqids; # $seqids{$sid} = ''; # to check if all seqids are uniques
my %tseq1; # $tseq{taxid}{sequence} = seqid # if identical sequences just keep one seqid keep all input sequences
my ($seqn_all, $seqn_uniq) = read_seq_tsv2(\%tseq1, $tsv1, \%seqids);
$stat{'1.1 Number of input taxids in tsv1: '} = scalar keys %tseq1;
$stat{'1.2 Number of input sequences in tsv1: '} = $seqn_all;
$stat{'1.3 Number of input unique sequences in tsv1: '} = $seqn_uniq;

my %tseq2; # $tseq{taxid}{sequence} = seqid # if identical sequences just keep one seqid keep all input sequences
($seqn_all, $seqn_uniq) = read_seq_tsv2(\%tseq2, $tsv2, \%seqids);
$stat{'2.1 Number of input taxids in tsv2: '} = scalar keys %tseq2;
$stat{'2.2 Number of input sequences in tsv2: '} = $seqn_all;
$stat{'2.3 Number of input unique sequences in tsv2: '} = $seqn_uniq;

print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

####
print "\n####\nSelect taxids common to the input files\n";
print LOG "\n####\nSelect taxids common to the input files\n";
my %tseq;
my $shared_sequences_all = 0;
my $shared_sequences_uniq = 0;
foreach my $taxid (keys %tseq1)
{
	if(exists $tseq2{$taxid})
	{
		($seqn_all, $seqn_uniq) = pool_taxids(\%tseq, $tseq1{$taxid}, $tseq2{$taxid}, $taxid);
		$shared_sequences_all += $seqn_all;
		$shared_sequences_uniq += $seqn_uniq;
		delete $tseq1{$taxid}; # delete taxids from the hashes specific to each of the input files
		delete $tseq2{$taxid};
	}
}
$stat{'3.1 Number of taxids shared berween tsv1 and tsv2 : '} = scalar keys %tseq;
$stat{'3.2 Number of sequences for shared taxids : '} = $shared_sequences_all;
$stat{'3.3 Number of unique sequences (after eliminating exact match) for shared taxids : '} = $shared_sequences_uniq;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####


#### Dereplication
print "\n####\nDereplication\n";
print LOG "\n####\nDereplication\n";
foreach my $taxid (sort {$a <=>$b} keys %tseq)
{
	if((scalar keys %{$tseq{$taxid}}) > 10000) # several sequences => start by clustering 100% id
	{
		derep_cluster(\%tseq, $taxid, $tmpdir);
	}
	elsif((scalar keys %{$tseq{$taxid}}) > 2) # not too many sequences => use motif search
	{
		my @seqs = keys %{$tseq{$taxid}};
		derep_string_match3(\@seqs, \%tseq, $taxid); # eliminate sequences from %tseq if substring of another sequence with the same taxid
	}

}
print LOG time - $t, " s\n";
$t = time;
#### 

####
# taxids shared between input files
print "\n####\nPrint output tsv file\n";
print LOG "\n####\nPrint output tsv file\n";
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "seqID	taxID	sequence\n";
foreach my $taxid (sort keys %tseq)
{
	++$stat{'4.1 Number of output taxids shared by tsv1 and tsv2 : '};
	foreach my $seq (sort keys %{$tseq{$taxid}})
	{
		++$stat{'4.2 Number of dereplicated sequences for shared taxids: '};
		print OUT "$tseq{$taxid}{$seq}	$taxid	$seq\n";
	}
}

# taxids unique to tsv1
foreach my $taxid (sort keys %tseq1)
{
	++$stat{'5.1 Number of taxids unique to tsv1: '};
	foreach my $seq (sort keys %{$tseq1{$taxid}})
	{
		++$stat{'5.2 Number of dereplicated sequences from taxids unique to tsv1: '};
		print OUT "$tseq1{$taxid}{$seq}	$taxid	$seq\n";
	}
}

# taxids unique to tsv2
foreach my $taxid (sort keys %tseq2)
{
	++$stat{'6.1 Number of taxids unique to tsv2: '};
	foreach my $seq (sort keys %{$tseq2{$taxid}})
	{
		++$stat{'6.2 Number of dereplicated sequences from taxids unique to tsv2: '};
		print OUT "$tseq2{$taxid}{$seq}	$taxid	$seq\n";
	}
}
close OUT;

$stat{'7.1 Total number of output taxids: '} =  $stat{'4.1 Number of output taxids shared by tsv1 and tsv2 : '} + $stat{'5.1 Number of taxids unique to tsv1: '} + $stat{'6.1 Number of taxids unique to tsv2: '};
$stat{'7.2 Total number of dereplicated sequences: '} = $stat{'4.2 Number of dereplicated sequences for shared taxids: '} + $stat{'5.2 Number of dereplicated sequences from taxids unique to tsv1: '} + $stat{'6.2 Number of dereplicated sequences from taxids unique to tsv2: '};

print LOG time - $t, " s\n";
$t = time;
####

#### delete tmp_dir
if(-e $tmpdir)
{
	system 'rm -r '.$tmpdir;
}

print LOG print_stat(\%stat, $t0);
close LOG;

exit;

###################################

sub pool_taxids
{
	my ($tseq, $tseq1_taxid, $tseq2_taxid, $taxid) = @_;
	#$tseq{taxid}{sequence} = seqid
	#$tseq1_taxid{sequence} = seqid
	
	my $all = 0; # count all sequences
	my $uniq = 0; # count unique sequences
	foreach my $seq (keys %$tseq1_taxid)
	{
		++$all;
		unless(exists $$tseq{$taxid}{$seq})
		{
			++$uniq;
			$$tseq{$taxid}{$seq} = $$tseq1_taxid{$seq};
		}
	}
	
	foreach my $seq (keys %$tseq2_taxid)
	{
		++$all;
		unless(exists $$tseq{$taxid}{$seq})
		{
			++$uniq;
			$$tseq{$taxid}{$seq} = $$tseq2_taxid{$seq};
		}
	}
	return($all, $uniq);
}

###################################

sub derep_cluster
{
	my ($tseq, $taxid, $tmpdir) = @_;
	
	# eliminate previous cluster dir and make a new
	my $clusterdir = $tmpdir.'cluster/';
	if(-e $clusterdir)
	{
		my $cmd = 'rm -r '.$clusterdir;
		system $cmd;
	}
	my $cmd = 'mkdir -p '.$clusterdir;
	system $cmd;

	my $fastmp = $tmpdir.'tmp.fas';
	make_fasta($$tseq{$taxid}, $fastmp);
	my $cluster = $vsearch_path.'vsearch --cluster_fast '.$fastmp.' --clusters '.$clusterdir.$taxid.'_'.' --id 1 --quiet'; # write output clusters to taxid_xxx files
	system $cluster;

	my $motif = '^'.$taxid.'_';
	my @files = get_file_list_from_folder($clusterdir, $motif);
	
	foreach my $file (@files) # derep by_string each cluster separatelly 
	{
		$file = $clusterdir.$file;
#print $file, "\n";
		open(IN, $file);
		my $i = -1;
		my @seqs;
		while(my $line = <IN>) # read sequences to @seqs (ids are not kept)
		{
			$line =~ s/\s*$//;
			if($line =~ /^>/)
			{
				++$i;
				$seqs[$i] = '';
			}
			else
			{
				$seqs[$i] .= uc $line; # uc is important
			}
		}
#		print join("\n", @seqs);
		derep_string_match3(\@seqs, $tseq, $taxid);
	}
	
	$cmd = 'rm -r '.$clusterdir;
	system $cmd;
}
###################################

sub derep_string_match
{
	my ($tseq, $taxid) = @_;
# $tseq{taxid}{sequence} = seqid # if identical sequnces just keep on seqid

	my @seq = keys %{$$tseq{$taxid}};
	my %delete;
	for(my $i = 0; $i< scalar @seq; ++$i)
	{
		for(my $j = $i+1; $j< scalar @seq; ++$j)
		{
			if($seq[$j] =~ /$seq[$i]/)
			{
				$delete{$seq[$i]} = '';
				last;
			}
			if($seq[$i] =~ /$seq[$j]/)
			{
				$delete{$seq[$j]} = '';
			}
		}
	}
	foreach my $seq (keys %delete)
	{
		delete $$tseq{$taxid}{$seq};
	}
}


###################################

sub derep_string_match2
{
	my ($seqs, $tseq, $taxid) = @_;
# $tseq{taxid}{sequence} = seqid # if identical sequences just keep on seqid

	my %delete;
#	print scalar keys %{$$tseq{$taxid}},"\n";
	for(my $i = 0; $i< scalar @$seqs; ++$i)
	{
		for(my $j = $i+1; $j< scalar @$seqs; ++$j)
		{
			if($$seqs[$j] =~ /$$seqs[$i]/)
			{
#				print "\n$$seqs[$j]\n$$seqs[$i]\n\n";
				++$delete{$$seqs[$i]};
				last;
			}
			if($$seqs[$i] =~ /$$seqs[$j]/)
			{
#							print "\n$$seqs[$j]\n$$seqs[$i]\n\n";
				++$delete{$$seqs[$j]};
			}
		}
	}
	foreach my $seq (keys %delete)
	{
		delete $$tseq{$taxid}{$seq};
	}
#	print scalar keys %{$$tseq{$taxid}},"\n";
#	print Dumper(\%delete);
}


###################################

sub derep_string_match3
{
	my ($seqs, $tseq, $taxid) = @_;
# $tseq{taxid}{sequence} = seqid # if identical sequences just keep on seqid

	my %l; #$l{length}{seq};
	foreach my $seq (@$seqs)
	{
		my $len = length $seq;
		$l{$len}{$seq} = '';
	}
	
	my @len_list = reverse sort {$a <=> $b} keys %l;
	my $max_len = shift @len_list;
	my %keep; # $keep{seq} ='';
	%keep = %{$l{$max_len}}; # Put the longest sequences to the %keep, since they cannot be a substring of another

	foreach my $seqlen (@len_list)# go though all other length in decreasing order
	{
		my %delete;
		foreach my $seq (keys %{$l{$seqlen}}) # check all sequences of a given length
		{
			foreach my $keep_seq (keys %keep) # go through all keep sequences (longer than $seq)
			{
				if($keep_seq =~ /$seq/) # put it on delete hash if it is a subsequence
				{
					$delete{$seq} = '';
					last;
				}
			}
		}
		# add remaining sequences to keep; 
		#it is better to do this step separatelly, to avoid adding
		#sequneces to a %keep that has the same length as the sequences to be tested
		foreach my $seq (keys %{$l{$seqlen}}) # go through all sequences of length $seqlen
		{
			unless(exists $delete{$seq}) # put it into %keep, if not on a delete list
			{
				$keep{$seq}= '';
			}
		}
	}
	
	# delete sequences that are not on the keep list
	foreach my $seq (@$seqs)
	{
		unless(exists $keep{$seq})
		{
			delete $$tseq{$taxid}{$seq};
		}
	}
}

###################################

sub read_seq_tsv2
{
	my ($hash, $file, $seqids) = @_;
	# $hash{taxid}{sequence} = seqid # if identical sequnces just keep on seqid
	#seqID	TaxID	sequence
	#%seqids; # $seqids{sid} = ''
	
	
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>;
	my $all = 0;
	my $uniq = 0;
	
	my %pb_seqids; 
	while(my $line = <IN>)
	{
		++$all;
		$line =~ s/\s*$//;
		$line =~ s/"//g;
		my @line = split("\t", $line);

		unless(exists $$hash{$line[1]}{$line[2]})
		{
			$$hash{$line[1]}{$line[2]} = $line[0];
			++$uniq; # count unique sequneces within taxon
			if(exists $$seqids{$line[0]}) # check if seqid is unique
			{
				$pb_seqids{$line[0]} = '';
			}
			$$seqids{$line[0]} = '';
		}
	}
	close IN;
	
	if(scalar keys %pb_seqids) # print error and stop if seqids are not uniques
	{
		print "ERROR: The following sequenceIDs are not unique:\n";
		print join("\n", sort keys %pb_seqids);
		print LOG "ERROR: The following sequenceIDs are not unique:\n";
		print LOG join("\n", sort keys %pb_seqids);
		exit;
	}
	
	return ($all, $uniq);
}

###################################

sub make_fasta
{
	my ($hash, $file) = @_;
	
	open(FAS, '>', $file) or die "Cannot open $file\n";
	foreach my $seq (keys %{$hash})
	{
		print FAS ">$$hash{$seq}\n$seq\n";
	}
	close FAS;

}

###################################
sub print_help
{

print '
usage: perl pool_and_dereplicate.pl [-options] -tsv1 INPUT_FILE -tsv2 INPUT_FILE -outdir OUTDIR -out OUTPUT_FILE

 ARGUMENTS
   -tsv1                    Input tsv file with seqID,taxID,sequence
   -tsv2                    Input tsv file with seqID,taxID,sequence
   -outdir                  Name of the otput directory
   -out                     Name of the output tsv file
 OPTIONS
   -vsearch_path            Path to the vsearch executables; Not necessarry if it is in the PATH
',  "\n";
  exit;
}
