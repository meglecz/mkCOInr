use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;


# INPUT
# tsv file with the following columns:
#	seqID	taxID	sequence


# ALGO
# foreach taxid, take all sequneces => is they are a substring of another sequence => delete
# If more than 10000 sequences with a taxid => cluster first sequences with 100 %id, than use the substring serch for each cluster

# OUTPUT
# tsv file with the following columns:
#	seqID	taxID	sequence


###############################
###############################
my %params = (
'tsv' => '', # seqID	taxID	sequence
'outdir' => '',
'out' => 'dereplicated_sequences.tsv',
'vsearch_path' => ''
);

modify_params_from_tags(\%params, \@ARGV);

my $tsv = $params{tsv};
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
my $log = $outdir.'dereplicate.log'; 
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
my %tseq; # $tseq{taxid}{sequence} = seqid # if identical sequences just keep one seqid keep all input sequences
my ($seqn_all, $seqn_uniq) = read_seq_tsv(\%tseq, $tsv);

$stat{'1. Number of input taxids: '} = scalar keys %tseq;
$stat{'2. Number of input sequences: '} = $seqn_all;
$stat{'3. Number of input unique sequences: '} = $seqn_uniq;
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
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
#### 

####
print "\n####\nPrint output tsv file\n";
print LOG "\n####\nPrint output tsv file\n";
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "seqID	taxID	sequence\n";
foreach my $taxid (sort keys %tseq)
{
	++$stat{'4. Number of output taxids: '};
	foreach my $seq (sort keys %{$tseq{$taxid}})
	{
		++$stat{'5. Number of output sequences: '};
		print OUT "$tseq{$taxid}{$seq}	$taxid	$seq\n";
	}
}
close OUT;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

#### delete tmp_dir
system 'rm -rf '.$tmpdir;

print LOG print_stat(\%stat, $t0);
close LOG;

exit;

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

sub read_seq_tsv
{
	my ($hash, $file) = @_;
	# $hash{taxid}{sequence} = seqid # if identical sequnces just keep on seqid
	#seqID	TaxID	sequence
	
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>;
	my $all = 0;
	my $uniq = 0;
	my %seqids; # $seqids{sid} = ''
	my %pb_seqids; 
	while(my $line = <IN>)
	{
		++$all;
		$line =~ s/\s*$//;
		$line=~ s/"//g;
		my @line = split("\t", $line);

		unless(exists $$hash{$line[1]}{$line[2]})
		{
			$$hash{$line[1]}{$line[2]} = $line[0];
			++$uniq; # count unique sequneces within taxon
			if(exists $seqids{$line[0]}) # check if seqid is unique
			{
				$pb_seqids{$line[0]} = '';
			}
			$seqids{$line[0]} = '';
		}
	}
	close IN;
	
	if(scalar keys %pb_seqids) # print error and stop if seqids are not uniques
	{
		print "ERROR: The following sequenceIDs are not unique in the input file:\n";
		print join("\n", sort keys %pb_seqids);
		print LOG "ERROR: The following sequenceIDs are not unique in the input file:\n";
		print LOG join("\n", sort keys %pb_seqids);
		exit;
	}
	%seqids = ();
	
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

############################################

sub print_help
{

print '
usage: perl dereplicate.pl [-options] -tsv INPUT_FILE -outdir OUTDIR

 ARGUMENTS
   -tsv                    Input tsv file with seqID,taxID,sequence
   -outdir                 Name of the otput directory
 OPTIONS
   -out                    Name of the output tsv file
   -vsearch_path           Path to the vsearch executables; Not necessarry if it is in the PATH
',  "\n";
  exit;
}
