use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkdb;
use Data::Dumper;

# INPUT
#	taxon_list: Text file with a list of taxa (one taxon per line) created from https://www.boldsystems.org/index.php/TaxBrowser_Home

# AIM
# Download all sequneces and lineages for all taxa on the taxon list from BOLD
#	if file exists already, skip the download
#	Check if the number of downloaded records corresponds to the extected one
#	if error, remove file and retry download $try_download times


#OUTPUT
# One tsv file for each taxa in the download subdirectory
# The taxa with failed downloads are listed at the end of the log file and their tsv files are deleted from the download subdirectory
# Columns are the following:
# processid	sampleid	recordID	catalognum	fieldnum	institution_storing	collection_code	bin_uri	phylum_taxID	phylum_name	class_taxID	class_name	order_taxID	order_name	family_taxID	family_name	subfamily_taxID	subfamily_name	genus_taxID	genus_name	species_taxID	species_name	subspecies_taxID	subspecies_name	identification_provided_by	identification_method	identification_reference	tax_note	voucher_status	tissue_type	collection_event_id	collectors	collectiondate_start	collectiondate_end	collectiontime	collection_note	site_code	sampling_protocol	lifestage	sex	reproduction	habitat	associated_specimens	associated_taxa	extrainfo	notes	lat	lon	coord_source	coord_accuracy	elev	depth	elev_accuracy	depth_accuracy	country	province_state	region	sector	exactsite	image_ids	image_urls	media_descriptors	captions	copyright_holders	copyright_years	copyright_licenses	copyright_institutions	photographers	sequenceID	markercode	genbank_accession	nucleotides	trace_ids	trace_names	trace_links	run_dates	sequencing_centers	directions	seq_primers	marker_codes

# Potential problems of downloading taken into account by the script:
	# Do not select for a marker when downloading
		# If a processID has at least one marker sequence, all sequences are downloaded, not just the marker, so a post download selection is necessary anyway
		# If taxon does not have COI, (or taxon name does not exists) a download of random sequneces starts, then stop with xml format lines at the end of the file. 
		# This pb is limited if there is no selection for the marker in the API command, and the taxon list is established from https://www.boldsystems.org/index.php/TaxBrowser_Home




###############################
###############################

my %params = (
	'taxon_list' => '',
	'outdir' => '',
	'try_download' => 3, # if pb download retry retry_download times
	'max_record_n' => 0, # if more than max_record_n records in json, cut up the taxon to subtaxa; 0 means Donot cut up taxa
	'taxonomy' => '' # necessary only if max_record_n != 0
);
modify_params_from_tags(\%params, \@ARGV);

my $taxon_list = $params{'taxon_list'};
my $try_download = $params{'try_download'};
my $outdir = $params{'outdir'};
my $taxonomy = $params{'taxonomy'};
my $max_record_n = $params{'max_record_n'};

$outdir = add_slash_to_dir($outdir);

###############################
###############################

#### define filenames
my $t = time;
my $date = get_date();
my $download_dir = $outdir.'files/';
my $log = $outdir.'download_bold.log';
unless(-e $download_dir)
{
	system 'mkdir -p '.$download_dir;
}
open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";


my %tax; #tax{taxid} = (parent_tax_id	rank	name_txt	tax_rank_index)
my %name_taxids; # $name_taxids{name}{$taxid} = '' # All names, including synonyms and homonyms
if($taxonomy)
{
	### Read taxonomy
	print "\n####\nRead taxonomy\n";
	print LOG "\n####\nRead taxonomy\n";
	read_taxonomy_to_tax_hash(\%tax, $taxonomy);
	my %taxid_names; # $taxid_names{taxid}{names} = '' # All names, including synonyms and homonyms
	my %name_taxids_par; #$name_taxids_par{$name}{taxid} = parent taxid # only scientific names
	read_taxonomy_names(\%taxid_names, \%name_taxids, \%name_taxids_par, $taxonomy);
	%taxid_names = (); # delete hash. It is produced automatically by the routine, but we do not need it.
	%name_taxids_par = (); # delete hash. It is produced automatically by the routine, but we do not need it.
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}


#### Read taxon names
print "\n####\nRead taxon names\n";
print LOG "\n####\nRead taxon names\n";
open(IN, $taxon_list) or die "Cannot open $taxon_list\n";
my @taxa = <IN>;
close IN;
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

#### Get exepected record count for each taxa (by downloading summary in json format)
# if expected record count in above $max_record_n => cut up taxon to subtaxa of the next major taxonomic level
print "\n####\nDownload summary files to get the expected number of records\n";
print LOG "\n####\nDownload summary files to get the expected number of reconrs\n";
my %updated_taxa; # %updated_taxa{taxon_name} = number of records based on json files # list of taxa is updated if with too many records, taxon have been cut up to sub-taxa
my %no_records; # $no_records{taxa} = ''; ŧaxa with no public record
my %fail_json; # $fail_json{taxa} = ''; error json download
download_json(\@taxa, $download_dir, $try_download, \%updated_taxa, \%no_records, \%fail_json);
print LOG "Runtime: ", time - $t, "s \n";
$t = time;


if($max_record_n) # nonZero, cut up taxa to subtaxa, that have more than $max_record_n records
{
	#### cut up taxa if too many records
	print "\n####\nCut up taxa to sub-taxa if more than $max_record_n records\n";
	print LOG "\n####\nCut up taxa to sub-taxa if more than $max_record_n records\n";
	my ($taxon, $max_record) = get_taxon_with_max_record_count(\%updated_taxa);
	while($max_record > $max_record_n)
	{
		# get lust of sub_taxa
		my @cut_up_taxon = cut_up_taxon($taxon, \%name_taxids, \%tax);

		# get expected record count for sub_taxa
		download_json(\@cut_up_taxon, $download_dir, $try_download, \%updated_taxa, \%no_records, \%fail_json);
		# eliminate taxa with too many records
		delete $updated_taxa{$taxon}; 
		# reinitialize the taxon with the maximum number of records
		($taxon, $max_record) = get_taxon_with_max_record_count(\%updated_taxa);
	}
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}

#### Download sequences by taxon
print "\n####\nDownload sequences taxon by taxon\n\n";
print LOG "\n####\nDownload sequences taxon by taxon\n";
my %fail; # $fail{taxon} = 1/2 1 if download started but did not finish correctly; 2 if connection loss; 3 empty file
@taxa = sort keys %updated_taxa; # download data only if public reports and less than $max_record_n records
for(my $i = 0; $i < $try_download; ++$i)
{
	print LOG "\n#### Download Iteration ", $i+1, "\nTaxon\tExpected Number of ProcessID\tDownloaded Number of ProcessID\n";
	print "\n#### Download Iteration ", $i+1, "\n";
	%fail = (); # empty the fail list at the begining of feach iteration.
	download_tsv(\@taxa, \%fail, $download_dir, \%updated_taxa);
	if(scalar keys %fail) # some downloads have failed
	{
		@taxa = sort keys %fail;
	}
	else # quit the loop
	{
		last;
	}
}

if(scalar keys %fail)
{
	print LOG "\n#####\nDownload failed for the following taxa:\n";
	print "\n#####\nDownload failed for the following taxa:\n";
	foreach my $taxon (sort keys %fail)
	{
		print "$taxon\t", join("\t", @{$fail{$taxon}}) ,"\n";
		print LOG "$taxon\t", join("\t", @{$fail{$taxon}}) ,"\n";
	}
}
else
{
	print LOG "\n#####\nSequneces have been downloaded sucessfully for all taxa\n";
	print "\n#####\nSequneces has been downloaded sucessfully for all taxa\n";
}
print LOG "Runtime: ", time - $t, "s \n";
$t = time;

close LOG;

exit;

###############################
sub cut_up_taxon
{
	my ($taxon, $name_taxids, $tax) = @_;
	
	my @cut_up_taxon;
	my @taxids;
	
	# get list of taxids that correspond to taxon name
	if(exists $$name_taxids{$taxon}) 
	{
		foreach my $taxid (keys %{$name_taxids{$taxon}})
		{
			push(@taxids, $taxid);
		}
	}
	
	if((scalar @taxids) == 0)
	{
		print "$taxon is not in the list of taxa in $taxonomy file\n";
		exit;
	}
	
	# get taxon list of the sub_taxa of the next major taxonomic rank
	foreach my $taxid (@taxids)
	{
		my $taxlevel = int($$tax{$taxid}[3] + 1); # get the taxlevel index for the next major taxonomic level
		foreach my $tid (keys %$tax) # go through all taxids if %tax
		{
			if($tax{$tid}[3] == $taxlevel) # correct taxlevel
			{
				my @lin = get_lineage_list($tid, \%tax, 0); # get complete lineages of taxids
				for(my $i=0; $i < scalar @lin; ++$i) # if input taxid is in the lineage of the taxon, put taxon to the subtaxon list
				{
					if($lin[$i] == $taxid)
					{
						push(@cut_up_taxon, $tax{$tid}[2]);
						last;
					}
				}
			}
		}
	}

	return @cut_up_taxon;
}
###############################
sub get_taxon_with_max_record_count
{
	my ($taxa) = @_;
	
	
	if(scalar keys %$taxa)
	{
		my @record_counts = sort {$a <=> $b} values %$taxa;
		my $max = $record_counts[-1];
		foreach my $tax (sort keys %$taxa)
		{
			if($$taxa{$tax} == $max)
			{
				return ($tax, $$taxa{$tax});
			}
		}
	}
	else
	{
		return ('NA', 0);
	}

}

###############################
sub download_json
{
	my ($taxa, $download_dir, $try_download, $updated_taxa, $no_records, $fail_json) = @_;
	

	foreach my $taxon (@$taxa)
	{
		### get filenames
		$taxon =~ s/\s*$//;
		$taxon =~ s/"//g;
		if($taxon eq 'taxon_name')
		{
			next;
		}
		my $tsv = get_file_name($taxon, $download_dir); # name of the download file
		my $json = $tsv;
		$json =~ s/\.tsv/.json/;
		
		my $expected_count;
		for(my $i = 0; $i < $try_download; ++$i)
		{
			### download stat file
			download_bold_stat($taxon, $json); # Download all sequences and metainfo for taxon
			$expected_count = get_count_json($json);
			if($expected_count eq 'NA') # pb with download
			{
				system 'rm '.$json;
			}
			else
			{
				last; # Download json OK
			}
		}

		if($expected_count eq 'NA') # pb with download even after $try_download times
		{
			$$fail_json{$taxon} = 'NA';
			print LOG "$taxon	NA	NA\n";
			print "$taxon	json download failed\n";
		}
		elsif($expected_count) # Public data exists
		{
			$$updated_taxa{$taxon} = $expected_count;
		}
		else # No public data for the taxon
		{
			$$no_records{$taxa} = '';
			print LOG "$taxon	No Public Data\n";
			print "$taxon	No Public Data\n";
		}
	}
}

###############################
sub download_tsv
{
	my ($taxa, $fail, $download_dir, $updated_taxa) = @_;
	
	foreach my $taxon (@$taxa)
	{
		### get filenames
		$taxon =~ s/\s*$//;
		$taxon =~ s/"//g;
		my $tsv = get_file_name($taxon, $download_dir); # name of the download file
		
		my $expected_count = $$updated_taxa{$taxon};

		download_bold_combined($taxon, $tsv); # Download all sequences and metainfo for taxon
		my $processID_count = get_count_tsv($tsv); # get the number of processIDs

		if( (($expected_count - $processID_count) / $expected_count) < 0.001) # less then 0.1% difference between expected count and real download
		{
			print LOG "$taxon	$expected_count	$processID_count\n";
			print "$taxon	OK\n";
		}
		else
		{
			@{$$fail{$taxon}} = ($expected_count, $processID_count);
			print "$taxon	FAILED\n";
			system 'rm '.$tsv;
		}
	}
}

###############################
sub download_and_check
{
	my ($taxa, $fail, $download_dir) = @_;
	
	foreach my $taxon (@$taxa)
	{
		### get filenames
		$taxon =~ s/\s*$//;
		$taxon =~ s/"//g;
		my $tsv = get_file_name($taxon, $download_dir); # name of the download file
		my $json = $tsv;
		$json =~ s/\.tsv/.json/;
		
		### download stat file
		download_bold_stat($taxon, $json); # Download all sequences and metainfo for taxon
		my $expected_count = get_count_json($json);
		
		if($expected_count eq 'NA') # pb in downloading stat file , try in the next round
		{
			@{$$fail{$taxon}} = ('NA', 'NA');
			print "$taxon	FAILED\n";
			system 'rm '.$json;
		}
		elsif($expected_count) # Public data exists
		{
			download_bold_combined($taxon, $tsv); # Download all sequences and metainfo for taxon
			my $processID_count = get_count_tsv($tsv); # get the number of processIDs
			
			if( (($expected_count - $processID_count) / $expected_count) < 0.001) # less then 0.1% difference between expected count and real download
			{
				print LOG "$taxon	$expected_count	$processID_count\n";
				print "$taxon	OK\n";
			}
			else
			{
				@{$$fail{$taxon}} = ($expected_count, $processID_count);
				print "$taxon	FAILED\n";
				system 'rm '.$tsv;
			}
		}
		else # No public data for the taxon
		{
			print LOG "$taxon	0	0\n";
			print "$taxon	No Public Data\n";
		}
	}
}

##################################################

sub get_count_tsv
{
	my ($tsv) = @_;
	
	my $cmd = 'cut -f1 '.$tsv.' | sort | uniq | wc -l';
	my $c = `$cmd`;
	$c =~ s/\s*$//;
	--$c; # do not count the title line
	return $c;
}

##################################################

sub get_count_json
{
	my ($json) = @_;
	# {"total_records":2139,"records_with_species_name":1421,"bins":{"count":184,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"countries":{"count":51,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"depositories":{"count":11,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"order":{"count":9,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"family":{"count":17,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"genus":{"count":42,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}},"species":{"count":107,"drill_down":{"entity":[{"name":null,"records":null},{"name":null,"records":null}]}}}
	
	open(IN, $json) or die "Cannot open $json\n";
	my $line = <IN>;
	close IN;
	
	if ($line =~ /^\{\"total_records\"\:([0-9]+)/)
	{
		return $1;
	}
	else
	{
		return 'NA';
	}
	
}


######################################################

sub download_bold_combined
{
	my ($taxon, $file) = @_;
	
	my $cmd = 'wget -q -nc -O '.$file.' "http://www.boldsystems.org/index.php/API_Public/combined?taxon='.$taxon.'&format=tsv"';
	system $cmd;
}

######################################################

sub download_bold_stat
{
	my ($taxon, $file) = @_;
	
	my $cmd = 'wget -q -nc -O '.$file.' "http://www.boldsystems.org/index.php/API_Public/stats?taxon='.$taxon.'"';
	system $cmd;
}


######################################################

sub get_file_name
{
	my ($taxon, $outdir) = @_;
	
	$taxon =~ s/\s*$//;

	my $file = $taxon;
	$file =~ s/[^a-z0-9_]/_/gi;
	$file = $outdir.$file.'.tsv';

	return $file;
}

###########################################

sub print_help
{

print '
usage: perl download_bold.pl [-options] -taxon_list TAXON_LIST -outdir OUTDIR

 ARGUMENTS
   -taxon_list             Input text file with one taxon per line
   -outdir                 Name of the otput directory
 OPTIONS
   -try_download           Integer; Default 3
                              Try to download files try_download times if some of the downloaded files are incomplete
   -max_record_n           If more than max_record_n records for a taxon, cut up the taxon to subtaxa;
                              Default: 0, means do not cut up the input taxa, not even if many records are in BOLD 
   -taxonomy               Input tsv with all taxids: tax_id,parent_tax_id,rank,name_txt,old_tax_id,taxlevel,synonyms; 
                              Not necessary if max_record_n is 0
',  "\n";
  exit;
}

