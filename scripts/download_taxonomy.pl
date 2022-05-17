use warnings;
use strict;
use mkdb;
use Data::Dumper;

# Downloads ncbi taxonomy dump files from https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/ 
# the decompressed taxonomy dmp files are in the ncbi_tax subdirectory

# prepares a taxonomy file (taxonomy_yyyy-mm-dd-hh-mm-ss.tsv) file with the following tab separated columns
# - tax_id
# - parent_tax_id
# - rank
# - name_txt
# - old_tax_id (taxIDs merged to another taxID, so not valid anymore, but occasionally still present in ncbi flatfiles)
# - taxlevel (8:species, 7:genus, 6:family, 5:order, 4:class, 3:phylum, 2:kingdom, 1:superkingdom; x.5 for levels between major taxonomic levels e.g. 5.5 for infra-orger or superfamily)
# - list of synonymes separated by ;



###############################
###############################
my %params = 
(
	'outdir' => '',
	'skip_download' => 0
);
modify_params_from_tags(\%params, \@ARGV);


my $skip_download = $params{'skip_download'};
my $outdir = $params{'outdir'};

$outdir = add_slash_to_dir($outdir);
###############################
###############################

my $t = time;
my $date = get_date();
#### define filenames
my $ncbitax_dir = $outdir.'ncbi_tax/';
unless(-e $ncbitax_dir)
{
	system 'mkdir -p '.$ncbitax_dir;
}
my $log = $outdir.'download_taxonomy.log';
open(LOG, '>', $log) or die "Cannot open $log\n";
my @parameters = print_params_hash_to_log(\%params, $date);
print LOG join("\n", @parameters), "\n";

#### Download ncbi taxonomy dmp file
unless($skip_download)
{
	print "\n####\nDownload ncbi taxonomy dmp files\n";
	print LOG "\n####\nDownload ncbi taxonomy dmp files\n";
	# download new_taxdump.tar.gz to $outdir
	my $cmd = 'wget -O '.$outdir.'new_taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz';
	system $cmd;
	# extract file
	$cmd = 'tar -zxvf '.$outdir.'new_taxdump.tar.gz -C '.$ncbitax_dir;
	system $cmd;
	# delete zipped file
	system 'rm '.$outdir.'new_taxdump.tar.gz';
	print LOG "Runtime: ", time - $t, "s \n";
	$t = time;
}
####

####Make taxonomy file from ncbi tax dmp files
print "\n####\nMake taxonomy file from ncbi tax dmp files\n";
print LOG "\n####\nMake taxonomy file from ncbi tax dmp files\n";
my $taxonomy = make_taxonomy_with_rank_levels_synonyms($ncbitax_dir, $outdir, '');
print LOG "Runtime: ", time - $t, "s \n";
$t = time;
####

close LOG;
exit;

###########################################

sub print_help
{

print '
usage: perl download_taxonomy.pl [-options] -outdir OUTDIR

 ARGUMENTS
   -outdir                 Name of the otput directory
 OPTIONS
   -skip_download          [0/1]; Default 0
                               If 1 skips download (dmp files has already been downloaded)
',  "\n";
  exit;
}

