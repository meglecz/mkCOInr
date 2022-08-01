use warnings;
use strict;
use FindBin qw( $RealBin );
use lib "$RealBin";
use Data::Dumper;

# run the same script with different parameters

my %params = (
'param_file' => '/home/meglecz/mkCOInr/no_git/benchmark_select_region/run_select_region.tsv',
'script_name' => 'run_select_region.pl'
);

modify_params_from_tags(\%params, \@ARGV);

my $param_file = $params{param_file};
my $script_name = $params{script_name};

open(IN, $param_file) or die  "Cannot open $param_file\n";
my $title = <IN>;
chomp $title;
$title =~ s/"//g;
my @title = split("\t", $title);

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/"//g;
	unless($line =~ /^#/)
	{
		my @line = split("\t", $line);
		my $cmd = 'perl '.$script_name;
		for(my $i = 0; $i < scalar @title; ++$i)
		{
			$cmd .= ' -'.$title[$i].' '.$line[$i];
		}
		print $cmd, "\n\n";
		system $cmd;
	}
}
close IN;


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
