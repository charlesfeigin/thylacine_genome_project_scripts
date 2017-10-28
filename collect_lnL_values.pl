#!/usr/bin/perl 

use strict;
use warnings;

# Collect each likelihood ratio into a tab-delimited file 
# with the first field being the phylip basename and the
# second being the lnL value. These will be concatenated
# together via the commandline afterward. The second script
# will them take both files as inputs, make hashes with the
# phylip basename as the key and the lnL as the value.
# Then, for each null, if the key exists in the alternative value
# feed them both to the likelihood ratio test.


my $pause;

my $rst_filename = "rst";

my $lnL_line;

my ($phylip_filename, $model) = @ARGV;

my @phylip_filename = split (/\./, $phylip_filename);

my $basename = $phylip_filename[0];

my $output_filename = "$basename" . "_lnL_$model" . ".txt";

unless (open (RST, $rst_filename) ) {

	print STDERR "Can't open rst for $phylip_filename\n";
	exit;

}

unless (open (OUTPUT, ">$output_filename") ) {

	print STDERR "Can't open output file $output_filename\n";
	exit;

}

while (my $line = <RST>) {

	if ($line =~ /^lnL/) {

		$lnL_line = $line;

		last;

	}

}


my @split_line = split (' ', $lnL_line);

my $lnL_value = $split_line[2];

print OUTPUT "$basename\t$lnL_value\n";

close RST;
close OUTPUT;

exit;



























		


