#!/usr/local/bin

use strict;
use warnings;

my ($res_one, $res_two) = @ARGV;

unless (open (RESONE, $res_one) ) {
	die "Can't open $res_one\n";
}

unless (open (RESTWO, $res_two) ) {
        die "Can't open $res_two\n";
}

my %res_one;
my %res_two;
my %res_combined;

# Collect all ensembl geneIDs from results file one
while (my $line = <RESONE>) {

	my @split_line = split ('_', $res_one);

	my $res_one{$split_line[0]} = 0;

}

# Collect all ensembl geneIDs from results file two
while (my $line = <RESTWO>) {

        my @split_line = split ('_', $res_two);

        my $res_two{$split_line[0]} = 0;

}

# Check keys from results one against results two, collect overlap
foreach (my $key (keys %res_one) {

	if (exists $res_two{$key}) {

		$res_combined{$key} = 0;

	}

}

# Check keys from results two against results one, collect overlap
foreach (my $key (keys %res_two) {

        if (exists $res_one{$key}) {

                $res_combined{$key} = 0;

        }

}


# Print the overlapping IDs
foreach (my $key (keys %res_combined) ) {

	print "$key\n";

}
