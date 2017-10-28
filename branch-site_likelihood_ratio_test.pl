#!/usr/bin/perl 

use strict;
use warnings;
use Data::Dumper;

my %alt_hash;
my %null_hash;

# alt_filename is a text file containing the ID and lnL value for a given sequence under the alternative model
# null_filename is a text file containing the ID and lnL value for the same sequence under the null model
# both files are created by script collect_lnL.pl which takes as its input the rst file produced by PAML

my ($alt_filename, $null_filename, $output_filename) = @ARGV;

# Open input files
unless (open (ALT, $alt_filename) ) {

	print STDERR "Can't open rst for $alt_filename\n";
	exit;

}

unless (open (NULL, $null_filename) ) {

	print STDERR "Can't open rst for $null_filename\n";
	exit;

}


unless (open (OUTPUT, ">$output_filename") ) {

	print STDERR "Can't open output file $output_filename\n";
	exit;

}


my @alt_lnL = <ALT>;

foreach my $line (@alt_lnL) {

	my @split_line = split(' ', $line);

	@split_line = grep { $_ ne '' } @split_line;

	if (scalar(@split_line == 2) ) {

		$alt_hash{$split_line[0]} = $split_line[1];

	}

}


my @null_lnL = <NULL>;

foreach my $line (@null_lnL) {

	my @split_line = split(' ', $line);

	@split_line = grep { $_ ne '' } @split_line;

	if (scalar(@split_line == 2) ) {

		$null_hash{$split_line[0]} = $split_line[1];

	}

}


#print Dumper(\%alt_hash);


foreach my $key (keys %alt_hash) {

	if (exists $null_hash{$key} ) {

		my $likelihood_ratio_test = 2 * (($alt_hash{$key}) - ($null_hash{$key}));

		if ($likelihood_ratio_test >= 5.41) {

			print OUTPUT "$key\t$likelihood_ratio_test\t1_percent\n";

		} elsif ($likelihood_ratio_test >= 2.71) {

			print OUTPUT "$key\t$likelihood_ratio_test\t5_percent\n";

		}

	}

}

		








