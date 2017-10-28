#!/usr/bin/perl

use strict;
use warnings;
use Bio::SimpleAlign;
use Bio::AlignIO;

# Prints a multifasta file's name if contained sequences have an
# overall_average_percentage_identity >= a specified minimum.
# Multifasta alignment and minimum overall_percentage_identity (BioPerl def)
my ($file, $min) = @ARGV;

chomp ($file);

my $alignment = Bio::AlignIO->new(-file => "$file", -format => "phylip");

while (my $aln = $alignment->next_aln()) {

	my $identity = $aln->overall_percentage_identity();

	if ($identity >= $min) {

		print "$file\n";

	}

}
