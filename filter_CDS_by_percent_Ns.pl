#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;


my $USAGE = 'This script filters a multifasta file to remove sequences with too many N characters.

USAGE: perl Filter_cds_by_percent_Ns.pl [cds_filename] [max_percent_Ns] [find_suffix] [replace_suffix]

cds_filename -	Name of multifasta containing Ensembl Biomart CDS sequences. 
max_percdent_Ns -	Maximum percent of "N" characters allowed in CDS sequence
find_suffix -	suffix in input file to replace when creating output file
replace_suffix -	new suffix for output file

';

#--------------------------Parse Arguments---------------------------#

# Check for arguments. If more or less than 4 are given print usage message and exit
if (scalar (@ARGV) != 4) {
		print $USAGE;
		exit;
}

# Collect arguments
my ($cds_filename, $max_N_percent, $find, $replace) = @ARGV;

# Remove metacharacters from find_suffix and create output filename
$find = quotemeta $find;
my $output_filename = $cds_filename;
$output_filename =~ s/$find/$replace/g;

# Open input file for reading, open output file for writing
my $multifasta = Bio::SeqIO->new(-file => "$cds_filename", -format => 'fasta');
my $output = Bio::SeqIO->new(-file => ">$output_filename", -format => 'fasta');


# Main loop. Read through each fasta and filter based on the number of Ns.
while (my $sequence_object = $multifasta->next_seq) {

	# Initial N count
	my $N_count = 0;

	# Get sequence
	my $sequence = $sequence_object->seq();

	my $sequence_length = length($sequence);

	my @split_seq = split('', $sequence);

	# Count Ns
	foreach my $base (@split_seq) {

		if ($base =~ /n/ig) {

			$N_count++;

		}

	}
	
	# Check N number, output or toss.
	my $N_percent = 100*($N_count/$sequence_length);

	unless ($N_percent > $max_N_percent) {

		$output->write_seq($sequence_object);

	}

}

exit;
	

		
