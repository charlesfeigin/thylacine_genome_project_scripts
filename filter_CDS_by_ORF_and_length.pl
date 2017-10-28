#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;


#--------------Define starting variables------------------#
my $USAGE = 'This script filters a multifasta file pick the best CDS for each gene. Before running
this script, the multifasta must be sorted. The first sort criterion is the ID, which must be
in the format >ENSG00000002746|ENST00000453890|ENSP00000407774 where the first field is the 
ensembl gene ID, the second is ensembl transcript ID, third is ensembl protein ID. The second sorting
criterion must be sequence length in descending order. 

USAGE: perl Filter_CDS_by_ORF_and_length_draft.pl [input_filename] [find_suffix] [replace_suffix]

input_filename -	Name of multifasta containing Ensembl Biomart CDS sequences filtered by N percent
			and sorted by header (format: ensembl gene ID, transcript ID, protein ID,  
find_suffix -	suffix in input file to replace when creating output file
replace_suffix -	new suffix for output file

';


#--------------Collect arguments------------------#

# Check for arguments. If more or less than 3 are given print usage message and exit
if (scalar (@ARGV) != 3) {
		print $USAGE;
		exit;
}


# input multifasta filename
my ($input_filename, $find, $replace) = @ARGV;

# Remove metacharacters from find_suffix and create output filename
$find = quotemeta $find;
my $file_out = $input_filename;
$file_out =~ s/$find/$replace/g;

my $width = 32766;

# open multifasta for reading
my $multifasta = Bio::SeqIO->new(-file => "$input_filename", -format => 'fasta');
my $filtered = Bio::SeqIO->new(-file => ">$file_out", -format => 'fasta', -width => "$width");


# Initialize variables for loops
my $current_gene;
my $current_longest_ORF = "0";
my $current_longest_INC = "0";
my $current_longest_CDS = "0";


# Loop through sequences
while (my $seq_obj = $multifasta->next_seq) {

	# Get the sequence's ID
	my $id = $seq_obj->id();

	# Unless the $current_gene already contains a value, set it to the first
	# segment of the $id, corresponding to the gene ID (leaving out transcript
	# and protein. This should only be triggered on the first loop iteration
	unless ($current_gene) { 
		my @id = split(/\|/, $id);
		$current_gene = $id[0]; 
	}

	# Check if the $id matches the current gene. It always will on the first iteration,
	# and when a given gene has multiple CDS/transcripts
	if ($id =~ /^$current_gene/) {
		
		#print "$current_gene\n";
		# When a new gene ID is found and a previous entry closed, this
		# brings it back into the loop in proper order
		CHECKS: 
		
		# Get the sequence itself, chomp it and make sure its all uppercase
		my $seq = $seq_obj->seq();
		#chomp($seq);
		$seq = uc $seq;

		#-----------Collect some info on the sequence for filtering-----------#

		# Collect, start codon, stop codon and remainder from division by 3
		# In order to test if the CDS is a complete ORF	
		my $first_codon = substr ($seq, 0, 3);
		my $last_codon = substr ($seq, -3);
		my $div_remainder = length($seq) % 3;
		my $trunc_check = 'false';

		#print "$current_gene\t$first_codon\t$last_codon\t$div_remainder\t$trunc_check\n$seq\n\n";

		# Check for early stop codons by spliting it at every 3 letters (codons),
		# Removing the first and last codons (indices 1 and 2)...
		my @seq = ($seq =~ m/.../g);
		shift @seq;
		pop @seq;

		# Then seeing if any in-frame codon is a stop
		foreach my $codon (@seq) {

			# If it is...
			if (($codon =~ /TAG/ig) or ($codon =~ /TAA/ig) or ($codon =~ /TGA/ig)) {
				# Then take note
				$trunc_check = 'true';
			}
		}


		#---------------------------------------------------------------------#
		
		# Now check to see if the sequence constitutes a single complete ORF
		if ( ( ($last_codon =~ /TAG/ig) or ($last_codon =~ /TAA/ig) or ($last_codon =~ /TGA/ig) ) and	    
		     ($div_remainder == 0) and 
		     ($trunc_check eq 'false') and 
		     ($first_codon =~ /ATG/ig) ) {
#			print "ORF\t$id\n";
			
			# If it starts with a start, ends with a stop, is divisible by 3 and
			# has no premature stops, then its considered an ORF.
			# So, now check to see if it is longer than the previous longest
			# ORF for the same gene. If always will be the fist iteration.
			# If it is...
			if ( length($seq) > length($current_longest_ORF) ) {

				# Then replace the previous longest ORF with the current sequence
				$current_longest_ORF = $seq;
			
			}
		
		# If on the other hand it doesn't constitute a complete ORF, it may still
		# be useful. If there isn't a complete ORF, but the CDS begins with a start
		# codon, then it is likely to be an incomplete ORF that is still potentially
		# easy to translate. This is slightly preferable to a long transcript with
		# an unknown reading frame. 
		} elsif ($first_codon =~ /ATG/ig) {
#			print "INC\t$id\n";
			#print "$seq\n";
			# See if this sequence is longer than the current longest
			# incomplete ORF. If it is..
			if (length($seq) > length($current_longest_INC)) {

				# Then replace the previous longest CDS with the current
			 	$current_longest_INC = $seq;

			}
			

		# If there are no complete or incomplete ORFs for this gene, we'll take 
		# the longest transcript instead. If we haven't found an ORF yet, then 
		# $current_longest_ORF should still equal "0". Once we've found an ORF
		# This loop will never be activated again for this gene
		} else {
#			print "DERP\t$id\n";
			# See if this sequence is longer than the current longest 
			# non-ORF sequence. If it is...
			if (length($seq) > length($current_longest_CDS)) {

				# Then replace the previous longest CDS with the current
			 	$current_longest_CDS = $seq;

			}

		}

	# If on the other hand the $id doesn't correspond to the $current_gene 
	# That means its time to determine which sequence we'll be going with, 
	# either the longest ORF if any ORFs were found, or the longest transcript
	# seen.
	} else {
		
		# See if any ORFs have been found. If they have...
		if ($current_longest_ORF ne "0") {
	
			# Make a new sequence object with the current gene as the id
			# and the longest ORF as the sequence
			my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_ORF");
		
			$filtered->write_seq($best_cds);			
	
		# Else, if no complete ORFs were found, but an incomplete one was
		} elsif ($current_longest_INC ne "0") {
			#print "$current_gene\n$current_longest_INC\n";

			# Make a sequence object with the longest observed incomplete transcript 
			# and write it to the output file for that gene.
			my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_INC");

			$filtered->write_seq($best_cds);	

		} elsif ($current_longest_CDS ne "0") {

			# Make a sequence object with the longest observed non-ORF transcript 
			# and write it to the output file for that gene.
			my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_CDS");

			$filtered->write_seq($best_cds);

		}

		# Now that we've printed our result for the previous gene, now its 
		# time to reset the loop variables

		# Update the current gene to the first section of the current id
		my @id = split(/\|/, $id);
		$current_gene = $id[0];
		# Reset ORF and CDS variables
		$current_longest_ORF = "0";
		$current_longest_INC = "0";
		$current_longest_CDS = "0";

		# And now I need to tell the script to go back into the loop
		goto CHECKS;

	}

}

# Process last entry
	
# See if any ORFs have been found. If they have...
if ($current_longest_ORF ne "0") {
	
# Make a new sequence object with the current gene as the id
# and the longest ORF as the sequence
my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_ORF");
	
$filtered->write_seq($best_cds);			
	
# Else, if no complete ORFs were found, but an incomplete one was
} elsif ($current_longest_INC ne "0") {

	# Make a sequence object with the longest observed incomplete transcript 
	# and write it to the output file for that gene.
	my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_INC");

	$filtered->write_seq($best_cds);	

} elsif ($current_longest_CDS ne "0") {

	# Make a sequence object with the longest observed non-ORF transcript 
	# and write it to the output file for that gene.
	my $best_cds = Bio::Seq->new(-id => "$current_gene", -seq => "$current_longest_CDS");

	$filtered->write_seq($best_cds);

}		


exit;	

		








