#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;



#--------------Define Misc. starting variables------------------#

# Hold onto special characters and get rid of metacharacter properties
my ($dash, $ast, $X) = ('-', '*', 'X');
$dash = quotemeta ($dash);
$ast = quotemeta ($ast);
$X = quotemeta ($X);

my $pause;

my @complete_convergent_results;
my @complete_parallel_results;

#--------------Collect Command Line Arguments------------------#

# Initialize command line arguments
my ($marsupial_node, $eutherian_node, $min_identity, $suffix) = @ARGV;

unless (@ARGV) {

	print STDERR "This script identifies convergent and parallel amino acids between two nodes of interest.
This version has some hard-coding of the phylogenies in our data sets and thus is not broadly useful to other projects.
\nperl find_parallel_and_convergent_amino_acids_draft4.pl [marsupial] [eutherian] [min_id] [suffix]\n\n";
	exit;

}

my $convergent_results_filename = "convergent_amino_acid_sites_$min_identity\_$suffix.txt";
my $parallel_results_filename = "parallel_amino_acid_sites_$min_identity\_$suffix.txt";

#-------------------------Main Body----------------------------#

# Collect all peptide alignment files in folder
my @files = grep ( -f ,<*AA_comparison.fasta>);

# Iterate through each file in the list of comparison fastas
foreach my $file (@files) {

	# Check the file for overall percentage identity to filter out bad alignments
	my $identity_test = &GetOverallIdentity($file, $min_identity);

	# If the identity test returns true, collect the actual peptide sequences
	if ($identity_test eq 'true') {

		# Split the file name to collect the devil geneID 
		my @file = split('_', $file);
		my $devil_geneID = $file[0];

		# Construct sequence object for marsupial species, canid hash and ancestral nodes hashes
		my ($thylacine_sequence_object, $canidae_sequence_object, $thylacine_MRAN, $canidae_MRAN) = &ConstructSearchStructures($file, $marsupial_node, $eutherian_node);

		# Only continue if suitable sequences were found
		if ($thylacine_sequence_object and $canidae_sequence_object and $thylacine_MRAN and $canidae_MRAN) {

			# Find Parallel amino acids and format result entries
			# Save them to an array.
			my ($convergent_results, $parallel_results) = &FindConvergentAndParallelAminos ($devil_geneID, $thylacine_sequence_object, $canidae_sequence_object, $thylacine_MRAN, $canidae_MRAN);

			#Dereference results arrays
			my @convergent_results = @$convergent_results;
			my @parallel_results = @$parallel_results;


			push(@complete_convergent_results, @convergent_results);

			push(@complete_parallel_results, @parallel_results);

		}

	}
		
}

open (CONV, ">$convergent_results_filename");
open (PARA, ">$parallel_results_filename");

foreach (@complete_convergent_results) { print CONV "$_\n"; }
foreach (@complete_parallel_results) { print PARA "$_\n"; }


#----------------------Subroutine Definitions-----------------------#
###################################################################
sub FindConvergentAndParallelAminos {

	#------------Initialize Subroutine variables----------------#

	my ($devil_geneID, $thylacine_sequence_object, $canidae_sequence_object, $thylacine_MRAN, $canidae_MRAN) = @_;

	# Get ID of MRANs
	my $thylacine_MRAN_id = $thylacine_MRAN->id();
	my $canidae_MRAN_id = $canidae_MRAN->id();

	# Get all protein sequences
	my $thylacine_protein = $thylacine_sequence_object->seq();
	my $thylacine_MRAN_protein = $thylacine_MRAN->seq();
	my $canidae_protein = $canidae_sequence_object->seq();
	my $canidae_MRAN_protein = $canidae_MRAN->seq();

	# Split them into arrays to iterate across
	my @thylacine_protein = split('', $thylacine_protein);
	my @thylacine_MRAN_protein = split('', $thylacine_MRAN_protein);
	my @canidae_protein = split('', $canidae_protein);
	my @canidae_MRAN_protein = split('', $canidae_MRAN_protein);

	# Initialize a count to keep track of the current amino acid
	my $current_amino_number = 0;

	# Initialize arrays to hold results;
	my @convergent_results;
	my @parallel_results;

	#------------Iterate across amino acid positions------------#

	# Parse the marsupial protein sequence one amino acid at a time
	foreach my $thylacine_amino (@thylacine_protein) {

		# Collect all corresponding amino acids in alignment column
		my $thylacine_MRAN_amino = $thylacine_MRAN_protein[$current_amino_number];
		my $canidae_amino = $canidae_protein[$current_amino_number];
		my $canidae_MRAN_amino = $canidae_MRAN_protein[$current_amino_number];

		# If thylacine amino acid is valid...
		if (($thylacine_amino !~ /$dash/) and ($thylacine_amino !~ /$ast/) and ($thylacine_amino !~ /$X/)) {

			# Perform first check: Thylacine == Canidae
			if (($thylacine_amino eq $canidae_amino) 
			# Perform second check: Thylacine =/= Thylacine_MRAN
			and ($thylacine_amino ne $thylacine_MRAN_amino)
			# Perform third check: Canidae =/= Canidae_MRAN
			and ($canidae_amino ne $canidae_MRAN_amino)) {

				

				# Now we know the position is either parallel or convergent
				# Perform the two alternate fourth checks to determine which

				# Convergent: Thylacine_MRAN =/= Canidae_MRAN
				if ($thylacine_MRAN_amino ne $canidae_MRAN_amino) {
		
					# Correct for Perl starting in index 0 by adding 1 to the position
					my $corrected_amino_number = $current_amino_number + 1;

					my @result;				

					push (@result, $devil_geneID);
					push (@result, $corrected_amino_number);
					push (@result, $thylacine_amino);
					push (@result, "$thylacine_MRAN_id");
					push (@result, "$thylacine_MRAN_amino");
					push (@result, "$canidae_MRAN_id");
					push (@result, "$canidae_MRAN_amino");
					my $result = join("\t", @result);

					push (@convergent_results, $result);
					
					
				# Parallel: Thylacine_MRAN == Canidae_MRAN (if they're not not equal, they're equal)
				} else {

					# Correct for Perl starting in index 0 by adding 1 to the position
					my $corrected_amino_number = $current_amino_number + 1;

					my @result;				

					push (@result, $devil_geneID);
					push (@result, $corrected_amino_number);
					push (@result, $thylacine_amino);
					push (@result, "$thylacine_MRAN_id");
					push (@result, "$thylacine_MRAN_amino");
					push (@result, "$canidae_MRAN_id");
					push (@result, "$canidae_MRAN_amino");
					my $result = join("\t", @result);

					push (@parallel_results, $result);
			
				}


			}
				
		}

		$current_amino_number++;

	}

	return (\@convergent_results, \@parallel_results);

}				






###################################################################
sub ConstructSearchStructures {

	my @no_result = undef;

	# Most Recent Ancestral Nodes (MRAN)
	my $thylacine_sequence_object;
	my $canidae_sequence_object;
	my $thylacine_MRAN;
	my $canidae_MRAN;

	#------------Initialize Subroutine variables----------------#
	# Collect arguments
	my $file = $_[0];
	my $marsupial_node = $_[1];
	my $eutherian_node = $_[2];

	# Define marsupial heirarchy
	my @marsupial_ranked_nodes = ("Tcyn", "Dasyuromorphia", "Australidelphia", "Marsupials");
	# Define eutherian heirarchy
	my @eutherian_ranked_nodes;
	if ($eutherian_node eq "Canis") { @eutherian_ranked_nodes = ("Canis", "Canidae", "Caniforms", "Carnivora", "CU", "CUC", "CUCE", "Eutherians") }
	elsif ($eutherian_node eq "Vulpes") { @eutherian_ranked_nodes = ("Vulpes", "Canidae", "Caniforms", "Carnivora", "CU", "CUC", "CUCE", "Eutherians") }
	else { @eutherian_ranked_nodes = ("Canidae", "Caniforms", "Carnivora", "CU", "CUC", "CUCE", "Eutherians") }

	# Determine which node the marsupial MRAN search will begin after 
	# based on user's choice
	my $count = 0;
	my $marsupial_index;

	foreach (@marsupial_ranked_nodes) {
		if ($marsupial_node eq $_) {
			$marsupial_index = $count;
		}
		$count++;
	}

	# Determine which node the eutherian MRAN search will begin after 
	# based on user's choice 
	$count = 0;
	my $eutherian_index;

	foreach (@eutherian_ranked_nodes) {
		if ($eutherian_node eq $_) {
			$eutherian_index = $count;
		}
		$count++;
	}
	

	# Initialize hashes to hold sequence objects
	my %canidae_MRAN_hash;
	my %thylacine_MRAN_hash;

	#----------- Construct hashes and sequence object-----------#

	# Open the multifasta peptide alignment	
	my $multifasta = Bio::SeqIO->new(-file => "$file", -format => 'fasta', -alphabet => 'protein');

	# Check each sequence in the alignment fasta
	while (my $sequence_object = $multifasta->next_seq) {

		# Get header
		my $header = $sequence_object->display_id();
		chomp $header;

		# Get sequence
#		my $sequence = $sequence_object->seq();

		# If the header matches thlyacine, grab the sequence for parsing later
		if ($header =~ /$marsupial_node/ig) {
			
			$thylacine_sequence_object = $sequence_object;

		# If the header matches canidae, grab the sequence for parsing later
		} elsif ($header =~ /$eutherian_node/ig) {

			$canidae_sequence_object = $sequence_object;

		} else {

			for (my $i = $marsupial_index + 1; $i < @marsupial_ranked_nodes; $i++) {

				if ($header =~ /$marsupial_ranked_nodes[$i]/) { $thylacine_MRAN_hash{$header} = $sequence_object }

			}
		
			for (my $i = $eutherian_index + 1; $i < @eutherian_ranked_nodes; $i++) {
			
				if ($header =~ /$eutherian_ranked_nodes[$i]/) { $canidae_MRAN_hash{$header} = $sequence_object }
			
			}

		}

	}

	# Unless a sequence for thylacine and canidae were found, return a blank line
	# for this iteration.
	unless ($thylacine_sequence_object) {

		print STDERR "$file: $marsupial_node sequences not present in alignment\n";
		last;
		return (@no_result);

	}
	
	unless ($canidae_sequence_object) {

		print STDERR "$file: $eutherian_node sequences not present in alignment\n";
		
		return (@no_result);

	}
		
	# Figure out which is the most recent ancestral node to thylacine
	if (exists $thylacine_MRAN_hash{"Dasyuromorphia"}) { $thylacine_MRAN = $thylacine_MRAN_hash{"Dasyuromorphia"} }
	elsif (exists $thylacine_MRAN_hash{"Australidelphia"}) { $thylacine_MRAN = $thylacine_MRAN_hash{"Australidelphia"} }
	elsif (exists $thylacine_MRAN_hash{"Marsupials"}) { $thylacine_MRAN = $thylacine_MRAN_hash{"Marsupials"} }
	else { print STDERR "$file: No suitable ancestor to Thylacine in alignment\n"; return (@no_result) }

	# Figure out which is the most recent ancestral node to canidae
	if (exists $canidae_MRAN_hash{"Canidae"}) { $canidae_MRAN = $canidae_MRAN_hash{"Canidae"} }
	elsif (exists $canidae_MRAN_hash{"Caniforms"}) { $canidae_MRAN = $canidae_MRAN_hash{"Caniforms"} }
	elsif (exists $canidae_MRAN_hash{"Carnivora"}) { $canidae_MRAN = $canidae_MRAN_hash{"Carnivora"} }
	elsif (exists $canidae_MRAN_hash{"CU"}) { $canidae_MRAN = $canidae_MRAN_hash{"CU"} }
	elsif (exists $canidae_MRAN_hash{"CUC"}) { $canidae_MRAN = $canidae_MRAN_hash{"CUC"} }
	elsif (exists $canidae_MRAN_hash{"CUCE"}) { $canidae_MRAN = $canidae_MRAN_hash{"CUCE"} }
	elsif (exists $canidae_MRAN_hash{"Eutherians"}) { $canidae_MRAN = $canidae_MRAN_hash{"Eutherians"} }
	else { print STDERR "$file: No suitable ancestor to Canidae in alignment\n"; return (@no_result) }

#	my $a = $thylacine_sequence_object->id();
#	my $b = $canidae_sequence_object->id();
#	my $c = $thylacine_MRAN->id(); 
#	my $d = $canidae_MRAN->id();

#	print STDOUT "$file\n$a\n$b\n$c\n$d\n\n";


	return ($thylacine_sequence_object, $canidae_sequence_object, $thylacine_MRAN, $canidae_MRAN);

}

sub GetOverallIdentity {

	# Collect arguments
	my ($file, $min_identity) = @_;

	# Initialize a variable to hold the overall alignment identity
	my $identity;

	# And results of the idenity test
	my $identity_test;
	
	# Open the file as a fasta alignment to check identity scores
	my $alignIO_object = Bio::AlignIO->new(-file => "$file", -format => "fasta");

	# Check each sequence in the alignment fasta
	while (my $simple_aln_object = $alignIO_object->next_aln) {

		# Collect the overall percentage idenity of the alignment file
		$identity = $simple_aln_object->overall_percentage_identity;

	}

	# See if the overall alignment identity is at or above the minimum
	if ($identity >= $min_identity) {
	
		$identity_test = 'true';

	} else {

		$identity_test = 'false';

	}

	# Return the results of the identity test
	return $identity_test;
	
}















