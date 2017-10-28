#!/usr/bin/perl

use strict;
use warnings; 

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;

#--------------Define starting variables------------------#

my $pause;
my $USAGE;

#-----------------Make output file names--------------------#
# From the phylip alignment used in CodeML, collect the basename
# and use it to assign sensible names to output files.

# Collect main phylip alignment file name
my @phylip_alignment_filename = grep ( -f ,<ENS*.phy>);
my $phylip_alignment_filename = $phylip_alignment_filename[0];

# Strip phylip alignment file name down to basename
@phylip_alignment_filename = split (/\./, $phylip_alignment_filename);
my $basename = $phylip_alignment_filename[0];


# Use basename to make output custom fixed tree structure filename
my $custom_fixed_tree_filename = "$basename" . "_custom_fixed_tree.tree";


#--------Parse phylip file to collect sequence IDs----------#

# Open phylip file
unless( open( PHYLIP, $phylip_alignment_filename ) ) {

	print "Can't open $phylip_alignment_filename\n";

}

# Make a hash to hold the list of ids from species in the phylip alignment file
# together with their respective sequences. Only some will be output at the end.
my %species;

# Collect non-header lines
#while (my $line = <PHYLIP>) {
my @phylip = <PHYLIP>;

foreach my $line (@phylip) {
	# If the line is not the numerical header...
	unless ($line =~ /^\d/) {

		chomp ($line);

		# Split up the ID and sequence
		my @split_line = split(' ', $line);

		# Add each id to the hash as a value and its sequence a the key
		$species{$split_line[0]} = $split_line[1];

	}

}

close PHYLIP;


#------------------Construct the branches of the tree----------------#
# All output trees will have 3 branches (trifurcating)
# Marsupial and Eutherian branches will always be present
# The third branch will vary depending on the other species present
# If the platypus is present (phylip file contains POG)
# Then it will be the third branch and no birds will be retained
# If the platypus isn't present but birds are (phyip file contains BOG)
# Then a tree for the birds will be  constructed and used as the third branch
# If neither birds or platypus are present (phylip file contains NOG)
# Then the lowest species on the eutherian tree will be droped down to third branch
# Must determine output before constructing eutherian tree
my $outgroup_subtree;
my $eutherian_subtree;
my $marsupial_subtree;

# First determine the outgroup 
#if ($basename =~ /POG/g) {$outgroup_subtree = "ENSOANG"}
#elsif ($basename =~ /BOG/g) {$outgroup_subtree = &ParseBirds(%species)}
#elsif ($basename =~ /NOG/g) {$outgroup_subtree = "false"}	
if (exists $species{"ENSOANG"}) {
	$outgroup_subtree = "ENSOANG";
	my @new_phylip;
	# If platypus is present as well as birds, we'll need to remove them from the file
	my $remove_count = 0;
	my $header;
	open (REDUCED, ">$phylip_alignment_filename");
	foreach my $line (@phylip) {
		chomp ($line);
		if ($line =~ /^\d/) {
			$header = $line; 
		} elsif (($line =~ /ENSGALG/g) or ($line =~ /ENSTGUG/g) or ($line=~ /ENSFALG/g)) {
			$remove_count++;
		} else {
			push(@new_phylip, $line);
		}
	}
	my @split_header = split(' ', $header);
	$split_header[0] = $split_header[0] - $remove_count;
	$header = join ("\t", @split_header);
	unshift (@new_phylip, $header);
	my $new_phylip = join ("\n", @new_phylip);
	print REDUCED $new_phylip;
	close REDUCED;

	# As well, we need to remove them from the species hash
}
elsif ((exists $species{"ENSGALG"}) or (exists $species{"ENSTGUG"}) or (exists $species{"ENSFALG"})) {
	$outgroup_subtree = &ParseBirds(%species);
	}
else {$outgroup_subtree = "false"}

#print "$outgroup_subtree\n";
chomp ($outgroup_subtree);

# Then construct the eutherian tree 
$eutherian_subtree = &ParseEutherians(%species);

# Then construct marsupial tree
$marsupial_subtree = &ParseMarsupials(%species);

# Assemble the custom fixed tree
my $custom_fixed_tree;

if ($outgroup_subtree eq "false") {
	
	my @split_tree = split('', $eutherian_subtree);
	pop @split_tree;
	shift@split_tree;
	$eutherian_subtree = join('', @split_tree);
	$custom_fixed_tree = "($marsupial_subtree,$eutherian_subtree);";

} else {

	$custom_fixed_tree = "($outgroup_subtree,$marsupial_subtree,$eutherian_subtree);";

}

open (OUTPUT, ">$custom_fixed_tree_filename");

print OUTPUT "$custom_fixed_tree";

#----------------------Subroutine Definitions-----------------------#

####################################################################
sub ParseMarsupials {
	
	# Collect species hash
	my %species_hash = @_;

	my @dasyuromorphia_subtree;
	my @australidelphia_subtree;
	my @marsupial_subtree;

	# dasyuromorphia collect
	if (exists $species_hash{"ENSSHAG"}) { my $x = "ENSSHAG"; push (@dasyuromorphia_subtree, $x) }
	if (exists $species_hash{"Tcyn"}) { my $x = "Tcyn"; push (@dasyuromorphia_subtree, $x) }
	# dasyuromorphia build
	my $dasyuromorphia_subtree;
	if (scalar(@dasyuromorphia_subtree) > 1) {
		$dasyuromorphia_subtree = join (',', @dasyuromorphia_subtree);
		$dasyuromorphia_subtree = "($dasyuromorphia_subtree)";
		push (@australidelphia_subtree, $dasyuromorphia_subtree);
	} elsif (scalar(@dasyuromorphia_subtree) == 1) {
		$dasyuromorphia_subtree = "@dasyuromorphia_subtree";
		#$dasyuromorphia_subtree = "($dasyuromorphia_subtree)";
		push (@australidelphia_subtree, $dasyuromorphia_subtree);
	}

	# Wallaby collect
	my $wallaby;
	if (exists $species_hash{"ENSMEUG"}) {
		$wallaby = "ENSMEUG";
		push (@australidelphia_subtree, $wallaby);
	}

	# Australidelphia build
	my $australidelphia_subtree;
	if (scalar(@australidelphia_subtree) > 1) {
		$australidelphia_subtree = join (',', @australidelphia_subtree);
		$australidelphia_subtree = "($australidelphia_subtree)";
		push (@marsupial_subtree, $australidelphia_subtree);
	} elsif (scalar(@australidelphia_subtree) == 1) {
		$australidelphia_subtree = "@australidelphia_subtree";
		#$australidelphia_subtree = "($australidelphia_subtree)";
		push (@marsupial_subtree, $australidelphia_subtree);
	}

	# opossum collect
	my $opossum;
	if (exists $species_hash{"ENSMODG"}) {
		$opossum = "ENSMODG";
		push (@marsupial_subtree, $opossum);
	}

	# marsupial build
	my $marsupial_subtree;
	if (scalar(@marsupial_subtree) > 1) {
		$marsupial_subtree = join (',', @marsupial_subtree);
		$marsupial_subtree = "($marsupial_subtree)";
	} elsif (scalar(@marsupial_subtree) == 1) {
		$marsupial_subtree = "@marsupial_subtree";
		#$marsupial_subtree = "($marsupial_subtree)";
	}


	# Initialize variable to build into the marsupial subtree
	# Thylacine and Devil will always be present
	# so this part of the tree can be constant
#	my $marsupial_subtree = "(Tcyn,ENSSHAG)";

#	if (exists $species_hash{"ENSMEUG"}) { $marsupial_subtree = '(ENSMEUG,' . "$marsupial_subtree" . ")" }
#	if (exists $species_hash{"ENSMODG"}) { $marsupial_subtree = '(ENSMODG,' . "$marsupial_subtree" . ")" }
#	chomp ($marsupial_subtree);
	return ($marsupial_subtree);

}

sub ParseEutherians {

	# Collect species hash
	my %species_hash = @_;

	# Subtrees within eutheria
#		my $laurasiatheria_subtree;
		my @CU_subtree;
				my @carnivora_subtree;
					my @caniform_subtree;
						my @canid_subtree;
							my @canis_subtree;
								my @WD_subtree;
								my @WDC_subtree;
							my @vulpes_subtree;
						my $panda_subtree;
					my @feliform_subtree;
				my @ungulate_subtree;
					my @artiodactyla_subtree;
					my $perissodactyla_subtree;
		my @CUC_subtree;
				my @chiroptera_subtree;
					# microbat
					# megabat
		my @CUCE_subtree;
			my @euarchontaglires_subtree;
				my @euarchonta_subtree;
				my @glires_subtree;
		my @eutherian_subtree;
			my $afrotheria_subtree;


	#----------------------Carnivora------------------------#
	# Carnivora collect
		# Caniform collect
			# Canid collect
				# Canis collect -> Canis phylogeny based on nuclear gene analysis in 
				#               -> (Koepfli et al. 2015), doi:10.1016/j.cub.2015.06.060
#				if (exists $species_hash{"Clup"}) { my $x = "(Clup,ENSCAFG)"; push (@canis_subtree, $x) }
#				if (exists $species_hash{"Clup"}) { $canis_subtree = "(Clup,ENSCAFG)" }
#				else { my $x = "ENSCAFG"; push (@canis_subtree, $x) }
#				else { $canis_subtree = "(ENSCAFG)" }
#				if (exists $species_hash{"Clat"}) { my $x = "Clat"; push (@canis_subtree, $x) }
#				if (exists $species_hash{"Clat"}) { $canis_subtree = "(Clat,$canis_subtree)" }
#				if (exists $species_hash{"Caur"}) { my $x = "Caur"; push (@canis_subtree, $x) }
#				if (exists $species_hash{"Caur"}) { $canis_subtree = "(Caur,$canis_subtree)" }
				# Canis build
#				my $canis_subtree = join (',', @canis_subtree);
#				$canis_subtree = "($canis_subtree)";
#				push (@canid_subtree, $canis_subtree);	

				# WD collect
				if (exists $species_hash{"ENSCAFG"}) { my $x = "ENSCAFG"; push (@WD_subtree, $x) }
				if (exists $species_hash{"Clup"}) { my $x = "Clup"; push (@WD_subtree, $x) }
				# WD build
				my $WD_subtree;
				if (scalar(@WD_subtree) > 1) {
					$WD_subtree = join (',', @WD_subtree);
					$WD_subtree = "($WD_subtree)";
					push (@WDC_subtree, $WD_subtree);
				} elsif (scalar(@WD_subtree) == 1) {
					$WD_subtree = "@WD_subtree";
					#$WD_subtree = "($WD_subtree)";
					push (@WDC_subtree, $WD_subtree);
				}

				# Coyote collect
				my $coyote;
				if (exists $species_hash{"Clat"}) {
					$coyote = "Clat";
					push (@WDC_subtree, $coyote);
				}

				# WDC build
				my $WDC_subtree;
				if (scalar(@WDC_subtree) > 1) {
					$WDC_subtree = join (',', @WDC_subtree);
					$WDC_subtree = "($WDC_subtree)";
					push (@canis_subtree, $WDC_subtree);
				} elsif (scalar(@WDC_subtree) == 1) {
					$WDC_subtree = "@WDC_subtree";
					#$WDC_subtree = "($WDC_subtree)";
					push (@canis_subtree, $WDC_subtree);
				}
				

				# Jackal collect
				my $jackal;
				if (exists $species_hash{"Caur"}) {
					$jackal = "Caur";
					push (@canis_subtree, $jackal);
				}

				# Canis build
				my $canis_subtree;
				if (scalar(@canis_subtree) > 1) {
					$canis_subtree = join (',', @canis_subtree);
					$canis_subtree = "($canis_subtree)";
					push (@canid_subtree, $canis_subtree);
				} elsif (scalar(@canis_subtree) == 1) {
					$canis_subtree = "@canis_subtree";
					#$canis_subtree = "($canis_subtree)";
					push (@canid_subtree, $canis_subtree);
				}

#print "@canid_subtree\n";
				# Vulpes collect
				if (exists $species_hash{"Vvul"}) { my $x = "Vvul"; push (@vulpes_subtree, $x) }
				if (exists $species_hash{"Vlag"}) { my $x = "Vlag"; push (@vulpes_subtree, $x) }
				# Vulpes build
				my $vulpes_subtree;
				if (scalar(@vulpes_subtree) > 1) {
					$vulpes_subtree = join (',', @vulpes_subtree);
					$vulpes_subtree = "($vulpes_subtree)";
					push (@canid_subtree, $vulpes_subtree);
				} elsif (scalar(@vulpes_subtree) == 1) {
					$vulpes_subtree = "@vulpes_subtree";
					#$vulpes_subtree = "($vulpes_subtree)";
					push (@canid_subtree, $vulpes_subtree);
				}
#print scalar(@canid_subtree),"\n";				
#print "$vulpes_subtree\n";

			# Canid build
			my $canid_subtree;
			if (scalar(@canid_subtree) > 1) {
				$canid_subtree = join (',', @canid_subtree);
				$canid_subtree = "($canid_subtree)";
				push (@caniform_subtree, $canid_subtree);
			} elsif (scalar(@canid_subtree) == 1) {
				$canid_subtree = "@canid_subtree";
				#$canid_subtree = "($canid_subtree)";
				push (@caniform_subtree, $canid_subtree);
			}
#print "$canid_subtree\n";			

			# Panda collect
			if (exists $species_hash{"ENSAMEG"}) {
				$panda_subtree = "ENSAMEG";
				# Panda build
				push (@caniform_subtree, $panda_subtree);
			}

		# Caniform build
		my $caniform_subtree;
		if (scalar(@caniform_subtree) > 1) {
#		if (($vulpes_subtree) and ($canis_subtree)) {$caniform_subtree = "($caniform_subtree)"}
#		$caniform_subtree = "($caniform_subtree)";
			$caniform_subtree = join (',', @caniform_subtree);
			$caniform_subtree = "($caniform_subtree)";
			push (@carnivora_subtree, $caniform_subtree);
		} elsif (scalar(@caniform_subtree) == 1) {
			$caniform_subtree = "@caniform_subtree";
			#$caniform_subtree = "($caniform_subtree)";
			push (@carnivora_subtree, $caniform_subtree);
		}


		# Feliform collect
		if (exists $species_hash{"ENSFCAG"}) { my $x = "ENSFCAG"; push (@feliform_subtree, $x) }
		if (exists $species_hash{"ENSMPUG"}) { my $x = "ENSMPUG"; push (@feliform_subtree, $x) }
		# Feliform build
		my $feliform_subtree;
		if (scalar(@feliform_subtree) > 1) {
			$feliform_subtree = join (',', @feliform_subtree);
			$feliform_subtree = "($feliform_subtree)";
			push (@carnivora_subtree, $feliform_subtree);
		} elsif (scalar(@feliform_subtree) == 1) {
			$feliform_subtree = "@feliform_subtree";
			#$feliform_subtree = "($feliform_subtree)";
			push (@carnivora_subtree, $feliform_subtree);
		}
#		$feliform_subtree = "($feliform_subtree)";

	# Carnivora build
#print scalar(@carnivora_subtree),"\n";
	my $carnivora_subtree;
	if (scalar(@carnivora_subtree) > 1) {
		$carnivora_subtree = join (',', @carnivora_subtree);
		$carnivora_subtree = "($carnivora_subtree)";
		push (@CU_subtree, $carnivora_subtree);
	} elsif (scalar(@carnivora_subtree) == 1) {
		$carnivora_subtree = "@carnivora_subtree";
		#$carnivora_subtree = "($carnivora_subtree)";
		push (@CU_subtree, $carnivora_subtree);
	}
#	my $carnivora_subtree = join (',', @carnivora_subtree);
#	if (($feliform_subtree) and ($caniform_subtree)) {$carnivora_subtree = "($carnivora_subtree)"}
#	$carnivora_subtree = "($carnivora_subtree)";

	#----------------------Ungulates------------------------#
	# Ungulate collect
		# artiodactyla collect
		if (exists $species_hash{"ENSBTAG"}) { my $x = "ENSBTAG"; push (@artiodactyla_subtree, $x) }
		if (exists $species_hash{"ENSOARG"}) { my $x = "ENSOARG"; push (@artiodactyla_subtree, $x) }
		# artiodactyla build
		my $artiodactyla_subtree;
		if (scalar(@artiodactyla_subtree) > 1) {
			$artiodactyla_subtree = join (',', @artiodactyla_subtree);
			$artiodactyla_subtree = "($artiodactyla_subtree)";
			push (@ungulate_subtree, $artiodactyla_subtree);
		} elsif (scalar(@artiodactyla_subtree) == 1) {
			$artiodactyla_subtree = "@artiodactyla_subtree";
			#$artiodactyla_subtree = "($artiodactyla_subtree)";
			push (@ungulate_subtree, $artiodactyla_subtree);
		}
	
		# perissodactyla (horse) collect
		if (exists $species_hash{"ENSECAG"}) {
			$perissodactyla_subtree = "ENSECAG";
			push (@ungulate_subtree, $perissodactyla_subtree);
		}

		

	# Ungulate build
	my $ungulate_subtree;
	if (scalar(@ungulate_subtree) > 1) {
		$ungulate_subtree = join (',', @ungulate_subtree);
		$ungulate_subtree = "($ungulate_subtree)";
		push (@CU_subtree, $ungulate_subtree);
	} elsif (scalar(@ungulate_subtree) == 1) {
		$ungulate_subtree = "@ungulate_subtree";
		#$ungulate_subtree = "($ungulate_subtree)";
		push (@CU_subtree, $ungulate_subtree);
	}

	#----------------------Chiropterans------------------------#
	# Chiropteran collect
	if (exists $species_hash{"ENSMLUG"}) { my $x = "ENSMLUG"; push (@chiroptera_subtree, $x) }
	if (exists $species_hash{"ENSPVAG"}) { my $x = "ENSPVAG"; push (@chiroptera_subtree, $x) }
	
	# Chiropteran build
	my $chiroptera_subtree;
	if (scalar(@chiroptera_subtree) > 1) {
		$chiroptera_subtree = join (',', @chiroptera_subtree);
		$chiroptera_subtree = "($chiroptera_subtree)";
		push (@CUC_subtree, $chiroptera_subtree);
	} elsif (scalar(@chiroptera_subtree) == 1) {
		$chiroptera_subtree = "@chiroptera_subtree";
		#$chiroptera_subtree = "($chiroptera_subtree)";
		push (@CUC_subtree, $chiroptera_subtree);
	}
	#----------------------Euarchontagliers------------------------#

		# Euarchonta collect
		if (exists $species_hash{"ENSG0"}) { my $x = "ENSG0"; push (@euarchonta_subtree, $x) }
		if (exists $species_hash{"ENSCJAG"}) { my $x = "ENSCJAG"; push (@euarchonta_subtree, $x) }
		# Euarchonta build
		my $euarchonta_subtree;
		if (scalar(@euarchonta_subtree) > 1) {
			$euarchonta_subtree = join (',', @euarchonta_subtree);
			$euarchonta_subtree = "($euarchonta_subtree)";
			push (@euarchontaglires_subtree, $euarchonta_subtree);
		} elsif (scalar(@euarchonta_subtree) == 1) {
			$euarchonta_subtree = "@euarchonta_subtree";
			#$euarchonta_subtree = "($euarchonta_subtree)";
			push (@euarchontaglires_subtree, $euarchonta_subtree);
		}

		# Gliers collect
		if (exists $species_hash{"ENSCPOG"}) { my $x = "ENSCPOG"; push (@glires_subtree, $x) }
		if (exists $species_hash{"ENSOCUG"}) { my $x = "ENSOCUG"; push (@glires_subtree, $x) }
		# Glires build
		my $glires_subtree;
		if (scalar(@glires_subtree) > 1) {
			$glires_subtree = join (',', @glires_subtree);
			$glires_subtree = "($glires_subtree)";
			push (@euarchontaglires_subtree, $glires_subtree);
		} elsif (scalar(@glires_subtree) == 1) {
			$glires_subtree = "@glires_subtree";
			#$glires_subtree = "($glires_subtree)";
			push (@euarchontaglires_subtree, $glires_subtree);
		}

	# Euarchontaglires build
	my $euarchontaglires_subtree;
	if (scalar(@euarchontaglires_subtree) > 1) {
		$euarchontaglires_subtree = join (',', @euarchontaglires_subtree);
		$euarchontaglires_subtree = "($euarchontaglires_subtree)";
		push (@CUCE_subtree, $euarchontaglires_subtree);
	} elsif (scalar(@euarchontaglires_subtree) == 1) {
		$euarchontaglires_subtree = "@euarchontaglires_subtree";
		#$euarchontaglires_subtree = "($euarchontaglires_subtree)";
		push (@CUCE_subtree, $euarchontaglires_subtree);
	}

	#----------------------Afrotherian------------------------#

	if (exists $species_hash{"ENSLAFG"}) {
		$afrotheria_subtree = "ENSLAFG";
		push (@eutherian_subtree, $afrotheria_subtree);
	}

	# Build CU tree
	my $CU_subtree;
	if (scalar(@CU_subtree) > 1) {
		$CU_subtree = join (',', @CU_subtree);
		$CU_subtree = "($CU_subtree)";
		push (@CUC_subtree, $CU_subtree);
	} elsif (scalar(@CU_subtree) == 1) {
		$CU_subtree = "@CU_subtree";
		#$CU_subtree = "($CU_subtree)";
		push (@CUC_subtree, $CU_subtree);
	}
	
	# Build CUC tree
	my $CUC_subtree;
	if (scalar(@CUC_subtree) > 1) {
		$CUC_subtree = join (',', @CUC_subtree);
		$CUC_subtree = "($CUC_subtree)";
		push (@CUCE_subtree, $CUC_subtree);
	} elsif (scalar(@CUC_subtree) == 1) {
		$CUC_subtree = "@CUC_subtree";
		#$CUC_subtree = "($CUC_subtree)";
		push (@CUCE_subtree, $CUC_subtree);
	}
	
	# Build CUCE tree
	my $CUCE_subtree;
	if (scalar(@CUCE_subtree) > 1) {
		$CUCE_subtree = join (',', @CUCE_subtree);
		$CUCE_subtree = "($CUCE_subtree)";
		push (@eutherian_subtree, $CUCE_subtree);
	} elsif (scalar(@CUCE_subtree) == 1) {
		$CUCE_subtree = "@CUCE_subtree";
		#$CUCE_subtree = "($CUCE_subtree)";
		push (@eutherian_subtree, $CUCE_subtree);
	}

	# Build eutherian tree
	my $eutherian_subtree;
	if (scalar(@eutherian_subtree) > 1) {
		$eutherian_subtree = join (',', @eutherian_subtree);
		$eutherian_subtree = "($eutherian_subtree)";

	} elsif (scalar(@eutherian_subtree) == 1) {
		$eutherian_subtree = "@eutherian_subtree";
		#$eutherian_subtree = "($eutherian_subtree)";

	}

	# Build laurasiatherian tree
	
	# There will always be carnivora in the tree, so start there
#	$eutherian_subtree = $carnivora_subtree;
	
	# Check to see if any ungulates are present. If they are, add branch
#	if ($ungulate_subtree) { $eutherian_subtree = "($ungulate_subtree,$eutherian_subtree)" }

	# Check to see if any chiropterans are present. If they are, add branch
#	if ($chiroptera_subtree) { $eutherian_subtree = "($chiroptera_subtree,$eutherian_subtree)" }

	# Check to see if any euarchontaglires are present. If they are, add branch
#	if ($euarchontaglires_subtree) { $eutherian_subtree = "($euarchontaglires_subtree,$eutherian_subtree)" }

	# Check to see if any afrotherians are present. If they are, add branch
#	if ($afrotheria_subtree) { $eutherian_subtree = "($afrotheria_subtree,$eutherian_subtree)" }

	chomp ($eutherian_subtree);
	
	return "$eutherian_subtree";

}



sub ParseBirds {

	
	# Collect species hash
	my %species_hash = @_;

	# Initialize variables to build into the bird subtree and
	# recursive subtrees
	my @bird_subtree;
	my @passeriform_subtree;
	my $galliform_subtree;

	# Passeriform collect
		# Passeriform subtree
		if (exists $species_hash{"ENSFALG"}) { my $x = "ENSFALG"; push (@passeriform_subtree, $x) }
		if (exists $species_hash{"ENSTGUG"}) { my $x = "ENSTGUG"; push (@passeriform_subtree, $x) }
		# Passeriform build
		my $passeriform_subtree;
		if (scalar(@passeriform_subtree) > 1) {
			$passeriform_subtree = join (',', @passeriform_subtree);
			$passeriform_subtree = "($passeriform_subtree)";
			push (@bird_subtree, $passeriform_subtree);
		} elsif (scalar(@passeriform_subtree) == 1) {
			$passeriform_subtree = "@passeriform_subtree";
			#$passeriform_subtree = "($passeriform_subtree)";
			push (@bird_subtree, $passeriform_subtree);
		}
	
		# Galliform subtree
		if (exists $species_hash{"ENSGALG"}) {
			$galliform_subtree = "ENSGALG";
			push (@bird_subtree, $galliform_subtree);
		}


		

	# Ungulate build
	my $bird_subtree;
	if (scalar(@bird_subtree) > 1) {
		$bird_subtree = join (',', @bird_subtree);
		$bird_subtree = "($bird_subtree)";
	} elsif (scalar(@bird_subtree) == 1) {
		$bird_subtree = "@bird_subtree";
		#$bird_subtree = "($bird_subtree)";
	}

#	print "$bird_subtree\n";
	return ($bird_subtree);

}



#	if (exists $species_hash{"PLACEHOLDER"}) {
#	if (exists $species_hash{"PLACEHOLDER"}) {
#	if (exists $species_hash{"PLACEHOLDER"}) {
#	if (exists $species_hash{"PLACEHOLDER"}) {



























