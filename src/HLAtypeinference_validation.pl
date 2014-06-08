#!/usr/bin/perl

BEGIN {
	push(@INC, 'Toolbox');
}

# Comments denoted by ALL-HLA refer to changes which would become necessary if we wanted
# to create a graph for all HLA loci at the same time.

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;
use Sys::Hostname;
use Commons;

use simpleHLA;

my $debug = 0;
my $all_2_dig = 0;
my $only_4_dig = 0;
my $T = 0.0;
my $reference_file;

my $panel_name;
my $training_snp_file;
GetOptions (
	'panel_name:s' => \$panel_name,
	
	'reference:s' => \$reference_file,
	'all_2_dig:s' => \$all_2_dig,
	'only_4_dig:s' => \$only_4_dig,
	'T:s' => \$T,
	'training_snp_file:s' => \$training_snp_file,
 );

 # read reference dataset
my %reference_data;
open(REFERENCE, "<", $reference_file) or die "Cannot open $reference_file";
my $headerLine = <REFERENCE>;
chomp($headerLine);
$headerLine =~ s/\n//g;
$headerLine =~ s/\r//g;
my @header_fields = split(/\s/, $headerLine);
while(<REFERENCE>)
{
	my $line = $_;
	chomp($line);
	
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	
	my @fields = split(/ /, $line);
	my %line = (mesh @header_fields, @fields);
	
	my $primary_key = $line{'IndividualID'};
	$reference_data{$primary_key} = \%line;
}
close(REFERENCE);

# read imputations
my %imputed_HLA;
my %imputed_HLA_Calls;

print qq(
Validation
	Panel name: $panel_name
	T = $T
	all_2_dig: $all_2_dig
	only_4_dig: $only_4_dig
	Reference: $reference_file
);

my $output_folder = '../tmp2/'.$panel_name;
unless(-e $output_folder)
{
	die "Cannot find expected panel output folder: $output_folder";
}

my $pattern = $output_folder.'/*';
my @files = glob($pattern);
@files = grep {$_ =~ /(^|\/|\\)PP_(.+)_bestguess.txt/} @files;
my %loci_present;
my @loci;
my %quality_measures;

foreach my $file (@files)
{
	die unless ($file =~ /PP_(\w+)_bestguess.txt/);
	my $locus = $1;
	
	#next unless (($locus eq 'A') || ($locus eq 'DQB') || ($locus eq 'DRB') || ($locus eq 'B'));
	if($loci_present{$locus})
	{
		die "More than one file for locus $locus?";
	}
	push(@loci, $locus);
	$loci_present{$locus} = 1;
	
	print "Read file $file\n";
	
	open(IMPUTATIONS, "<", $file) or die "Cannot open $file";
	my $header = <IMPUTATIONS>;
	chomp($header);
	my @header_fields = split(/\s+/, $header);
	while(<IMPUTATIONS>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/\s+/, $line);
		my %linehash = (mesh @header_fields, @fields);
		
		my $indivID = $linehash{'IndividualID'};
		my $chromo = $linehash{'Chromosome'};
	
		my $hla_val = $linehash{'HLA'.$locus};
		unless($hla_val)
		{
			$hla_val = $linehash{'KIR'.$locus};
		}
		($linehash{'HLA'.$locus} xor $linehash{'KIR'.$locus}) or die;
		$hla_val or die "Value for $locus not present in $file, line $line -- expected field name $locus";
		
		die "No quality field?" unless ((exists $linehash{'Q2'} ) or  (exists $linehash{'Q'} ))  ;
		my $Q;
		if(exists $linehash{'Q2'})
		{
			$Q = $linehash{'Q2'};
		}
		else
		{
			$Q = $linehash{'Q'};			
		}
		
		if(exists $linehash{'P'})
		{
			$quality_measures{$indivID}{P} = $linehash{'P'}+0.0;
		}
		
		# $imputed_HLA_Calls{$locus}{sum}++;
		if($Q < $T)
		{
			$imputed_HLA{$locus}{$indivID}{$chromo} = '0000';			
		}
		else
		{
			$imputed_HLA{$locus}{$indivID}{$chromo} = $hla_val;
			# $imputed_HLA_Calls{$locus}{called}++;
		}
	}
	close(IMPUTATIONS);
}

my %training_allele_frequencies;
if($training_snp_file)
{
	open(TRAINING, "<", $training_snp_file) or die "Cannot open $training_snp_file";
	my $headerLine = <TRAINING>;
	chomp($headerLine);
	$headerLine =~ s/\n//g;
	$headerLine =~ s/\r//g;
	my @header_fields = split(/\s/, $headerLine);
	my %complete_locus_names;
	foreach my $locus (keys %loci_present)
	{
		my $f = 0;
		if('KIR'.$locus ~~ \@header_fields)
		{
			$complete_locus_names{$locus} = 'KIR'.$locus;
			$f++;
		}
		if('HLA'.$locus ~~ \@header_fields)
		{
			$complete_locus_names{$locus} = 'HLA'.$locus;
			$f++;
		}
		die unless($f == 1);
	}
	while(<TRAINING>)
	{
		my $line = $_;
		chomp($line);
		
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		
		my @fields = split(/ /, $line);
		my %line = (mesh @header_fields, @fields);
		
		foreach my $locus (keys %loci_present)
		{
			my @alleles = split(/\//, $line{$complete_locus_names{$locus}});
			foreach my $allele (@alleles)
			{
				next if ($allele =~ /\?/);
				if($all_2_dig)
				{
					$a = &simpleHLA::HLA_2digit($a);
				}
				$training_allele_frequencies{$locus}{$allele}++;
			}
		}
	}
	close(TRAINING);
}

foreach my $locus (keys %training_allele_frequencies)
{
	my $a_sum = 0;
	foreach my $allele (keys %{$training_allele_frequencies{$locus}})
	{
		$a_sum += $training_allele_frequencies{$locus}{$allele};
	}
	foreach my $allele (keys %{$training_allele_frequencies{$locus}})
	{
		$training_allele_frequencies{$locus}{$allele} = $training_allele_frequencies{$locus}{$allele}/$a_sum;
	}	
}

my %lls_correct_individuals;
my %lls_incorrect_individuals;

my $process_quality_measures = sub {
	my $locus = shift;
	my $individual_q_measures = shift;
	my $num_correct = shift;
	my $num_total = shift;
	
	if(exists $individual_q_measures->{P})
	{
		my $P = $individual_q_measures->{P};
		
		if($num_correct == $num_total)
		{
			push(@{$lls_correct_individuals{$locus}},  $P);
		}	
		else
		{
			push(@{$lls_incorrect_individuals{$locus}},  $P);	
		}
	}
		
};

# compare predictions
my %problem_haplo_counter;
my %problem_haplo_detail;
my %problem_locus_detail;
my %problem_locus_examined;
my %reference_predictions;
my %imputations_predictions;
my %locus_expected_OK;

my $comparisons;
my $compare_problems;

foreach my $locus (@loci)
{
	my $arbitraty_indiv = (keys %reference_data)[0];
	next unless((defined $reference_data{$arbitraty_indiv}{'HLA'.$locus}) or (defined $reference_data{$arbitraty_indiv}{'KIR'.$locus}));
	
	$problem_locus_detail{$locus} = 0;
	my @indivIDs = keys %{$imputed_HLA{$locus}};
	INDIV: foreach my $indivID (@indivIDs)
	{	
		$debug = 0;
		
		my @imputed_hla_values = map { $imputed_HLA{$locus}{$indivID}{$_} } keys %{$imputed_HLA{$locus}{$indivID}};
		my @reference_hla_values;
		
		next INDIV unless($#imputed_hla_values == 1);
		next INDIV unless(exists $reference_data{$indivID});
		
		($reference_data{$indivID}{'HLA'.$locus} xor $reference_data{$indivID}{'KIR'.$locus}) or die;
		
		if($reference_data{$indivID}{'HLA'.$locus})
		{
			@reference_hla_values = split(/\//, $reference_data{$indivID}{'HLA'.$locus});
		}
		else
		{
			@reference_hla_values = split(/\//, $reference_data{$indivID}{'KIR'.$locus});
			
		}
		die Dumper($reference_data{$indivID}, \@reference_hla_values) unless($#reference_hla_values == 1);
		
		@reference_hla_values = grep {! &simpleHLA::is_missing($_)} @reference_hla_values;
		
		if($only_4_dig)
		{
			next unless (&simpleHLA::HLA_is4digit($reference_hla_values[0]) and (($#reference_hla_values == 0) || (&simpleHLA::HLA_is4digit($reference_hla_values[1]))));
		}
		
		$imputed_HLA_Calls{$locus}{sum} += scalar(@imputed_hla_values);		
		@imputed_hla_values = grep {! &simpleHLA::is_missing($_)} @imputed_hla_values;
		$imputed_HLA_Calls{$locus}{called} += scalar(@imputed_hla_values);
		
		if($all_2_dig)
		{
			@reference_hla_values = map {&simpleHLA::HLA_2digit($_)} @reference_hla_values;
			@imputed_hla_values = map {&simpleHLA::HLA_2digit($_)} @imputed_hla_values;
		}
		
		
		if($#imputed_hla_values == 1)
		{
			if($#reference_hla_values > -1)
			{
				if($#reference_hla_values == 0)
				{
					if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatible($locus, $reference_hla_values[0], $imputed_hla_values[1]))
					{
						$comparisons++;
						$problem_locus_examined{$locus}++;
						
						if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
						}
						elsif(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[1]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;						
						}
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);
						
					}	
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect} += 0.5;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect} += 0.5;
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);
					}
					
					if($training_snp_file)
					{
						my $f_a1 = $training_allele_frequencies{$locus}{$reference_hla_values[0]} ? $training_allele_frequencies{$locus}{$reference_hla_values[0]} : 0;
						my $expected_OK = 2*$f_a1*(1-$f_a1)+($f_a1**2);
						die unless($expected_OK <= 1);
						$locus_expected_OK{$locus} += $expected_OK;
					}
				}
				elsif($#reference_hla_values == 1)
				{
					if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]) and &compatible($locus, $reference_hla_values[1], $imputed_hla_values[1]))
					{
						$comparisons += 2;
						$problem_locus_examined{$locus} += 2;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;	
						$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;	
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);	
					}		
					elsif(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[1]) and &compatible($locus, $reference_hla_values[1], $imputed_hla_values[0]))
					{
						$comparisons += 2;
						$problem_locus_examined{$locus} += 2;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
						$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;			
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);					
					}
					else
					{
						if(
							&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatible($locus, $reference_hla_values[1], $imputed_hla_values[1]) or
							&compatible($locus, $reference_hla_values[1], $imputed_hla_values[0]) or &compatible($locus, $reference_hla_values[0], $imputed_hla_values[1])
						)
						{
							$comparisons += 2;
							$compare_problems++;
							$problem_haplo_counter{$indivID}++;
							$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 2)';						
							$problem_locus_detail{$locus}++;
							$problem_locus_examined{$locus} += 2;
							
							if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;							
							}
							elsif(&compatible($locus, $reference_hla_values[1], $imputed_hla_values[1]))
							{
								$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;									
							}
							elsif(&compatible($locus, $reference_hla_values[1], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;						
							}
							elsif(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[1]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
								$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;								
							}
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 2);
						}
						else
						{
							$comparisons += 2;
							$compare_problems += 2;
							$problem_haplo_counter{$indivID} += 2;
							$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (2 problems of 2)';						
							$problem_locus_detail{$locus} += 2;
							$problem_locus_examined{$locus} += 2;
							
							$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;					
							$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;	
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 2);						
						}
					}
					
					if($training_snp_file)
					{
						my $expected_accuracy = 0;
						
						my $f_a1 = $training_allele_frequencies{$locus}{$reference_hla_values[0]} ? $training_allele_frequencies{$locus}{$reference_hla_values[0]} : 0;
						my $f_a2 = $training_allele_frequencies{$locus}{$reference_hla_values[1]} ? $training_allele_frequencies{$locus}{$reference_hla_values[1]} : 0;
						
						if($reference_hla_values[0] ne $reference_hla_values[1])
						{
							$expected_accuracy += 2 * ($f_a1) * ($f_a2);
						}
						else
						{
							$expected_accuracy += ($f_a1) * ($f_a2);
						}
						
						$expected_accuracy *= 2;
						
						die unless($expected_accuracy <= 2);
						
						$locus_expected_OK{$locus} += $expected_accuracy;
					}
										
				}
				else
				{
					die;
				}
			}
		}
		elsif($#imputed_hla_values == 0)
		{
			if($#reference_hla_values > -1)
			{
				if($#reference_hla_values == 0)
				{
					if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]))
					{
						$comparisons++;
						$problem_locus_examined{$locus}++;
						$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);
						
					}	
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);									
					}
					
					if($training_snp_file)
					{
						my $f_a1 = $training_allele_frequencies{$locus}{$reference_hla_values[0]} ? $training_allele_frequencies{$locus}{$reference_hla_values[0]} : 0;
						my $expected_OK = $f_a1;
						die unless($expected_OK <= 1);
						$locus_expected_OK{$locus} += $expected_OK;
					}
					
				}
				elsif($#reference_hla_values == 1)
				{
					if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatible($locus, $reference_hla_values[1], $imputed_hla_values[0]))
					{
						$comparisons += 1;
						$problem_locus_examined{$locus} += 1;
						
						if(&compatible($locus, $reference_hla_values[0], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;						
						}
						elsif(&compatible($locus, $reference_hla_values[1], $imputed_hla_values[0]))
						{
							$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;							
						}
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 1, 1);

					}		
					else
					{
						$comparisons++;
						$compare_problems++;
						$problem_haplo_counter{$indivID}++;
						$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 1)';
						$problem_locus_detail{$locus}++;
						$problem_locus_examined{$locus}++;
						
						$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect} += 0.5;
						$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect} += 0.5;
						$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;	
						
						$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);								
					}
					
					if($training_snp_file)
					{
						my $f_a1 = $training_allele_frequencies{$locus}{$reference_hla_values[0]} ? $training_allele_frequencies{$locus}{$reference_hla_values[0]} : 0;
						my $f_a2 = $training_allele_frequencies{$locus}{$reference_hla_values[1]} ? $training_allele_frequencies{$locus}{$reference_hla_values[1]} : 0;
						
						my $expected_OK = $f_a1+$f_a2;
						die unless($expected_OK <= 1);
						$locus_expected_OK{$locus} += $expected_OK;
					}
				}
				else
				{
					die;
				}
			}		
		}
		
	}
}

print "\n\n EXAMINED: $comparisons -- PROBLEMS: $compare_problems\n\n\n";

open(TMP_OUTPUT, '>', 'temp_validation_output.txt') or die;
print TMP_OUTPUT join("\t", qw/Locus N CallRate Accuracy ExpectationFromGuessing/), "\n";
print "PER-LOCUS SUMMARY:\n";
foreach my $key (sort keys %problem_locus_detail)
{
	my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$key} != 0) ? $problem_locus_detail{$key}/$problem_locus_examined{$key} : 0));
	my $CR = sprintf("%.2f", $imputed_HLA_Calls{$key}{sum} ? ($imputed_HLA_Calls{$key}{called}/$imputed_HLA_Calls{$key}{sum}) : 0);
	my $expectedAcc = -1;
	print "\t", $key, ": ", $problem_locus_detail{$key}, " of ", $problem_locus_examined{$key}, ", accuracy ", $accuracy, "\n";
	print "\t\t Call rate: ", $CR,  "\n";
	if($training_snp_file)
	{
		$expectedAcc = sprintf("%.2f", $problem_locus_examined{$key} ? ($locus_expected_OK{$key}/$problem_locus_examined{$key}) : 0);
		print "\t\t Expected accuracy from guessing:: ", $expectedAcc,  "\n";
	}
	
	my @fields = ($key, $problem_locus_examined{$key}, $CR, $accuracy, $expectedAcc);
	print TMP_OUTPUT join("\t", @fields), "\n";
}

close(TMP_OUTPUT);

#print Dumper($imputations_predictions{A});


print "\n\nPER-HAPLOTYPE SUMMARY:\n";
foreach my $key (keys %problem_haplo_counter)
{
	#print "\t", $key, ": ", $problem_haplo_counter{$key}, "\n";
}


# print "\n\nCHECKED FIELDS:\n\n".Dumper(\%field_to_check)."\n\n";

# print "\n\n";

# print "HLA in REFERENCE en detail:\n";
# print Dumper($reference_predictions{DRB}), "\n\n";

# print "HLA in PREDICTIONS in detail\n";
# print Dumper($imputations_predictions{DRB});

 foreach my $key (keys %problem_haplo_detail)
 {
	 foreach my $fieldName (keys %{$problem_haplo_detail{$key}})
	 {
		#next unless(($fieldName eq 'A') or ($fieldName eq 'B'));
		#next unless($fieldName eq 'DQB');
		#print $key, "\t", $fieldName, ": ", $problem_haplo_detail{$key}{$fieldName}, "\n";
	 }
 }

if(1 == 0)
{
	print "\n\nLikelihood analysis:\n";
	foreach my $locus (sort keys %problem_locus_detail)
	{
		my $mean_correct;
		my $mean_incorrect;
		my $variance_correct;
		my $variance_incorrect;
		
		foreach my $v (@{$lls_correct_individuals{$locus}})
		{
			$mean_correct += $v;
		}
		$mean_correct = $mean_correct / scalar(@{$lls_correct_individuals{$locus}});
		
		foreach my $v (@{$lls_incorrect_individuals{$locus}})
		{
			$mean_incorrect += $v;
		}
		$mean_incorrect = $mean_incorrect / scalar(@{$lls_incorrect_individuals{$locus}});
			
		foreach my $v (@{$lls_correct_individuals{$locus}})
		{
			$variance_correct += ($v - $mean_correct)**2;
		}
		$variance_correct = $variance_correct / scalar(@{$lls_correct_individuals{$locus}});
		
		foreach my $v (@{$lls_incorrect_individuals{$locus}})
		{
			$variance_incorrect += ($v - $mean_incorrect)**2
		}
		$variance_incorrect = $variance_incorrect / scalar(@{$lls_incorrect_individuals{$locus}});
		
		$lls_correct_individuals{$locus} = [sort {$a <=> $b} @{$lls_correct_individuals{$locus}}];
		$lls_incorrect_individuals{$locus} = [sort {$a <=> $b} @{$lls_incorrect_individuals{$locus}}];
		
		my $lower_correct = int((scalar(@{$lls_correct_individuals{$locus}})-1)*0.1);
		my $upper_correct = int((scalar(@{$lls_correct_individuals{$locus}})-1)*0.9);
		my $lower_incorrect = int((scalar(@{$lls_incorrect_individuals{$locus}})-1)*0.1);
		my $upper_incorrect = int((scalar(@{$lls_incorrect_individuals{$locus}})-1)*0.9);
				
		my $correct_lower_value = $lls_correct_individuals{$locus}[$lower_correct];
		my $correct_upper_value = $lls_correct_individuals{$locus}[$upper_correct];
		
		my $incorrect_lower_value = $lls_incorrect_individuals{$locus}[$lower_incorrect];
		my $incorrect_upper_value = $lls_incorrect_individuals{$locus}[$upper_incorrect];
			
		print "\t$locus\n";
		
		print "\t\t ", scalar(@{$lls_correct_individuals{$locus}}), " correct, mean likelihood $mean_correct\n";
		print "\t\t\t sd ", sqrt($variance_correct), "\n";
		print "\t\t\t from $correct_lower_value to $correct_upper_value \n";	
		print "\t\t ", scalar(@{$lls_incorrect_individuals{$locus}}), " false, mean likelihood $mean_incorrect\n";
		print "\t\t\t sd ", sqrt($variance_incorrect), "\n";
		print "\t\t\t from $incorrect_lower_value to $incorrect_upper_value \n\n";	
	}
}

sub compatible
{
	my $locus = shift;	
	my $four_or_two_dig = shift;
	my $four_dig = shift;
	
	if($all_2_dig)
	{
		my $four_dig_reduced = &simpleHLA::HLA_2digit($four_dig);
		# print &simpleHLA::HLA_2digit($four_or_two_dig), "\t", $four_dig_reduced, "\t", (&simpleHLA::HLA_2digit($four_or_two_dig) eq $four_dig_reduced) ? 1 : 0, "\n" if ($locus eq 'DRB');
		return (&simpleHLA::HLA_2digit($four_or_two_dig) eq $four_dig_reduced) ? 1 : 0;		
	}
	else
	{
		if(&simpleHLA::HLA_is4digit($four_or_two_dig))
		{
			return (&simpleHLA::HLA_4digit($four_or_two_dig) eq &simpleHLA::HLA_4digit($four_dig)) ? 1 : 0;
		}
		else
		{
			my $four_dig_reduced = &simpleHLA::HLA_2digit($four_dig);
			return (&simpleHLA::HLA_2digit($four_or_two_dig) eq $four_dig_reduced) ? 1 : 0;
		}
	}
}	

