#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;
use simpleHLA;
use File::Basename;

my $kMer_size = 55;  


# my @testCases = (
	# [[qw/A A/], [qw/A A/]],
	# [[qw/? A/], [qw/A A/]],
	# [[qw/A ?/], [qw/A A/]],
	# [[qw/A T/], [qw/A A/]],
	# [[qw/A A/], [qw/T A/]],
	# [[qw/A C/], [qw/G T/]],
	# [[qw/A C/], [qw/? T/]],
	# [[qw/? C/], [qw/G T/]],
	# [[qw/? ?/], [qw/G T/]],
	# [[qw/? T/], [qw/? T/]],
	# [[qw/? C/], [qw/? T/]],
	
# );
# foreach my $testCase (@testCases)
# {
	# print join(' vs ', map {join('/', @$_)} @$testCase), "   ", join(' ', compatibleStringAlleles($testCase->[0], $testCase->[1])), "\n";
# }
# exit;

# input parameters

my $graph = 'hla';   
my $sampleIDs = '';
my $BAMs = '';
my $actions;
my $trueHLA;
my $trueHaplotypes;
#my $validation_round = 'R1';
my $T = 0;
my $minCoverage = 0;
my $all_2_dig = 0;
my $only_4_dig = 1;
my $HiSeq250bp = 0;
my $fastExtraction = 0;

my $fromPHLAT = 0;

my $referenceGenome;

my @loci_for_check = qw/A B C DQA1 DQB1 DRB1/;

GetOptions ('graph:s' => \$graph,
 'sampleIDs:s' => \$sampleIDs, 
 'BAMs:s' => \$BAMs, 
 'actions:s' => \$actions, 
 'trueHLA:s' => \$trueHLA,
 'trueHaplotypes:s' => \$trueHaplotypes, 
 'referenceGenome:s' => \$referenceGenome, 
 #'validation_round:s' => \$validation_round,
 'T:s' => \$T,
 'minCoverage:s' => \$minCoverage,
 'all_2_dig:s' => \$all_2_dig,
 'only_4_dig:s' => \$only_4_dig,
 'HiSeq250bp:s' => \$HiSeq250bp, 
 'fastExtraction:s' => \$fastExtraction, 
 'fromPHLAT:s' => \$fromPHLAT,
);         

if($minCoverage)
{
	print "Minimum coverage threshold in place: $minCoverage\n";
}

if($fastExtraction)
{
	$HiSeq250bp = 1;
}

my $genome_graph_file = qq(../tmp2/GS_nextGen/hla/derived/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph_k25.ctx);
unless(-e $genome_graph_file)
{
	die "Please set variable \$genome_graph_file to an existing file - the current value $genome_graph_file is not accessible.";
}

my $expected_kMer_file = qq(../tmp2/GS_nextGen/${graph}/requiredkMers_graph.txt.kmers_25);
unless(-e $expected_kMer_file)
{
	die "Please set variable \$expected_kMer_file to an existing file - the current value $expected_kMer_file is not accessible.";
}

my $exon_folder = qq(../tmp2/GS_nextGen/${graph}/);
unless(-e $expected_kMer_file)
{
	die "Please provide a kMerified graph -- exepcted file $expected_kMer_file not there!";
}

my $normal_bin = qq(../bin/MHC-PRG);
my $cluster3_bin = qq(../bin/MHC-PRG);
my $use_bin = ((hostname() =~ /cluster3/) or (hostname() =~ /^comp[AB]\d+$/)) ? $cluster3_bin : $normal_bin;
unless(-e $use_bin)
{
	die "Cannot find expected binary: $use_bin";
}
	
my @BAMs = split(/,/, $BAMs);
my @sampleIDs = split(/,/, $sampleIDs);
if(@sampleIDs)
{
	foreach my $sampleID (@sampleIDs)
	{
		unless($sampleID =~ /^[\w]+$/)
		{
			die "Please provide only sample IDs with normal characters.";
		}
	}
}	

my $sample_IDs_abbr;
if($sampleIDs =~ /^allSimulations(_\w+)?/)
{
	my $addFilter = $1;
	my @dirs;
	if($addFilter)
	{
		@dirs = grep {$_ =~ /I\d+_simulations${addFilter}/} grep {-d $_} glob('../tmp/hla/*');
	}
	else
	{
		@dirs = grep {$_ =~ /I\d+_simulations/} grep {-d $_} glob('../tmp/hla/*');
	}
	
	@sampleIDs = map {die "Can't parse $_" unless($_ =~ /tmp\/hla\/(.+)/); $1} @dirs;
	
	if($sampleIDs =~ /^all_simulations_I(\d+)/i)
	{
		my $iteration = $1;
		@sampleIDs = grep {$_ =~ /^I${iteration}_/i} @sampleIDs;
	}
	
	my $debug = 1;
	if($debug)
	{
		@sampleIDs = grep {die unless($_ =~ /sample(\d+)$/); ($1 < 30)} @sampleIDs;
	}
	
	$sample_IDs_abbr = $sampleIDs;
}
elsif($sampleIDs =~ /^all/)
{
	my @dirs = grep {$_ !~ /simulations/} grep {-d $_} glob('../tmp/hla/*');
	@sampleIDs = map {die "Can't parse $_" unless($_ =~ /tmp\/hla\/(.+)/); $1} @dirs;
	
	if($sampleIDs =~ /^all_I(\d+)/i)
	{
		my $iteration = $1;
		@sampleIDs = grep {$_ =~ /^I${iteration}_/i} @sampleIDs;
	}
	
	if($actions =~ /v|w/)
	{
		@sampleIDs = @sampleIDs[5 .. $#sampleIDs];
		@sampleIDs = grep {$_ !~ /AA02O9Q_Z2/} @sampleIDs;
		@sampleIDs = grep {$_ !~ /AA02O9R/} @sampleIDs;
		
		warn "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!\n\nRemove first five and NA12878 and high-coverage sample for validation!\n\n!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
		
		die unless(($actions eq 'v') or ($actions eq 'w'));
	}
	
	$sample_IDs_abbr = $sampleIDs;
}
else
{
	$sample_IDs_abbr = join('_', @sampleIDs);
	if(length($sample_IDs_abbr) > 50)
	{
		$sample_IDs_abbr = substr($sample_IDs_abbr, 0, 50);
	}
}

if(scalar(@sampleIDs) > 5)
{
	#@sampleIDs = @sampleIDs[0 .. 4];
	#warn "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!\n\nLimited samples!\n\n!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
}

#@sampleIDs = $sampleIDs[0]; # todo remove
#warn "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!\n\nLimited samples:\n".join("\n", @sampleIDs)."\n\n!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";

if($actions =~ /p/)
{
	unless(@BAMs)
	{
		die "Please provide --BAMs for positive filtering";
	}
	unless($#BAMs == $#sampleIDs)
	{
		die "Please provide an equal number of --BAMs and --sampleIDs";
	}
	
	for(my $bI = 0; $bI <= $#BAMs; $bI++)
	{
		my $BAM = $BAMs[$bI];
		my $sampleID = $sampleIDs[$bI];
		
		unless(-e $BAM)
		{
			die "Specified BAM $BAM (in --BAMs) does not exist!\n";
		}
		
		my $output_file = '../tmp/hla/'.$sampleID.'/reads.p';
		unless(-e '../tmp/hla/'.$sampleID)
		{
			mkdir('../tmp/hla/'.$sampleID) or die "Cannot mkdir ".'../tmp/hla/'.$sampleID;
		}
		
		my $command = qq($use_bin domode filterReads --input_BAM $BAM --positiveFilter $expected_kMer_file --output_FASTQ $output_file );
		
		if($referenceGenome)
		{
			$command .= qq( --referenceGenome $referenceGenome);
		}	
		
		if($HiSeq250bp)
		{
			$command .= qq( --HiSeq250bp);
		}
		
		print "Now executing command:\n$command\n\n";
		
		system($command);
	}
}

if($actions =~ /n/)
{
	unless(scalar(@sampleIDs))
	{
		die "Please provide some --sampleIDs for negative filtering.";
	}
		
	my @fastQ_files;
	my @output_files;
	foreach my $sampleID (@sampleIDs)
	{
		my $fastQ_file = '../tmp/hla/'.$sampleID.'/reads.p';
		my $fastQ_file_1 = $fastQ_file.'_1';
		my $fastQ_file_2 = $fastQ_file.'_2';
		unless(-e $fastQ_file_1)
		{
			die "Expected file $fastQ_file_1 not found";
		}
		unless(-e $fastQ_file_2)
		{
			die "Expected file $fastQ_file_2 not found";
		}		
		my $output_file = '../tmp/hla/'.$sampleID.'/reads.p.n';
		
		push(@fastQ_files, $fastQ_file);
		push(@output_files, $output_file);
	}
	
	my $fastQ_files = join(',', @fastQ_files);
	my $output_files = join(',', @output_files);
	
	my $command = qq($use_bin domode filterReads --input_FASTQ $fastQ_files --negativeFilter $genome_graph_file --output_FASTQ $output_files --negativePreserveUnique --uniqueness_base ${expected_kMer_file} --uniqueness_subtract ${genome_graph_file});
	
	print "Now executing command:\n$command\n\n";
	
	system($command);
}


if($actions =~ /a/)
{
	unless(scalar(@sampleIDs))
	{
		die "Please provide some --sampleIDs for alignment.";
	}
		
	my @fastQ_files;
	foreach my $sampleID (@sampleIDs)
	{
		my $fastQ_file = '../tmp/hla/'.$sampleID.'/reads.p.n';
		my $fastQ_file_1 = $fastQ_file.'_1';
		my $fastQ_file_2 = $fastQ_file.'_2';
		unless(-e $fastQ_file_1)
		{
			die "Expected file $fastQ_file_1 not found";
		}
		unless(-e $fastQ_file_2)
		{
			die "Expected file $fastQ_file_2 not found";
		}		
		my $output_file = '../tmp/hla/'.$sampleID.'/reads.p.n';
		
		push(@fastQ_files, $fastQ_file);
	}
	
	my $fastQ_files = join(',', @fastQ_files);
	
	my $pseudoReferenceGenome = qq(../tmp2/GS_nextGen/${graph}/pseudoReferenceGenome.txt);
	unless(-e $pseudoReferenceGenome)
	{
		die "Pseudo-reference file $pseudoReferenceGenome not existing.";
	}
	my $command = qq($use_bin domode alignShortReadsToHLAGraph --input_FASTQ $fastQ_files --graphDir ../tmp2/GS_nextGen/${graph} --referenceGenome ${pseudoReferenceGenome});
	
	print "Now executing command:\n$command\n\n";
	
	system($command);
}


if($actions =~ /i/)
{
	unless(scalar(@sampleIDs))
	{
		die "Please provide some --sampleIDs for HLA type inference.";
	}
		
	my @aligned_files;
	my @stdout_files;
	foreach my $sampleID (@sampleIDs)
	{
		my $aligned_file = '../tmp/hla/'.$sampleID.'/reads.p.n.aligned';
		unless(-e $aligned_file)
		{
			die "Expected file $aligned_file not found";
		}
	
		push(@aligned_files, $aligned_file);
		
		my $stdout_file = '../tmp/hla/'.$sampleID.'/inference.stdout';
		push(@stdout_files, $stdout_file);
		
		foreach my $validation_round (qw/R1 R2/)
		{
			my $bestguess_file = '../tmp/hla/'.$sampleID.'/'.$validation_round.'_bestguess.txt';
			if(-e $bestguess_file)
			{
				warn "Delete existing best-guess file $bestguess_file";
				unlink($bestguess_file) or die "Cannot delete $bestguess_file";
			}
		}
	}
		
	for(my $sI = 0; $sI <= $#aligned_files; $sI++)
	{
		my $sampleID = $sampleIDs[$sI];
		my $aligned_file = $aligned_files[$sI];
		my $stdout_file = $stdout_files[$sI];
		
		
		my ($aligned_file_name, $aligned_file_path) = fileparse($aligned_file);
					
		my $command = qq($use_bin domode HLATypeInference --input_alignedReads $aligned_file --graphDir ../tmp2/GS_nextGen/${graph} --sampleID $sampleID > $stdout_file);
	
		print "Now executing command:\n$command\n\n";
		
		my $ret = system($command);	
		
		unless($ret == 0)
		{
			die "When executing $command, got return code $ret";
		}
		
		my $expected_bestguess_file = $aligned_file_path . '/' . 'R1_bestguess.txt';
		unless(-e $expected_bestguess_file)
		{
			die "File $expected_bestguess_file not existing!";			
		}
		
		my %l_counter;
		open(F, '<', $expected_bestguess_file) or die "Cannot open $expected_bestguess_file";
		<F>;
		while(<F>)
		{
			my $l = $_;
			die unless($l =~ /^(\w+)\t/);
			$l_counter{$1}++;
		}
		close(F);

		foreach my $locus (@loci_for_check)
		{
			unless($l_counter{$locus} == 2)
			{
				die Dumper("Wrong bestguess count", $locus, $expected_bestguess_file, \%l_counter);
			}
			
			my $bestguess_haplotypes_file =  $aligned_file_path . '/' . 'R2_haplotypes_bestguess_' . $locus . '.txt';
			unless(-e $bestguess_haplotypes_file)
			{
				#die "Expected best-guess haplotypes file cannot be found : ". $bestguess_haplotypes_file;
			}	
		}
	}
}  

if($actions =~ /v/)
{
	my $validation_round = 'R1';
	die "Please specify --trueHLA for validation" unless($trueHLA);
			
	# read reference dataset
	my %reference_data;
	open(REFERENCE, "<", $trueHLA) or die "Cannot open $trueHLA";
	my $headerLine = <REFERENCE>;
	chomp($headerLine);
	$headerLine =~ s/\n//g;
	$headerLine =~ s/\r//g;
	my @header_fields = split(/[\t ]/, $headerLine);
	@header_fields = map {if($_ =~ /HLAD((QA)|(QB)|(RB))$/){$_ .= '1';} $_} @header_fields;	
	while(<REFERENCE>)
	{
		my $line = $_;
		chomp($line);
		
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		
		my @fields = split(/[\t ]/, $line);
		my %line = (mesh @header_fields, @fields);
		
		my $primary_key = $line{'IndividualID'};
		$reference_data{$primary_key} = \%line;
	}
	close(REFERENCE);
	
	# die Dumper(\%reference_data);
	
	my %imputed_HLA;
	my %imputed_HLA_Q;
	my %imputed_HLA_avgCoverage;
	my %imputed_HLA_lowCoverage;
	my %imputed_HLA_minCoverage;

	my %sample_noI_toI;
	
	my $summary_file = 'temp/summary_' . $sample_IDs_abbr . '.txt';	
	open(SUMMARY, '>', $summary_file) or die "Cannot open $summary_file";
	print SUMMARY $sample_IDs_abbr, "\n";
	print SUMMARY "\t", join("\t", qw/Locus N CallRate Accuracy/), "\n";
	foreach my $sampleID (@sampleIDs)
	{
		my $sampleID_noI = $sampleID;
		$sampleID_noI =~ s/^I\d+_//g;
		
		
		my $bestGuess_file;
		
		if($fromPHLAT)
		{
			$bestGuess_file = '/gpfs1/well/gsk_hla/PHLAT/'.$sampleID.'/'.$validation_round.'_bestguess.txt';	
			unless(-e $bestGuess_file)
			{
				warn "Best-guess file $bestGuess_file not existing";
				next;
			}		
			
		}
		else
		{
			$bestGuess_file = '../tmp/hla/'.$sampleID.'/'.$validation_round.'_bestguess.txt';
			unless(-e $bestGuess_file)
			{
				warn "Best-guess file $bestGuess_file not existing";
				next;
			}		
		}
		  
		open(BESTGUESS, '<', $bestGuess_file) or die "Cannot open $bestGuess_file";
		my $bestguess_header_line = <BESTGUESS>;
		chomp($bestguess_header_line);
		my @bestguess_header_files = split(/\t/, $bestguess_header_line);
		while(<BESTGUESS>)
		{
			my $line = $_;
			chomp($line);
			my @line_fields = split(/\t/, $line);
			my %line_hash = (mesh @bestguess_header_files, @line_fields);
			
			my $Q = $line_hash{'Q1'};
			$imputed_HLA_Q{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = $Q;
			if($Q < $T)
			{
				$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = '??:??';			
			}
			else
			{
				die unless(defined $line_hash{'Allele'});
				$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = $line_hash{'Allele'};			
			}
			
			
			if($line_hash{'Chromosome'} eq '1')
			{
				$imputed_HLA_avgCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'AverageCoverage'};
				$imputed_HLA_lowCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'CoverageFirstDecile'};
				$imputed_HLA_minCoverage{$line_hash{'Locus'}}{$sampleID_noI} = $line_hash{'MinimumCoverage'};
				
				if($minCoverage and ($line_hash{'MinimumCoverage'} < $minCoverage))
				{
					$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = '??:??';								
				}
			}
		}	
		close(BESTGUESS);
		
		die if(exists $sample_noI_toI{$sampleID_noI});
		$sample_noI_toI{$sampleID_noI} = $sampleID;
	}
		
	my $debug = 0;
	my $comparisons = 0;
	my $compare_problems = 0;
	my %locus_avgCoverages;
	my %locus_lowCoverages;
	my %locus_minCoverages;
	
	my %problem_locus_detail;
	my %problem_locus_examined;
	my %problem_haplo_counter;
	my %problem_haplo_detail;
	my %imputations_predictions;
	my %reference_predictions;
	my %imputed_HLA_Calls;
	my %quality_measures; # not used
	my $pileup_href = {};
	
	my @loci = sort keys %imputed_HLA;
	
	my $process_quality_measures = sub {};
	
	my $PP_to_basket = sub {
		my $PP = shift;
		die unless(($PP >= 0) && ($PP <= 1));
		my $basket = int($PP * 10);
		$basket = 9 if($basket == 10);
		return $basket;			
	};
		
	foreach my $locus (@loci)
	{
		my $arbitraty_indiv = (keys %reference_data)[0];
		next unless((defined $reference_data{$arbitraty_indiv}{'HLA'.$locus}));
		
		my %calibration_baskets;
		my %coverage_over_samples;
		my %coverage_over_samples_individualValues;
		my $coverage_over_samples_nSamples = 0;
		
		my $add_to_calibration_basket = sub {
			my $str_correct = shift;
			my $PP = shift;
			my $weight = shift;
			
			die unless(($str_correct eq 'correct') or ($str_correct eq 'incorrect'));
			die unless(defined $PP);
			die unless(defined $weight);
			
			push(@{$calibration_baskets{$PP_to_basket->($PP)}{$str_correct}}, {PP => $PP, weight => $weight});
		};	
		
		$problem_locus_examined{$locus} = 0;
		$problem_locus_detail{$locus} = 0;
		my @indivIDs = keys %{$imputed_HLA{$locus}};
		my $indivI_processedForEvaluation = -1;
		INDIV: foreach my $indivID (@indivIDs)
		{	
			$| = 1;
						
			$debug = 0;
			
			my @imputed_hla_values = map { $imputed_HLA{$locus}{$indivID}{$_} } keys %{$imputed_HLA{$locus}{$indivID}};
			my @imputed_hla_values_q = map { my $r = $imputed_HLA_Q{$locus}{$indivID}{$_}; die unless(defined $r); $r } keys %{$imputed_HLA{$locus}{$indivID}};
			
			die "Undefined HLA ".join(', ', @imputed_hla_values) unless(scalar(grep {defined $_} @imputed_hla_values) == scalar(@imputed_hla_values));
					
			my @reference_hla_values;
			
			next INDIV unless($#imputed_hla_values == 1);
			unless(exists $reference_data{$indivID})
			{
				warn "No reference data for $locus $indivID";
			}
			next INDIV unless(exists $reference_data{$indivID});
			
			$reference_data{$indivID}{'HLA'.$locus} or die;
			
			if($reference_data{$indivID}{'HLA'.$locus})
			{
				@reference_hla_values = split(/\//, $reference_data{$indivID}{'HLA'.$locus});
			}

			die Dumper($reference_data{$indivID}, \@reference_hla_values) unless($#reference_hla_values == 1);				
						
			die "Undefined HLA ".join(', ', @reference_hla_values) unless(scalar(grep {defined $_} @reference_hla_values) == scalar(@reference_hla_values));			
			@reference_hla_values = grep {! &simpleHLA::is_missing($_)} @reference_hla_values;

			next if($#reference_hla_values == -1);
						
			if($only_4_dig)
			{
				next unless (&simpleHLA::HLA_is4digit($reference_hla_values[0]) and (($#reference_hla_values == 0) || (&simpleHLA::HLA_is4digit($reference_hla_values[1]))));
			}
					
			$imputed_HLA_Calls{$locus}{sum} += scalar(@imputed_hla_values);		
			$indivI_processedForEvaluation++;
			
			my @imputed_present = map {(! &simpleHLA::is_missing($_)) ? 1 : 0} @imputed_hla_values;
			my @imputed_hla_values_q_new;
			for(my $i = 0; $i <= $#imputed_hla_values; $i++)
			{
				my $Q = $imputed_hla_values_q[$i];
				die unless(defined $Q);
				die unless(defined $imputed_present[$i]);
				if($imputed_present[$i])
				{
					push(@imputed_hla_values_q_new, $Q);
				}
			}
			
			@imputed_hla_values = grep {! &simpleHLA::is_missing($_)} @imputed_hla_values;
			@imputed_hla_values_q = @imputed_hla_values_q_new;
			die unless($#imputed_hla_values == $#imputed_hla_values_q);
			
			$imputed_HLA_Calls{$locus}{called} += scalar(@imputed_hla_values);
			
			if($all_2_dig)
			{
				@reference_hla_values = map {&simpleHLA::HLA_2digit($_)} @reference_hla_values;
				@imputed_hla_values = map {join(';', map {&simpleHLA::HLA_2digit($_)} split(/;/, $_))} @imputed_hla_values;
			}
						
			# print Dumper(\@reference_hla_values, @imputed_hla_values), "\n";
		
			my $comparisons_before = $comparisons;
			my $problem_locus_detail_before = $problem_locus_detail{$locus};
			
			if($#imputed_hla_values == 1)
			{
				if($#reference_hla_values > -1)
				{
					if($#reference_hla_values == 0)
					{
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
						{
							$comparisons++;
							$problem_locus_examined{$locus}++;
							
							if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
								
								$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);														
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;

								$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);								
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
							
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 0.5);
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 0.5);
							
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);
						}
					}
					elsif($#reference_hla_values == 1)
					{
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) and &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]))
						{
							$comparisons += 2;
							$problem_locus_examined{$locus} += 2;
							
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;	
							$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;	
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);	
							
							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
							$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
							
						}		
						elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]) and &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
						{
							$comparisons += 2;
							$problem_locus_examined{$locus} += 2;
							
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
							$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;			
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 2, 2);		

							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
							$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
							
						}
						else
						{
							if(
								&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]) or
								&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1])
							)
							{
								$comparisons += 2;
								$compare_problems++;
								$problem_haplo_counter{$indivID}++;
								$problem_haplo_detail{$indivID}{$locus} = 'Reference: '.join('/', @reference_hla_values).' VS Imputed: '.join('/', @imputed_hla_values).'   (1 problem of 2)';						
								$problem_locus_detail{$locus}++;
								$problem_locus_examined{$locus} += 2;
								
								if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
								{
									$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;		

									$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
									$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);									
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]))
								{
									$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;		
									
									$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
									$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);	
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
								{
									$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;	
									
									$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
									$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);	
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
								{
									$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[1]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;		

									$add_to_calibration_basket->('correct', $imputed_hla_values_q[1], 1);
									$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);										
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
								
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);
								$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[1], 1);									
								
								$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 2);						
							}
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
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
						{
							$comparisons++;
							$problem_locus_examined{$locus}++;
							$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
							$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;
							
							$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
								
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
							
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);
							
							$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);									
						}					
					}
					elsif($#reference_hla_values == 1)
					{
						if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]) or &compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
						{
							$comparisons += 1;
							$problem_locus_examined{$locus} += 1;
							
							if(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[0]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;

								$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);
								
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
							{
								$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
								$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;	

								$add_to_calibration_basket->('correct', $imputed_hla_values_q[0], 1);								
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
							
							$add_to_calibration_basket->('incorrect', $imputed_hla_values_q[0], 1);

							$process_quality_measures->($locus, $quality_measures{$indivID}, 0, 1);								
						}
					}
					else
					{
						die;
					}
				}		
			}
			
			my $thisIndiv_problems = $problem_locus_detail{$locus} - $problem_locus_detail_before;
				
			my $avgCoverage = $imputed_HLA_avgCoverage{$locus}{$indivID};
			my $lowCoverage = $imputed_HLA_lowCoverage{$locus}{$indivID};
			my $minCoverage = $imputed_HLA_minCoverage{$locus}{$indivID};

			# average coverages
			if(($thisIndiv_problems > 0))
			{
				push(@{$locus_avgCoverages{$locus}{problems}}, $avgCoverage);
				push(@{$locus_lowCoverages{$locus}{problems}}, $lowCoverage);
				push(@{$locus_minCoverages{$locus}{problems}}, $minCoverage);
			}
			else
			{
				push(@{$locus_avgCoverages{$locus}{ok}}, $avgCoverage);
				push(@{$locus_lowCoverages{$locus}{ok}}, $lowCoverage);
				push(@{$locus_minCoverages{$locus}{ok}}, $minCoverage);
				
			}
			
			# print "\t", $thisIndiv_problems, "\n";
			
			# just for debugging - deactivated
			if($thisIndiv_problems == 0)
			{
				if($locus eq 'A')
				{
					# print join(' vs ', join('/', @reference_hla_values), join('/', @imputed_hla_values)), "\n";
				}	
			}
			
			my $thisIndiv_comparions = $comparisons - $comparisons_before;
			my $thisIndiv_OK = $thisIndiv_comparions - $thisIndiv_problems;
			
			my $indivID_withI = $sample_noI_toI{$indivID};
			die unless(defined $indivID_withI);				
			my $pileup_file = qq(../tmp/hla/$indivID_withI/${validation_round}_pileup_${locus}.txt);
				
			my $coverages_href = load_coverages_from_pileup($pileup_file);

			my @k_coverages_existing;
			if($indivI_processedForEvaluation > 0)
			{
				foreach my $exon (keys %coverage_over_samples)
				{
					foreach my $exonPos (keys %{$coverage_over_samples{$exon}})
					{
						push(@k_coverages_existing, $exon . '-/-' . $exonPos);
					}
				}					
			}
			

			my @k_coverages_thisSample;
			foreach my $exon (keys %$coverages_href)
			{
				foreach my $exonPos (keys %{$coverages_href->{$exon}})
				{
					$coverage_over_samples{$exon}{$exonPos} += $coverages_href->{$exon}{$exonPos};
					push(@{$coverage_over_samples_individualValues{$exon}{$exonPos}}, $coverages_href->{$exon}{$exonPos});
					push(@k_coverages_thisSample, $exon . '-/-' . $exonPos);
				}
			}	
			$coverage_over_samples_nSamples++;
			
			if($indivI_processedForEvaluation > 0)
			{
				my ($n_shared, $aref_l1_excl, $aref_l2_excl) = list_comparison(\@k_coverages_existing, \@k_coverages_thisSample);
				if(($#{$aref_l1_excl} == -1) and ($#{$aref_l2_excl} == -1))
				{	
					die unless($n_shared == scalar(@k_coverages_existing));
					die unless($n_shared == scalar(@k_coverages_thisSample));					
				}
				else
				{
					die Dumper("There is a problem with per-exon coverage numbers.", $indivI_processedForEvaluation, scalar(@k_coverages_existing), scalar(@k_coverages_thisSample), scalar(@$aref_l1_excl), scalar(@$aref_l2_excl));
				}
			}

			if(($thisIndiv_problems > 0) and (not $all_2_dig) and not ($fromPHLAT))
			{
				my %readIDs;
				
				unless(-e $pileup_file)
				{
					die "Can't find pileup file $pileup_file";
				}
				
				load_pileup($pileup_href, $pileup_file, $indivID_withI);

			
				my $output_fn = '../tmp/hla_validation/pileup_'.$validation_round.'_'.$indivID_withI.'_'.$locus.'.txt';
				open(my $output_fh, '>', $output_fn) or die "Cannot open $output_fn";
				print $output_fh join("\t", $indivID_withI, $locus, $thisIndiv_OK), "\n";
				
				unless(scalar(@imputed_hla_values) == 2)
				{
					warn "Can't produce pileup for $locus / $indivID";
					next;
				}
				my $inferred = \@imputed_hla_values;
				
				my $truth = \@reference_hla_values;
				
				print {$output_fh} $inferred->[0], "\t\t\t", $truth->[0], "\n";
				print {$output_fh} $inferred->[1], "\t\t\t", $truth->[1], "\n\n";
				
				print {$output_fh} join("\t", "Inferred", "", "", "Apparently true", ""), "\n";
				
				my @exons = print_which_exons($locus);				
				my @inferred_trimmed = twoClusterAlleles($inferred);
				
				foreach my $exon (@exons)
				{
					my $file = find_exon_file($locus, $exon);
					my $sequences = read_exon_sequences($file);
										
					my @validated_extended = twoValidationAlleles_2_proper_names($truth, $locus, $sequences);
								
					my $oneAllele = (keys %$sequences)[0];
					my $length = length($sequences->{$oneAllele});
										
					# print "-", $sequences->{$oneAllele}, "-\n\n";
					
					die unless(scalar(grep {$sequences->{$_}} @validated_extended) == 2);
					die unless(scalar(grep {$sequences->{$_}} @inferred_trimmed) == 2);
					
					print {$output_fh} join("\t",
						$inferred_trimmed[0].' Exon '.$exon,
						$inferred_trimmed[1].' Exon '.$exon,
						"",
						$validated_extended[0].' Exon '.$exon,
						$validated_extended[1].' Exon '.$exon,
						), "\n";
					
					# print Dumper($locus, [keys %{$pileup_href->{$indivID_withI}}]);

					for(my $i = 0; $i < $length; $i++)
					{
						my @chars_for_print;
						
						
						push(@chars_for_print, map {substr($sequences->{$_}, $i, 1)} @inferred_trimmed);
						push(@chars_for_print, '');
						push(@chars_for_print, map {substr($sequences->{$_}, $i, 1)} @validated_extended);
						
						if($chars_for_print[0] ne $chars_for_print[1])
						{
							$chars_for_print[2] = '!';
						}
						
						if($chars_for_print[0] ne $chars_for_print[3])
						{
							$chars_for_print[5] = '!';
						}
						else
						{
							$chars_for_print[5] = '';
						}
						
						if($chars_for_print[1] ne $chars_for_print[4])
						{
							$chars_for_print[6] = '!';
						}			
						else
						{
							$chars_for_print[6] = '';				
						}
						
				
						if(exists $pileup_href->{$indivID_withI}{'HLA'.$locus})
						{
							my $pileUpString = $pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2][$i];
							die unless(defined $pileup_href->{$indivID_withI});
							die unless(defined $pileup_href->{$indivID_withI}{'HLA'.$locus});
							die unless(defined $pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2]);
							# die "Problem with pileup for $locus / $indivID / $exon / $i " unless(defined $pileUpString);
							# next unless(defined $pileUpString);
							
							
							# die Dumper("Pileup too short", $length, scalar(@{$pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2]})) unless(defined $pileUpString);
							unless(($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!'))
							{
								# $pileUpString =~ s/\[.+?\]//g;
							}
							
							my $printPileUpDetails = (1 or (($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!')));
			
							my $getShortRead = sub {
								my $read = shift;
								my $rE = shift;
								
								my @readIDs = split(/ /, $read);
								die unless($#readIDs == 1);
								
								
								
								my $rI; 
								if(exists $readIDs{$readIDs[0]})
								{
									$rI = $readIDs{$readIDs[0]};
								}
								else
								{
									my $nR = scalar(keys %readIDs);
									$nR++;
									$rI = "Read${nR}X";
									$readIDs{$readIDs[0]} = $rI;
									$readIDs{$readIDs[1]} = $rI;
								}
								
								if($printPileUpDetails)
								{
									return $rI.$rE;
								}
								else
								{
									return $rE;
								}
							};
							
							$pileUpString =~ s/(\@\@.+?)(\](\,|$))/$getShortRead->($1, $2);/ge;
							
							if(($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!'))
							{
								push(@chars_for_print, $pileUpString);
							}
						}
						print {$output_fh} join("\t", @chars_for_print), "\n";
					}	
				}
				
				foreach my $readID (keys %readIDs)
				{
					print {$output_fh} $readIDs{$readID}, "\t", $readID, "\n";
				}
				close($output_fh);					
			}
		}
		
		if($T == 0)
		{
			my $totalAlleles = 0;
			
			my $calibration_file = 'temp/calibration_' . $locus . '_' . $sample_IDs_abbr . '.txt';	
			open(CALIBRATION, ">", $calibration_file) or die "Cannot open $calibration_file";
			print CALIBRATION join("\t", qw/Bin MeanPP PercCorrect NCorrect NIncorrect/), "\n";
			for(my $i = 0; $i <= 9; $i++)
			{
				my $meanPP = 0;
				my $percCorrect = 0;
				
				my $nCorrect = 0;
				my $nIncorrect = 0;
				
				my $mean_normalize = 0;
				foreach my $elem (@{$calibration_baskets{$i}{correct}})			
				{
					$nCorrect += $elem->{weight};
					$meanPP += $elem->{PP}*$elem->{weight};
					$mean_normalize += $elem->{weight};
				}
				foreach my $elem (@{$calibration_baskets{$i}{incorrect}})			
				{
					$nIncorrect += $elem->{weight};
					$meanPP += $elem->{PP}*$elem->{weight};
					$mean_normalize += $elem->{weight};				
				}			
				
				if(($nCorrect+$nIncorrect) != 0)
				{
					$percCorrect = $nCorrect/($nCorrect+$nIncorrect);
				}
				
				if($mean_normalize != 0)
				{
					$meanPP = $meanPP / $mean_normalize;
				}
				
				print CALIBRATION join("\t",
					$i,
					sprintf("%.2f", $meanPP),
					sprintf("%.2f", $percCorrect),
					sprintf("%.2f", $nCorrect),
					sprintf("%.2f", $nIncorrect),
				), "\n";
				
				$totalAlleles += ($nCorrect + $nIncorrect);
			}
			close(CALIBRATION);		
			
			die "Problem $totalAlleles vs $imputed_HLA_Calls{$locus}{called}" unless($totalAlleles == $imputed_HLA_Calls{$locus}{called});
					
			my $spatial_coverage_file = 'temp/spatialCoverage_' . $locus . '_' . $sample_IDs_abbr . '.txt';	
			open(SPATIALCOVERAGE, ">", $spatial_coverage_file) or die "Cannot open $spatial_coverage_file";			
			foreach my $exon (sort {$a <=> $b} keys %coverage_over_samples)
			{
				foreach my $exonPos (sort {$a <=> $b} keys %{$coverage_over_samples{$exon}})
				{
					my @individualValues = @{$coverage_over_samples_individualValues{$exon}{$exonPos}};
					@individualValues = sort {$a <=> $b} @individualValues;
					my $idx_10 = int($#individualValues * 0.1 + 0.5);
					my $idx_90 = int($#individualValues * 0.9 + 0.5);
					print SPATIALCOVERAGE join("\t", $exon, $exonPos, $exon . '-' . $exonPos, $coverage_over_samples{$exon}{$exonPos} / $coverage_over_samples_nSamples, $individualValues[$idx_10], $individualValues[$idx_90]), "\n";
				}
			}				
			close(SPATIALCOVERAGE)		
		}
				
		my $CR = sprintf("%.2f", $imputed_HLA_Calls{$locus}{sum} ? ($imputed_HLA_Calls{$locus}{called}/$imputed_HLA_Calls{$locus}{sum}) : 0);
		my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$locus} != 0) ? $problem_locus_detail{$locus}/$problem_locus_examined{$locus} : 0));
				
		print SUMMARY "\t", join("\t", $locus, $problem_locus_examined{$locus}, $CR, $accuracy), "\n";
	}

	my $comparions_OK = $comparisons - $compare_problems;
	print "\nComparisons: $comparisons -- OK: $comparions_OK\n";
		
	open(TMP_OUTPUT, '>', '../tmp/hla_validation/validation_summary.txt') or die;
	print "\nPER-LOCUS SUMMARY:\n";
	foreach my $key (sort keys %problem_locus_detail)
	{
		my $OK_locus = $problem_locus_examined{$key} - $problem_locus_detail{$key};
		my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$key} != 0) ? $problem_locus_detail{$key}/$problem_locus_examined{$key} : 0));
		my $CR = sprintf("%.2f", $imputed_HLA_Calls{$key}{sum} ? ($imputed_HLA_Calls{$key}{called}/$imputed_HLA_Calls{$key}{sum}) : 0);
		print "\t", $key, ": ", $OK_locus, " of ", $problem_locus_examined{$key}, ",\tAccuracy ", $accuracy, " ";
		print "\tCall rate: ", $CR,  "\n";
		
		my @fields = ($key, $problem_locus_examined{$key}, $CR, $accuracy);
		print TMP_OUTPUT join("\t", @fields), "\n";
	
	}
	close(TMP_OUTPUT);	
		
	print "\nCorrect vs incorrect coverages per locus:\n";
	foreach my $locus (sort keys %problem_locus_detail)
	{
		my @avg_minMax_ok = min_avg_max(@{$locus_avgCoverages{$locus}{ok}});
		my @low_minMax_ok = min_avg_max(@{$locus_lowCoverages{$locus}{ok}});
		my @min_minMax_ok = min_avg_max(@{$locus_minCoverages{$locus}{ok}});
			
		my @avg_minMax_problems = min_avg_max(@{$locus_avgCoverages{$locus}{problems}});
		my @low_minMax_problems = min_avg_max(@{$locus_lowCoverages{$locus}{problems}});
		my @min_minMax_problems = min_avg_max(@{$locus_minCoverages{$locus}{problems}});
		
		print "\t", $locus, "\n";

		print "\t\tAverage ", join(' / ', @avg_minMax_ok), " vs ", join(' / ', @avg_minMax_problems), " [problems]", "\n";
		print "\t\tLow ", join(' / ', @low_minMax_ok), " vs ", join(' / ', @low_minMax_problems), " [problems]", "\n";
		print "\t\tMin ", join(' / ', @min_minMax_ok), " vs ", join(' / ', @min_minMax_problems), " [problems]", "\n";
	}
	
	print "\n";
	close(SUMMARY);
}


if($actions =~ /w/)
{
	my $validation_round = 'R2';
	die "Please specify --trueHaplotypes for validation" unless($trueHaplotypes);
	
	my $debug = 0;
	
	# read reference dataset
	my %reference_data;
	open(REFERENCE, "<", $trueHaplotypes) or die "Cannot open $trueHaplotypes";
	my $headerLine = <REFERENCE>;
	chomp($headerLine);
	$headerLine =~ s/\n//g;
	$headerLine =~ s/\r//g;
	my @header_fields = split(/\t/, $headerLine);
	@header_fields = map {if($_ =~ /HLAD((QA)|(QB)|(RB))$/){$_ .= '1';} $_} @header_fields;	
	while(<REFERENCE>)
	{
		my $line = $_;
		chomp($line);
		
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		
		my @fields = split(/\t/, $line);
		my %line = (mesh @header_fields, @fields);
		
		my $primary_key = $line{'IndividualID'};
		$reference_data{$primary_key} = \%line;
	}
	close(REFERENCE);
	
	# perturbed positions
	my %reference_data_perturbedWhere;	
	my $trueHaplotypes_perturbed = $trueHaplotypes . '.perturbed';
	if(-e $trueHaplotypes_perturbed)
	{
		open(REFERENCE, "<", $trueHaplotypes_perturbed) or die "Cannot open $trueHaplotypes_perturbed";
		my $headerLine = <REFERENCE>;
		chomp($headerLine);
		$headerLine =~ s/\n//g;
		$headerLine =~ s/\r//g;
		my @header_fields = split(/\t/, $headerLine);
		@header_fields = map {if($_ =~ /HLAD((QA)|(QB)|(RB))$/){$_ .= '1';} $_} @header_fields;	
		while(<REFERENCE>)
		{
			my $line = $_;
			chomp($line);
			
			$line =~ s/\n//g;
			$line =~ s/\r//g;
			
			my @fields = split(/\t/, $line);
			my %line = (mesh @header_fields, @fields);
			
			my $primary_key = $line{'IndividualID'};
			$reference_data_perturbedWhere{$primary_key} = \%line;
		}
		close(REFERENCE);
	}
	else
	{
		warn "No perturbation information ($trueHaplotypes_perturbed) found.";
	}	
	
	
	my %imputed_haplotypes;
	my %sample_noI_toI;
	
	my %imputed_haplotypes_lines;
	
	$debug = 1;
	
	foreach my $sampleID (@sampleIDs)
	{
		my $sampleID_noI = $sampleID;
		$sampleID_noI =~ s/^I\d+_//g;
		
		my @bestGuess_files = glob('../tmp/hla/'.$sampleID.'/'.$validation_round.'_haplotypes_bestguess_*.txt');
		
		foreach my $bestGuess_file (@bestGuess_files)
		{
			die unless($bestGuess_file =~ /haplotypes_bestguess_(.+)\.txt/);
			my $locus = $1;
			
			# next unless($locus eq 'A'); # todo remove   
			
			unless(-e $bestGuess_file)
			{
				warn "Best-guess file $bestGuess_file not existing";
				next;
			}		
			  
			open(BESTGUESS, '<', $bestGuess_file) or die "Cannot open $bestGuess_file";
			my @bestguess_header_fields = qw/Position GraphLevel Type TypeNumber PositionInType C1 C2 GTConfidence PruningInterval PruningC1 PruningC1HaploConfidence PruningC2 PruningC2HaploConfidence NReads ReadAlleles ConsideredGTPairs ConsideredH1 ConsideredH2 AA1 AA1HaploConf AA2 AA2HaploConf/;
			while(<BESTGUESS>)
			{
				my $line = $_;
				chomp($line);
				my @line_fields = split(/\t/, $line, -1);
				die Dumper("Line fields problem", $#line_fields, $#bestguess_header_fields) unless($#line_fields == $#bestguess_header_fields);
				my %line_hash = (mesh @bestguess_header_fields, @line_fields);
				
				my $Q = $line_hash{'GTConfidence'};
				if($Q < $T)
				{
					$imputed_haplotypes{$locus}{$sampleID_noI}[0] .= '?;';	
					$imputed_haplotypes{$locus}{$sampleID_noI}[1] .= '?;';	
				}
				else
				{
					$imputed_haplotypes{$locus}{$sampleID_noI}[0] .= $line_hash{'C1'}.';';	
					$imputed_haplotypes{$locus}{$sampleID_noI}[1] .= $line_hash{'C2'}.';';
				}
				
				if($debug)
				{
					push(@{$imputed_haplotypes_lines{$locus}{$sampleID_noI}}, $line);
				}
			}	
			close(BESTGUESS);
			
			die if((exists $sample_noI_toI{$sampleID_noI}) and ($sample_noI_toI{$sampleID_noI} ne $sampleID));
			$sample_noI_toI{$sampleID_noI} = $sampleID;
		}
	}
		
	my @loci = sort keys %imputed_haplotypes;
	
	my $process_quality_measures = sub {};
	
	
	foreach my $locus (@loci)
	{
		my $arbitraty_indiv = (keys %reference_data)[0];
		next unless((defined $reference_data{$arbitraty_indiv}{'HLA'.$locus}));
		
		my $locus_agree = 0;
		my $locus_disagree = 0;
		my $locus_missing = 0;

		my $locus_gt_agree = 0;
		my $locus_gt_disagree = 0;
		my $locus_gt_missing = 0;
					
		my $locus_pt_agree = 0;
		my $locus_pt_disagree = 0;
		my $locus_pt_missing = 0;

		my $locus_pt_gt_agree = 0;
		my $locus_pt_gt_disagree = 0;
		my $locus_pt_gt_missing = 0;
		
		my @indivIDs = keys %{$imputed_haplotypes{$locus}};
		INDIV: foreach my $indivID (@indivIDs)
		{	
			$debug = 0;
			
			# my @imputed_hla_values = map { $imputed_HLA{$locus}{$indivID}{$_} } keys %{$imputed_HLA{$locus}{$indivID}};
			my @imputed_haplotypes = @{ $imputed_haplotypes{$locus}{$indivID} };
								
			next INDIV unless($#imputed_haplotypes == 1);
			next INDIV unless(exists $reference_data{$indivID});
			next INDIV unless(exists $reference_data{$indivID}{'HLA'.$locus});
			
			$reference_data{$indivID}{'HLA'.$locus} or die;
			
			my @reference_haplotypes = split(/\//, $reference_data{$indivID}{'HLA'.$locus});
			if(not $reference_data_perturbedWhere{$indivID}{'HLA'.$locus})
			{
				warn "No reference haplotype for $locus";
				$reference_data_perturbedWhere{$indivID}{'HLA'.$locus} = '';
			}
			my @reference_haplotypes_perturbed = split(/\//, $reference_data_perturbedWhere{$indivID}{'HLA'.$locus});
			
			die Dumper($reference_data{$indivID}, \@reference_haplotypes) unless($#reference_haplotypes == 1);
				
			die unless($#imputed_haplotypes == 1);
			die unless($#reference_haplotypes == 1);
			
			my @reference_haplotypes_split = (map {[split(/;/, $_)]} @reference_haplotypes);
			my @reference_haplotypes_perturbed_split = (map {split(/;/, $_)} @reference_haplotypes_perturbed);
			
			my %perturbedWhere = map {$_ => 1} @reference_haplotypes_perturbed_split;
			
			# die Dumper(\%perturbedWhere);

			my @imputed_haplotypes_split = (map {[split(/;/, $_)]} @imputed_haplotypes);
			
			die unless($#reference_haplotypes_split == 1);
			die unless($#imputed_haplotypes_split == 1);
			
			die Dumper($#reference_haplotypes_split, $#imputed_haplotypes_split, [$#{$reference_haplotypes_split[0]}, $#{$reference_haplotypes_split[1]}], [$#{$imputed_haplotypes_split[0]}, $#{$imputed_haplotypes_split[1]}]) unless($#{$reference_haplotypes_split[0]} == $#{$reference_haplotypes_split[1]});
			
			die Dumper($#reference_haplotypes_split, $#imputed_haplotypes_split, [$#{$reference_haplotypes_split[0]}, $#{$reference_haplotypes_split[1]}], [$#{$imputed_haplotypes_split[0]}, $#{$imputed_haplotypes_split[1]}]) unless($#{$imputed_haplotypes_split[0]} == $#{$imputed_haplotypes_split[1]});

			die Dumper($#reference_haplotypes_split, $#imputed_haplotypes_split, [$#{$reference_haplotypes_split[0]}, $#{$reference_haplotypes_split[1]}], [$#{$imputed_haplotypes_split[0]}, $#{$imputed_haplotypes_split[1]}]) unless($#{$reference_haplotypes_split[0]} == $#{$imputed_haplotypes_split[1]});
			
			my @alleles_agree = (0, 0);
			my @alleles_missing = (0, 0);;
			my @alleles_disagree = (0, 0);;

			my @alleles_gt_agree = (0, 0);
			my @alleles_gt_missing = (0, 0);;
			my @alleles_gt_disagree = (0, 0);;
			
			my @alleles_pt_agree = (0, 0);
			my @alleles_pt_missing = (0, 0);;
			my @alleles_pt_disagree = (0, 0);;

			my @alleles_pt_gt_agree = (0, 0);
			my @alleles_pt_gt_missing = (0, 0);;
			my @alleles_pt_gt_disagree = (0, 0);;
			
			
			for(my $invertImputations = 0; $invertImputations <= 1; $invertImputations++)
			{
				my @imputed_haplotypes_split_forAnalysis = ($invertImputations) ? reverse(@imputed_haplotypes_split) : @imputed_haplotypes_split;
				
				for(my $i = 0; $i <= $#{$reference_haplotypes_split[0]}; $i++)
				{
					for(my $j = 0; $j <= 1; $j++)
					{
						die unless((defined $reference_haplotypes_split[$j][$i]) and (defined $imputed_haplotypes_split_forAnalysis[$j][$i]));
						if($imputed_haplotypes_split_forAnalysis[$j][$i] eq "?")
						{
							$alleles_missing[$invertImputations]++;
							if($perturbedWhere{$i})
							{
								$alleles_pt_missing[$invertImputations]++;
							}
						}
						else
						{
							if($reference_haplotypes_split[$j][$i] eq $imputed_haplotypes_split_forAnalysis[$j][$i])
							{
								$alleles_agree[$invertImputations]++;
								if($perturbedWhere{$i})
								{
									$alleles_pt_agree[$invertImputations]++;
								}								
							}
							else
							{
								$alleles_disagree[$invertImputations]++;
								if($perturbedWhere{$i})
								{
									$alleles_pt_disagree[$invertImputations]++;
								}									
							}
						}
					}
					
					my ($thisPosition_gt_agree, $thisPosition_gt_disagree, $thisPosition_gt_missing) = compatibleStringAlleles([$reference_haplotypes_split[0][$i], $reference_haplotypes_split[1][$i]], [$imputed_haplotypes_split_forAnalysis[0][$i], $imputed_haplotypes_split_forAnalysis[1][$i]]);
					
					$alleles_gt_agree[$invertImputations] += $thisPosition_gt_agree;
					$alleles_gt_disagree[$invertImputations] += $thisPosition_gt_disagree;
					$alleles_gt_missing[$invertImputations] += $thisPosition_gt_missing;
					
					if($perturbedWhere{$i})
					{					
						$alleles_pt_gt_agree[$invertImputations] += $thisPosition_gt_agree;
						$alleles_pt_gt_disagree[$invertImputations] += $thisPosition_gt_disagree;
						$alleles_pt_gt_missing[$invertImputations] += $thisPosition_gt_missing;
					
						if(($invertImputations == 0) and ($thisPosition_gt_agree != 2))
						{
							print "Position $i -- agreement: $thisPosition_gt_agree\n";
							print "\tTrue genotypes: ", join('/', $reference_haplotypes_split[0][$i], $reference_haplotypes_split[1][$i]), "\n";
							print "\tImputed genotypes: ", join('/', $imputed_haplotypes_split_forAnalysis[0][$i], $imputed_haplotypes_split_forAnalysis[1][$i]), "\n";
							print "\t\tLine: ",	$imputed_haplotypes_lines{$locus}{$indivID}[$i], "\n";
							print "\n";
						}					
					}
					
					if(($invertImputations == 0) and ($thisPosition_gt_agree != 2))
					{
						# print "Position $i -- agreement: $thisPosition_gt_agree\n";
						# print "\tTrue genotypes: ", join('/', $reference_haplotypes_split[0][$i], $reference_haplotypes_split[1][$i]), "\n";
						# print "\tImputed genotypes: ", join('/', $imputed_haplotypes_split_forAnalysis[0][$i], $imputed_haplotypes_split_forAnalysis[1][$i]), "\n";
						# print "\t\tLine: ",	$imputed_haplotypes_lines{$locus}{$indivID}[$i], "\n";
						# print "\n";
					}
				}
			}
			
			die unless($alleles_gt_agree[0] == $alleles_gt_agree[1]);
			die unless($alleles_gt_disagree[0] == $alleles_gt_disagree[1]);
			die unless($alleles_gt_missing[0] == $alleles_gt_missing[1]);
			
			# print "Individual $indivID Locus $locus\n";
			my @invertImputations_notOK;
			for(my $invertImputations = 0; $invertImputations <= 1; $invertImputations++)
			{	
				my $alleles_sum = $alleles_agree[$invertImputations] + $alleles_disagree[$invertImputations] + $alleles_missing[$invertImputations];
				die unless($alleles_sum > 0);
				
				# print "\t", "Inversion ", $invertImputations, "\n";
				# print "\t\tOK:       ", $alleles_agree[$invertImputations], " ", sprintf("%.2f", $alleles_agree[$invertImputations]/$alleles_sum * 100), "%\n";
				# print "\t\tNOT OK:   ", $alleles_disagree[$invertImputations], " ", sprintf("%.2f", $alleles_disagree[$invertImputations]/$alleles_sum * 100), "%\n";
				# print "\t\tMISSING:  ", $alleles_missing[$invertImputations], " ", sprintf("%.2f", $alleles_missing[$invertImputations]/$alleles_sum * 100), "%\n";
				# print "\n";
				
				$invertImputations_notOK[$invertImputations] = $alleles_disagree[$invertImputations];
				
				my $alleles_pt_sum = $alleles_pt_agree[$invertImputations] + $alleles_pt_disagree[$invertImputations] + $alleles_pt_missing[$invertImputations];
				
				# print "\t", "Inversion ", $invertImputations, " (perturbed alleles only)\n";
				# if($alleles_pt_sum > 0)		
				# {
					# print "\t\tOK:       ", $alleles_pt_agree[$invertImputations], " ", sprintf("%.2f", $alleles_pt_agree[$invertImputations]/$alleles_pt_sum * 100), "%\n";
					# print "\t\tNOT OK:   ", $alleles_pt_disagree[$invertImputations], " ", sprintf("%.2f", $alleles_pt_disagree[$invertImputations]/$alleles_pt_sum * 100), "%\n";
					# print "\t\tMISSING:  ", $alleles_pt_missing[$invertImputations], " ", sprintf("%.2f", $alleles_pt_missing[$invertImputations]/$alleles_pt_sum * 100), "%\n";
				# }
				# print "\n";
				
			}
			
			my $invertImputations_optimal = ($invertImputations_notOK[1] < $invertImputations_notOK[0]) ? 1 : 0;
			
			my $alleles_sum = $alleles_agree[$invertImputations_optimal] + $alleles_disagree[$invertImputations_optimal] + $alleles_missing[$invertImputations_optimal];
			die unless($alleles_sum > 0);
						
			
			# print "\t", "Haplotypes - inverted ", $invertImputations_optimal, "\n";
			# print "\t\tOK:       ", $alleles_agree[$invertImputations_optimal], " ", sprintf("%.2f", $alleles_agree[$invertImputations_optimal]/$alleles_sum * 100), "%\n";
			# print "\t\tNOT OK:   ", $alleles_disagree[$invertImputations_optimal], " ", sprintf("%.2f", $alleles_disagree[$invertImputations_optimal]/$alleles_sum * 100), "%\n";
			# print "\t\tMISSING:  ", $alleles_missing[$invertImputations_optimal], " ", sprintf("%.2f", $alleles_missing[$invertImputations_optimal]/$alleles_sum * 100), "%\n";
			# print "\n";
				
			
			my $alleles_gt_sum = $alleles_gt_agree[0] + $alleles_gt_disagree[0] + $alleles_gt_missing[0];
			# print "\t", "Genotypes ", "\n";
			# print "\t\tOK:       ", $alleles_gt_agree[0], " ", sprintf("%.2f", $alleles_gt_agree[0]/$alleles_gt_sum * 100), "%\n";
			# print "\t\tNOT OK:   ", $alleles_gt_disagree[0], " ", sprintf("%.2f", $alleles_gt_disagree[0]/$alleles_gt_sum * 100), "%\n";
			# print "\t\tMISSING:  ", $alleles_gt_missing[0], " ", sprintf("%.2f", $alleles_gt_missing[0]/$alleles_gt_sum * 100), "%\n";
			# print "\n";			
			
			my $alleles_pt_gt_sum = $alleles_pt_gt_agree[0] + $alleles_pt_gt_disagree[0] + $alleles_pt_gt_missing[0];
			# print "\t", "Genotypes at perturbed positions", "\n";
			# if($alleles_pt_gt_sum > 0)
			# {
				# print "\t\tOK:       ", $alleles_pt_gt_agree[0], " ", sprintf("%.2f", $alleles_pt_gt_agree[0]/$alleles_pt_gt_sum * 100), "%\n";
				# print "\t\tNOT OK:   ", $alleles_pt_gt_disagree[0], " ", sprintf("%.2f", $alleles_pt_gt_disagree[0]/$alleles_pt_gt_sum * 100), "%\n";
				# print "\t\tMISSING:  ", $alleles_pt_gt_missing[0], " ", sprintf("%.2f", $alleles_pt_gt_missing[0]/$alleles_pt_gt_sum * 100), "%\n";
			# }
			# print "\n";					
								
			$locus_agree += $alleles_agree[$invertImputations_optimal];
			$locus_disagree += $alleles_disagree[$invertImputations_optimal];
			$locus_missing += $alleles_missing[$invertImputations_optimal];
			
			$locus_gt_agree += $alleles_gt_agree[0];
			$locus_gt_disagree += $alleles_gt_disagree[0];
			$locus_gt_missing += $alleles_gt_missing[0];	

			$locus_pt_agree += $alleles_pt_agree[$invertImputations_optimal];
			$locus_pt_disagree += $alleles_pt_disagree[$invertImputations_optimal];
			$locus_pt_missing += $alleles_pt_missing[$invertImputations_optimal];
			
			$locus_pt_gt_agree += $alleles_pt_gt_agree[0];
			$locus_pt_gt_disagree += $alleles_pt_gt_disagree[0];
			$locus_pt_gt_missing += $alleles_pt_gt_missing[0];	
		}
		
		my $locus_sum = $locus_agree + $locus_disagree + $locus_missing;
		my $locus_gt_sum = $locus_gt_agree + $locus_gt_disagree + $locus_gt_missing;
		
		die if($locus_sum == 0);
		print "LOCUS SUMMARY $locus\n";
		print "\tHaplotypes:\n";
		print "\t\tOK:       ", $locus_agree, " ", sprintf("%.2f", $locus_agree/$locus_sum * 100), "%\n";
		print "\t\tNOT OK:   ", $locus_disagree, " ", sprintf("%.2f", $locus_disagree/$locus_sum * 100), "%\n";
		print "\t\tMISSING:  ", $locus_missing, " ", sprintf("%.2f", $locus_missing/$locus_sum * 100), "%\n\n";
		print "\tGenotypes:\n";
		print "\t\tOK:       ", $locus_gt_agree, " ", sprintf("%.2f", $locus_gt_agree/$locus_gt_sum * 100), "%\n";
		print "\t\tNOT OK:   ", $locus_gt_disagree, " ", sprintf("%.2f", $locus_gt_disagree/$locus_gt_sum * 100), "%\n";
		print "\t\tMISSING:  ", $locus_gt_missing, " ", sprintf("%.2f", $locus_gt_missing/$locus_gt_sum * 100), "%\n";
		print "\n";
		
		my $locus_pt_sum = $locus_pt_agree + $locus_pt_disagree + $locus_pt_missing;
		my $locus_pt_gt_sum = $locus_pt_gt_agree + $locus_pt_gt_disagree + $locus_pt_gt_missing;
		
		print "LOCUS SUMMARY $locus (perturbed only)\n";
		if($locus_pt_sum > 0)
		{
			print "\tHaplotypes:\n";
			print "\t\tOK:       ", $locus_pt_agree, " ", sprintf("%.2f", $locus_pt_agree/$locus_pt_sum * 100), "%\n";
			print "\t\tNOT OK:   ", $locus_pt_disagree, " ", sprintf("%.2f", $locus_pt_disagree/$locus_pt_sum * 100), "%\n";
			print "\t\tMISSING:  ", $locus_pt_missing, " ", sprintf("%.2f", $locus_pt_missing/$locus_pt_sum * 100), "%\n\n";
			print "\tGenotypes:\n";
			print "\t\tOK:       ", $locus_pt_gt_agree, " ", sprintf("%.2f", $locus_pt_gt_agree/$locus_pt_gt_sum * 100), "%\n";
			print "\t\tNOT OK:   ", $locus_pt_gt_disagree, " ", sprintf("%.2f", $locus_pt_gt_disagree/$locus_pt_gt_sum * 100), "%\n";
			print "\t\tMISSING:  ", $locus_pt_gt_missing, " ", sprintf("%.2f", $locus_pt_gt_missing/$locus_pt_gt_sum * 100), "%\n";
		}	
	
		# exit; # todo remove
	}

	# my $comparions_OK = $comparisons - $compare_problems;
	# print "\nComparisons: $comparisons -- OK: $comparions_OK\n";
		
	open(TMP_OUTPUT, '>', '../tmp/hla_validation/validation_haplotypes_summary.txt') or die;
	# print "\nPER-LOCUS SUMMARY:\n";
	# foreach my $key (sort keys %problem_locus_detail)
	# {
		# my $OK_locus = $problem_locus_examined{$key} - $problem_locus_detail{$key};
		# my $accuracy = sprintf("%.2f", 1 - (($problem_locus_examined{$key} != 0) ? $problem_locus_detail{$key}/$problem_locus_examined{$key} : 0));
		# my $CR = sprintf("%.2f", $imputed_HLA_Calls{$key}{sum} ? ($imputed_HLA_Calls{$key}{called}/$imputed_HLA_Calls{$key}{sum}) : 0);
		# print "\t", $key, ": ", $OK_locus, " of ", $problem_locus_examined{$key}, ",\tAccuracy ", $accuracy, " ";
		# print "\tCall rate: ", $CR,  "\n";
		
		# my @fields = ($key, $problem_locus_examined{$key}, $CR, $accuracy);
		# print TMP_OUTPUT join("\t", @fields), "\n";
	# }
	close(TMP_OUTPUT);	

	print "\n";
}


sub compatibleStringAlleles
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my ($OK1, $NOTOK1, $MISSING1) = compatibleStringAlleles_noFlip($alleles_validation, $alleles_inference);
	
	my $alleles_validation_flipped = [reverse(@$alleles_validation)];
	my ($OK2, $NOTOK2, $MISSING2) = compatibleStringAlleles_noFlip($alleles_validation_flipped, $alleles_inference);
	
	# die unless($MISSING1 == $MISSING2);
	
	if(($NOTOK1 < $NOTOK2))
	{
		return ($OK1, $NOTOK1, $MISSING1);
	}
	elsif($NOTOK1 == $NOTOK2)
	{
		return (($OK1 > $OK2) ? ($OK1, $NOTOK1, $MISSING1) : ($OK2, $NOTOK2, $MISSING2));		
	}
	else
	{
		return ($OK2, $NOTOK2, $MISSING2);
	}
	
	
}


sub compatibleStringAlleles_noFlip
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my $alleles_OK = 0;
	my $alleles_NOTOK = 0;
	my $alleles_MISSING = 0;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{	
		my ($OK, $NOTOK, $MISSING) = (compatibleStringAlleles_individual(undef, $alleles_validation->[$aI], $alleles_inference->[$aI]));
		$alleles_OK += $OK;
		$alleles_NOTOK += $NOTOK;
		$alleles_MISSING += $MISSING;
	}
		
	return ($alleles_OK, $alleles_NOTOK, $alleles_MISSING);
}


sub compatibleStringAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;	
	my $allele_inference = shift;

	if(($allele_validation =~ /\?/) || ($allele_inference =~ /\?/))
	{
		return (0, 0, 1);
	}
	else
	{
		if($allele_validation eq $allele_inference)
		{
			return (1, 0, 0);
		}
		else
		{
			return (0, 1, 0);
		}
	}
}



sub compatibleAlleles
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my $r1 = compatibleAlleles_noFlip($alleles_validation, $alleles_inference);
	
	my $alleles_validation_flipped = [reverse(@$alleles_validation)];
	my $r2 = compatibleAlleles_noFlip($alleles_validation_flipped, $alleles_inference);
	
	return (($r1 > $r2) ? $r1 : $r2);
}

sub compatibleAlleles_noFlip
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	
	my $alleles_compatible = 0;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{	
		$alleles_compatible += (compatibleAlleles_individual(undef, $alleles_validation->[$aI], $alleles_inference->[$aI]));
	}
		
	return $alleles_compatible;
}

sub compatibleAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;
		my $allele_validation_original = $allele_validation;
		
	my $allele_inference = shift;
	
	die unless(length($allele_validation) >= 4); die unless($allele_validation =~ /^\d\d/);
	
	my $components_validation;
	if($allele_validation !~ /\:/)
	{
		$allele_validation = substr($allele_validation, 0, 2).':'.substr($allele_validation, 2);
		$components_validation = 2;
	}
	else
	{
		my @_components = split(/\:/, $allele_validation);
		$components_validation = scalar(@_components);
	}
	die unless(defined $components_validation);
	
	
	my $true_allele = $allele_validation;
	$true_allele =~ s/g//;
		
	my $inferred_alleles = $allele_inference;
	my @inferred_alleles = split(/;/, $inferred_alleles);
	
	@inferred_alleles = map {
		# die "Can't parse allele $_" unless($_ =~ /^\s*\w+\*(\d\d\d?)\:(\d\d\d?N?).*$/);
		# my $allele = $1.':'.$2;
		# $allele
		
		die "Can't parse allele $_" unless($_ =~ /^\s*(\w+)\*([\d\:N]+Q?)L?S?$/);
		
		my $allele = $2;
		
		my @components_allele = split(/:/, $allele);
		my @components_allele_rightLength;
		
		for(my $i = 0; $i < $components_validation; $i++)
		{
			if($i <= $#components_allele)
			{
				push(@components_allele_rightLength, $components_allele[$i]);
			}
			else
			{
				push(@components_allele_rightLength, '00');
			}
		}
		$allele = join(':', @components_allele_rightLength);
		
	} @inferred_alleles;
	
	my %inferred = map {$_ => 1} @inferred_alleles;
	
	if($inferred{$true_allele})
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 1), "\n" if ($locus eq 'DQB1');	
		
		return 1;  
	}		
	else
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 0), "\n"  if ($locus eq 'DQB1');	
	
		return 0;
	}
}

sub intersection
{
	my $aref1 = shift;
	my $aref2 = shift;

	my %h1 = map {$_ => 1} @$aref1;
	my %h2 = map {$_ => 1} @$aref2;
	
	my %c = map {$_ => 1} ((keys %h1), (keys %h2));
	
	my @ret = grep {$h1{$_} and $h2{$_}} keys %c;
	
	return @ret;
}

sub print_which_exons
{
	my $locus = shift;
	if((length($locus) == 1) or ($locus =~ /HLA\w/))
	{
		return (2, 3);
	}
	else
	{
		return (2);
	}
}	


sub load_pileup
{
	my $r_href = shift;
	my $file = shift;
	my $indivID = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line, -1);
		die unless(($#f == 3) or ($#f == 2));
		if($#f == 3)
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = $f[3];		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		}
		else
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = '';		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		}
	}	
	close(F);
}

sub load_coverages_from_pileup
{
	my $file = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	my %forReturn;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die "Cannot parse pileup coverage line" unless($line =~ /^(\d+)\s(\d+)\s(\d+)/);
		my $exon = $1;
		my $exonPos = $2;
		my $coverage = $3;				
		
		die "Pos twice in $file ? $exon $exonPos" if(defined $forReturn{$exon}{$exonPos});
		$forReturn{$exon}{$exonPos} = $coverage;
	}	
	close(F);
	
	return \%forReturn;
}


sub twoValidationAlleles_2_proper_names
{
	my $alleles_validation = shift;
	my $locus = shift;
	my $exon_sequences = shift;
	
	die unless($#{$alleles_validation} == 1);
	
	my @forReturn;
	
	$locus =~ s/HLA//;
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my $validation_allele = $alleles_validation->[$aI];
		
		if($validation_allele =~ /\:/)
		{
			$validation_allele = $locus . '*' . $validation_allele;
			# $validation_allele =~ s/g//;
			# my @components = split(/\:/, $validation_alleles),
			
			# die unless(length($validation_allele) >= 4);
			
			# $validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
			# $validation_allele = $locus . '*' . $validation_allele;
		
		}
		else
		{
			$validation_allele =~ s/g//;
			die unless(length($validation_allele) >= 4);
			
			$validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
			$validation_allele = $locus . '*' . $validation_allele;
		}
		my @validation_alleles = ($validation_allele);
		
		my $extensions = 0;
		my $notFoundFirstRound = 0;
		while( not exists $exon_sequences->{$validation_allele} )
		{
			$extensions++;
		
			$validation_allele .= ':01';
			if($extensions > 3)
			{
				$notFoundFirstRound = 1;
				last;
			}
		}
		
		if($notFoundFirstRound)
		{
			my $extensions = 0;
			while(scalar(grep {exists $exon_sequences->{$_}} @validation_alleles) == 0)
			{
				$extensions++;
				my $viMax = $#validation_alleles;
				for(my $vI = 0; $vI <= $viMax; $vI++)
				{
					my $a1 = $validation_alleles[$vI];
					$validation_alleles[$vI] .= ':01';
					push(@validation_alleles, $a1.':02');
				}
				
				$validation_allele .= ':01';
				if($extensions > 3)
				{
					print Dumper([grep {$_ =~ /DQB1/} keys %$exon_sequences]);
					die Dumper("Can't identify (II) exon alleles for $alleles_validation->[$aI]", \@validation_alleles, $locus, $alleles_validation, $extensions, [(keys %$exon_sequences)[0 .. 10]]);
				}
			}		
			
			my @foundAlleles = grep {exists $exon_sequences->{$_}} @validation_alleles;
			die unless(scalar(@foundAlleles) > 0);
			push(@forReturn, $foundAlleles[0]);		
			
		}
		else
		{
			push(@forReturn, $validation_allele);		
		}
	}
			
	return @forReturn;
}

sub list_comparison
{
	my $list1_aref = shift;
	my $list2_aref = shift;
	
	my %l1 = map {$_ => 1} @$list1_aref;
	my %l2 = map {$_ => 1} @$list2_aref;
	
	unless(scalar(keys %l1) == scalar(@$list1_aref))
	{
		die "List 1 non-unique elements";
	}
	
	unless(scalar(keys %l2) == scalar(@$list2_aref))
	{
		die "List 2 non-unique elements";
	}	
	
	my %combined = map {$_ => 1} (@$list1_aref, @$list2_aref);
	
	my @l1_exclusive;
	my @l2_exclusive;
	my $n_shared = 0;
	foreach my $e (keys %combined)
	{
		if($l1{$e} and $l2{$e})
		{
			$n_shared++;
		}
		elsif($l1{$e})
		{
			push(@l1_exclusive, $e);
		}
		elsif($l2{$e})
		{
			push(@l2_exclusive, $e);
		}
		else
		{
			die;
		}
				
	}
	
	die unless(($n_shared + scalar(@l1_exclusive) + scalar(@l2_exclusive)) == scalar(keys %combined));
	
	return ($n_shared, \@l1_exclusive, \@l2_exclusive);
}

sub twoClusterAlleles
{
	my $alleles_inference = shift;
	die Dumper("Don't have two inferred alleles", $alleles_inference) unless($#{$alleles_inference} == 1);
	
	my @forReturn;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my $inferred_alleles = $alleles_inference->[$aI];
		my @inferred_alleles = split(/;/, $inferred_alleles);
	
		$inferred_alleles[0] =~ s/^\s+//;
		
		push(@forReturn, $inferred_alleles[0]);
	}
	
	return @forReturn;
}


sub read_exon_sequences
{
	my $file = shift;
	my %r;
	open(F, '<', $file) or die "Cannot open $file";
	my $firstLine = <F>;
	chomp($firstLine);
	my @header_fields = split(/ /, $firstLine);
	die unless($header_fields[0] eq 'IndividualID');
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_fields = split(/ /, $line);
		my $allele = shift(@line_fields);
		$r{$allele} = join('', @line_fields);
	}
	close(F);
	
	return \%r;
}


sub find_exon_file
{
	my $locus = shift;
	my $exon = shift;
	$locus =~ s/HLA//;
	
	opendir(my $dh, $exon_folder) or die;
	my @exon_folder_files = readdir($dh);
	closedir($dh);
	
	@exon_folder_files = grep {$_ =~ /^${locus}_\d+_exon_${exon}\.txt$/} @exon_folder_files;
	if($#exon_folder_files != 0)
	{
		die Dumper("Can't find exon file for $locus // $exon", @exon_folder_files);
	}
	else
	{
		return $exon_folder.'/'.$exon_folder_files[0];
	}	
}


sub min_avg_max
{
	my @v = @_;
	@v = sort {$a <=> $b} @v;
	if($#v >= 1)
	{
		die unless($v[0] <= $v[1]);
	}
	
	if(scalar(@v) == 0)
	{
		return ('', '', '');
	}
	if(scalar(@v) == 1)
	{
		return($v[0], $v[0], '');
	}
	
	my $min = $v[0];
	my $max = $v[$#v];
	  
	my $sum = 0;
	foreach my $vE (@v)
	{
		$sum += $vE;
	}
	
	my $avg = '';
	if(scalar(@v) > 0)
	{
		$avg = $sum / scalar(@v);
	}
	
	return ($min, $avg, $max);
}
