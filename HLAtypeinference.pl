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

my $kMer_size = 55;  

# input parameters

my $graph = 'hla';   
my $sampleIDs = '';
my $BAMs = '';
my $actions;
my $trueHLA;
my $validation_round = 'R1';
my $T = 0;
my $all_2_dig = 0;
my $only_4_dig = 1;

GetOptions ('graph:s' => \$graph,
 'sampleIDs:s' => \$sampleIDs, 
 'BAMs:s' => \$BAMs, 
 'actions:s' => \$actions, 
 'trueHLA:s' => \$trueHLA,
 'validation_round:s' => \$validation_round,
 'T:s' => \$T,
);         


my $genome_graph_file = qq(../tmp2/GS_nextGen/hla/derived/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph_k25.ctx);
unless(-e $genome_graph_file)
{
	die "Please set variable \$genome_graph_file to an existing file - the current value $genome_graph_file is not accessible.";
}
my $expected_kMer_file = qq(../tmp2/GS_nextGen/${graph}/requiredkMers_graph.txt.kmers_25);
my $exon_folder = qq(../tmp2/GS_nextGen/${graph}/);
unless(-e $expected_kMer_file)
{
	die "Please provide a kMerified graph -- exepcted file $expected_kMer_file not there!";
}

my $normal_bin = qq(../bin/MHC-PRG);
my $cluster3_bin = qq(../bin_cluster3/MHC-PRG);
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

if($sampleIDs =~ /^all/)
{
	my @dirs = glob('../tmp/hla/*');
	@sampleIDs = map {die "Can't parse $_" unless($_ =~ /tmp\/hla\/(.+)/); $1} @dirs;
	
	if($sampleIDs =~ /^all_I(\d+)/i)
	{
		my $iteration = $1;
		@sampleIDs = grep {$_ =~ /^I${iteration}_/i} @sampleIDs;
	}
}

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
		
		my $command = qq($use_bin domode filterReads --input_BAM $BAM --positiveFilter $expected_kMer_file --output_FASTQ $output_file);
		
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
	
	my $command = qq($use_bin domode filterReads --input_FASTQ $fastQ_files --negativeFilter $genome_graph_file --output_FASTQ $output_files);
	
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
	}
		
	for(my $sI = 0; $sI <= $#aligned_files; $sI++)
	{
		my $sampleID = $sampleIDs[$sI];
		my $aligned_file = $aligned_files[$sI];
		my $stdout_file = $stdout_files[$sI];
		
		my $command = qq($use_bin domode HLATypeInference --input_alignedReads $aligned_file --graphDir ../tmp2/GS_nextGen/${graph} --sampleID $sampleID > $stdout_file);
	
		print "Now executing command:\n$command\n\n";
		
		system($command);		
	}
}  

if($actions =~ /v/)
{
	die "Please specify --trueHLA for validation" unless($trueHLA);
		
	# read reference dataset
	my %reference_data;
	open(REFERENCE, "<", $trueHLA) or die "Cannot open $trueHLA";
	my $headerLine = <REFERENCE>;
	chomp($headerLine);
	$headerLine =~ s/\n//g;
	$headerLine =~ s/\r//g;
	my @header_fields = split(/\t/, $headerLine);
	@header_fields = map {if($_ =~ /HLAD((QA)|(QB)|(RB))/){$_ .= '1';} $_} @header_fields;	
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
	
	# die Dumper(\%reference_data);
	
	my %imputed_HLA;
	my %sample_noI_toI;
	
	foreach my $sampleID (@sampleIDs)
	{
		my $sampleID_noI = $sampleID;
		$sampleID_noI =~ s/^I\d+_//g;
		
		my $bestGuess_file = '../tmp/hla/'.$sampleID.'/'.$validation_round.'_bestguess.txt';
		unless(-e $bestGuess_file)
		{
			warn "Best-guess file $bestGuess_file not existing";
			next;
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
			if($Q < $T)
			{
				$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = '0000';			
			}
			else
			{
				$imputed_HLA{$line_hash{'Locus'}}{$sampleID_noI}{$line_hash{'Chromosome'}} = $line_hash{'Allele'};			
			}
		}	
		close(BESTGUESS);
		
		die if(exists $sample_noI_toI{$sampleID_noI});
		$sample_noI_toI{$sampleID_noI} = $sampleID;
	}
	
	my $debug = 0;
	my $comparisons = 0;
	my $compare_problems = 0;
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
	
	foreach my $locus (@loci)
	{
		my $arbitraty_indiv = (keys %reference_data)[0];
		next unless((defined $reference_data{$arbitraty_indiv}{'HLA'.$locus}));
		
		$problem_locus_examined{$locus} = 0;
		$problem_locus_detail{$locus} = 0;
		my @indivIDs = keys %{$imputed_HLA{$locus}};
		INDIV: foreach my $indivID (@indivIDs)
		{	
			$debug = 0;
			
			my @imputed_hla_values = map { $imputed_HLA{$locus}{$indivID}{$_} } keys %{$imputed_HLA{$locus}{$indivID}};
			my @reference_hla_values;
			
			next INDIV unless($#imputed_hla_values == 1);
			next INDIV unless(exists $reference_data{$indivID});
			
			$reference_data{$indivID}{'HLA'.$locus} or die;
			
			if($reference_data{$indivID}{'HLA'.$locus})
			{
				@reference_hla_values = split(/\//, $reference_data{$indivID}{'HLA'.$locus});
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
				@imputed_hla_values = map {join(';', map {&simpleHLA::HLA_2digit($_)} split(/;/, $_))} @imputed_hla_values;
			}
			
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
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
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
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[1]))
								{
									$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{incorrect}++;									
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
								{
									$reference_predictions{$locus}{$reference_hla_values[1]}{correct}++;
									$imputations_predictions{$locus}{$imputed_hla_values[0]}{correct}++;					
									$reference_predictions{$locus}{$reference_hla_values[0]}{incorrect}++;
									$imputations_predictions{$locus}{$imputed_hla_values[1]}{incorrect}++;						
								}
								elsif(&compatibleAlleles_individual($locus, $reference_hla_values[0], $imputed_hla_values[1]))
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
							}
							elsif(&compatibleAlleles_individual($locus, $reference_hla_values[1], $imputed_hla_values[0]))
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
					}
					else
					{
						die;
					}
				}		
			}
			
			my $thisIndiv_problems = $problem_locus_detail{$locus} - $problem_locus_detail_before;
			
			my $thisIndiv_comparions = $comparisons - $comparisons_before;
			my $thisIndiv_OK = $thisIndiv_comparions - $thisIndiv_problems;
			
			if(1 or ($thisIndiv_problems > 0))
			{
				my $indivID_withI = $sample_noI_toI{$indivID};
				die unless(defined $indivID_withI);
				
				my $pileup_file = qq(../tmp/hla/$indivID_withI/${validation_round}_pileup_${locus}.txt);
				
				unless(-e $pileup_file)
				{
					die "Can't find pileup file $pileup_file";
				}
				
				load_pileup($pileup_href, $pileup_file, $indivID_withI);

			
				my $output_fn = '../tmp/hla_validation/pileup_'.$validation_round.'_'.$indivID_withI.'_'.$locus.'.txt';
				open(my $output_fh, '>', $output_fn) or die "Cannot open $output_fn";
				print $output_fh join("\t", $indivID_withI, $locus, $thisIndiv_OK), "\n";
				
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
							# next unless(defined $pileUpString);
							
							# die Dumper("Pileup too short", $length, scalar(@{$pileup_href->{$indivID_withI}{'HLA'.$locus}[$exon-2]})) unless(defined $pileUpString);
							unless(($chars_for_print[5] eq '!') or ($chars_for_print[6] eq '!'))
							{
								$pileUpString =~ s/\[.+?\]//g;
							}
							push(@chars_for_print, $pileUpString);
						}
						print {$output_fh} join("\t", @chars_for_print), "\n";
					}	
				}
				
				close($output_fh);					
			}
		}
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

	print "\n";
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
	my $allele_inference = shift;

	die unless(length($allele_validation) >= 4); die unless($allele_validation =~ /^\d\d/);
	$allele_validation = substr($allele_validation, 0, 2).':'.substr($allele_validation, 2);
	
	my $true_allele = $allele_validation;
	$true_allele =~ s/g//;
		
	my $inferred_alleles = $allele_inference;
	my @inferred_alleles = split(/;/, $inferred_alleles);
	
	@inferred_alleles = map {
		die "Can't parse allele $_" unless($_ =~ /^\s*\w+\*(\d\d)\:(\d\dN?).*$/);
		my $allele = $1.':'.$2;
		$allele
	} @inferred_alleles;
	
	my %inferred = map {$_ => 1} @inferred_alleles;
	
	if($inferred{$true_allele})
	{
		# print Dumper($allele_validation, $allele_inference, \@inferred_alleles, 1), "\n" if ($locus eq 'DQB1');	
		
		return 1;  
	}		
	else
	{
		# print Dumper($allele_validation, $allele_inference, \@inferred_alleles, 0), "\n"  if ($locus eq 'DQB1');	
	
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
		my @f = split(/\t/, $line);
		die unless($#f == 2);
		my $exon = $f[0];
		my $exonPos = $f[1];
		my $pileUp = $f[2];
		
		$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		
	}	
	close(F);
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
		$validation_allele =~ s/g//;
		die unless(length($validation_allele) >= 4);
		
		$validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
		$validation_allele = $locus . '*' . $validation_allele;
		
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

sub twoClusterAlleles
{
	my $alleles_inference = shift;
	die unless($#{$alleles_inference} == 1);
	
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
