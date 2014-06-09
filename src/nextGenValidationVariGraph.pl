#!/usr/bin/perl -w

use strict;
use 5.010;
use List::Util qw/min max/;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;
use File::Path qw/mkpath rmtree/;

my $kMer_size = 31;  

# input parameters

my $graph;
my $sample;
my $vcfPos;
my $VCF_for_comparison = '';
my $classical_VCF;
my $contigs_file;
my $cluster3 = 0;

my $localOxford = 0;
if(hostname() =~ /(sequoia)|(elm)|(birch)|(banyan)|(cluster3)/) 
{
	$localOxford = 1;
}

# Paths

my $original_alignment = '../data/mhc_ref_8_haplotypes/alignment/all_aligned.fasta';
my $referenceGenome =  '../data/GRCh37.60/fasta/Homo_sapiens.GRCh37.60.dna.chromosome.6.fa';
my $sample_path = qq(../data/samples/CortexGraphs);
my $path_to_PGF_haplotype = qq(../data/mhc_ref_8_haplotypes/pgf_ref.fasta);

if($localOxford)
{
	$original_alignment = '/gpfs1/well/gsk_hla/shared/mhc_ref_8_haplotypes/alignment/all_aligned.fasta';
	$referenceGenome =  '/gpfs1/well/gsk_hla/GRCh37.60/fasta/Homo_sapiens.GRCh37.60.dna.chromosome.6.fa';
	$sample_path = qq(/gpfs1/well/gsk_hla/CortexGraphs/);
}

# stop modifications here

GetOptions (
	'graph:s' => \$graph,
	'sample:s' => \$sample, 
	'kmer:s' => \$kMer_size,
	'vcfPos:s' => \$vcfPos,
	'original_alignment:s' => \$original_alignment,
	'referenceGenome:s' => \$referenceGenome,	
	'classical_VCF:s' => \$classical_VCF,
	'contigs_file:s' => \$contigs_file,
	'cluster3:s' => \$cluster3,
);         
 
if($cluster3)
{
	unless($localOxford)
	{
		die "You activated switch --cluster3, which parallelizes contig alignment but will only work in Oxford microenvironment - if you want to implement parallelization on your system, manually modify the code activated by the cluster3 switch (not too difficult).";
	}
}
die "Alignment file $original_alignment not existing" unless (-e $original_alignment);
die "No classical VCF file specified or file not existing" unless (-e $classical_VCF);
die "Reference genome $referenceGenome not existing" unless(-e $referenceGenome);
die "Please provide file with contigs!\n" unless(-e $contigs_file);

## Haplotype graph

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);
my $expected_normal_graph = $graph_dir.'/graph.txt';
die "Normal graph $expected_normal_graph not there!" unless(-e $expected_normal_graph);
my $expected_kMer_graph_file = $expected_normal_graph.".kmers_".$kMer_size;
die "kMerified graph $expected_kMer_graph_file not there" unless(-e $expected_kMer_graph_file);
my $expected_segments_file = $graph_dir.'/segments.txt';
die "Segments file $expected_segments_file not there" unless(-e $expected_segments_file);

## de Bruijn graph

my $sample_deBruijn_graph_file = $sample_path.'/'.$sample.'_'.$kMer_size.'.ctx';
die "de Bruijn graph $sample_deBruijn_graph_file not there" unless(-e $sample_deBruijn_graph_file);

## Coordinate bookkeeping & get data

my $xMHC_reference = $path_to_PGF_haplotype;
die "Cannot access $path_to_PGF_haplotype" unless (-e $path_to_PGF_haplotype);
 
my $pgf_start = 28702185;
my $remapping_flank_length = 100;

## get data on original alignment

my $graph_loci_data;
my @loci_in_graph;
my %pgf_to_graphAlignment;
my %PGF_within_HLA;
my $original_alignment_complete;
my $prefix_for_cache = 'temp/validation_temp_'.$graph.'_';
if(-e $prefix_for_cache.'graph_loci.data')
{
	warn "Use cached data ${prefix_for_cache}*! If you don't want this, delete files from temp/!";
	
	$graph_loci_data = retrieve $prefix_for_cache.'graph_loci.data';
	@loci_in_graph = @{retrieve $prefix_for_cache.'loci_in_graph.data'};
	%pgf_to_graphAlignment = %{retrieve $prefix_for_cache.'pgf_to_graphAlignment.data'};
	%PGF_within_HLA = %{retrieve $prefix_for_cache.'PGF_within_HLA.data'};
	$original_alignment_complete = retrieve $prefix_for_cache.'original_alignment_complete.data';
	
	print "\tloading data done\n";
}
else
{
	$graph_loci_data = get_graph_loci($graph_dir);
	@loci_in_graph = @{$graph_loci_data->[0]};

	# create index into graph PGF positions

	my $last_used_pgf = -1;
	my $first_position_PGF = -1;
	my $last_position_PGF = -1;
	foreach my $segment_file (keys %{$graph_loci_data->[1]})
	{
		if($segment_file =~ /HLA\w+\.txt/)
		{
			foreach my $f (@{$graph_loci_data->[1]{$segment_file}})
			{
				my @parts = split(/_/, $f);
				die unless($#parts == 2);
				$PGF_within_HLA{$parts[2]} = 1;
			}
		}
	}
	  
	for(my $i = 0; $i <= $#loci_in_graph; $i++)
	{
		my @parts = split(/_/, $loci_in_graph[$i]);
		my $pgfPos;
		
		if($#parts == 2)
		{
			$pgfPos = $parts[2];	
			if($i == 0)
			{
				$first_position_PGF	= $pgfPos;
			}
			if($i == $#loci_in_graph)
			{
				$last_position_PGF	= $pgfPos;		
			}
		}
		else
		{
			die $loci_in_graph[$i];
			die unless($last_used_pgf > 0);
			$pgfPos = $last_used_pgf + 1;
		}
		
		if(not exists $pgf_to_graphAlignment{$pgfPos})
		{
			$pgf_to_graphAlignment{$pgfPos} = [];
		}
		push(@{$pgf_to_graphAlignment{$pgfPos}}, $i);
		$last_used_pgf = $pgfPos;	
	}
	die unless($first_position_PGF != -1);
	die unless($last_position_PGF != -1);

	# die Dumper($pgf_to_graphAlignment{1877}, $pgf_to_graphAlignment{1878}, $pgf_to_graphAlignment{1879});

	print "Extracting from graph: $graph\n";
	print "Original alignment: $original_alignment\n";

	$original_alignment_complete = read_alignment($original_alignment);

	store $graph_loci_data, $prefix_for_cache.'graph_loci.data';
	store \@loci_in_graph, $prefix_for_cache.'loci_in_graph.data';
	store \%pgf_to_graphAlignment, $prefix_for_cache.'pgf_to_graphAlignment.data';
	store \%PGF_within_HLA, $prefix_for_cache.'PGF_within_HLA.data';
	store $original_alignment_complete, $prefix_for_cache.'original_alignment_complete.data';	
}

## output directory for aligned contigs
my $graphForFileName = $graph_dir;
$graphForFileName =~ s/^.+tmp2//;
$graphForFileName =~ s/\W/_/g;

my $contig_general_dir = qq(../tmp/alignedContigs/);
unless(-e $contig_general_dir)
{
	mkdir($contig_general_dir) or die "Cannot mkdir $contig_general_dir";
}
my $contig_dir_graph = $contig_general_dir . join('_', $graphForFileName, $sample, $kMer_size);
unless(-e $contig_dir_graph)
{
	mkdir($contig_dir_graph) or die "Cannot mkdir $contig_dir_graph";
}
my $kMer_validation_output_dir = $contig_dir_graph . '/kMerValidation';
unless(-e $kMer_validation_output_dir)
{
	mkdir($kMer_validation_output_dir) or die "Cannot mkdir $kMer_validation_output_dir";
}
die unless((-e $contig_general_dir) and (-e $contig_dir_graph) and (-e $kMer_validation_output_dir));
my $contigFile_noSpecialCharacters = (fileparse($contigs_file))[0];
$contigFile_noSpecialCharacters =~ s/\W/_/g;
my $contig_specific_dir = $contig_dir_graph . '/' . $contigFile_noSpecialCharacters;
unless(-e $contig_specific_dir)
{
	mkdir($contig_specific_dir) or die "Cannot mkdir $contig_specific_dir";
	mkdir($contig_specific_dir.'/toReference') or die "Cannot mkdir a dir";
	mkdir($contig_specific_dir.'/toViterbiChromotypes') or die "Cannot mkdir a dir";
	mkdir($contig_specific_dir.'/toAmendedChromotypes') or die "Cannot mkdir a dir";
	mkdir($contig_specific_dir.'/toVCF') or die "Cannot mkdir a dir";
}

## read chromotypes

my $kMer_count_sample_required = qq(../tmp/kMerCount_).join('_', $graphForFileName, $sample, $kMer_size,  'required');
my $kMer_count_sample = $kMer_count_sample_required.'.binaryCount';
my $expected_output_haplotypes = $kMer_count_sample.'.viterbiHaplotypes';
my $expected_output_chromotypes = $kMer_count_sample.'.viterbiGenomeString';
my $haplotypes_aref = read_haplotypes($expected_output_haplotypes, \@loci_in_graph);
my $expected_chromotype_length = length($haplotypes_aref->[0][0]);
	
my @chromotypes_withTrailingTwoBuffers;
my $chromotype_length = 0;
open(CHROMOTYPES, '<', $expected_output_chromotypes) or die "Cannot open $expected_output_chromotypes";
while(<CHROMOTYPES>)
{
	my $line = $_;
	chomp($line);
	my @parts = split(/,/, $line, -1);
	die unless($#parts == 1);
	my $compartment = [$parts[0]];
	if(length($parts[1]) > 0)
	{
		push(@$compartment, $parts[1]);
		die unless(length($parts[1]) == length($parts[0]));
	}
	$chromotype_length += length($parts[0]);
	push(@chromotypes_withTrailingTwoBuffers, $compartment);
}
close(CHROMOTYPES);
$chromotype_length -= 2;
die "Wrong chromotype length: expect $expected_chromotype_length, but got $chromotype_length" unless($chromotype_length == $expected_chromotype_length);

# we need to find the right coordinates for the comparison

my $classical_VCF_xMHC_boundaries = VCF_get_xMHC_minmax($classical_VCF);
my $classical_VCF_min_xMHC = $classical_VCF_xMHC_boundaries->[0];
my $classical_VCF_max_xMHC = $classical_VCF_xMHC_boundaries->[1];
my $chromotypes_min_xMHC = min(keys %pgf_to_graphAlignment);
my $chromotypes_max_xMHC = max(keys %pgf_to_graphAlignment);
my $chromotypes_min_xMHC_inChromotypes = $pgf_to_graphAlignment{$chromotypes_min_xMHC}[0];
my $chromotypes_max_xMHC_inChromotypes = $pgf_to_graphAlignment{$chromotypes_max_xMHC}[-1];

print "Determined the following available coordinates:\n\tclassical VCF: $classical_VCF_min_xMHC - $classical_VCF_max_xMHC\n\tchromotypes: $chromotypes_min_xMHC - $chromotypes_max_xMHC\n\t\tchromotype coordinates: $chromotypes_min_xMHC_inChromotypes - $chromotypes_max_xMHC_inChromotypes\n\n";

my $consensus_min_PGF = ($classical_VCF_min_xMHC > $chromotypes_min_xMHC) ? $classical_VCF_min_xMHC : $chromotypes_min_xMHC;
my $consensus_max_PGF = ($classical_VCF_max_xMHC > $chromotypes_max_xMHC) ? $chromotypes_max_xMHC : $classical_VCF_max_xMHC;
die unless($consensus_max_PGF > $consensus_min_PGF);

my $consensus_min_genomic = $consensus_min_PGF + $pgf_start;
my $consensus_max_genomic = $consensus_max_PGF + $pgf_start;

my $consensus_min_inChromotypes = $pgf_to_graphAlignment{$consensus_min_PGF}[0];
my $consensus_max_inChromotypes = $pgf_to_graphAlignment{$consensus_max_PGF}[-1];
die unless($consensus_max_inChromotypes > $consensus_min_inChromotypes);

print "Consensus coordinates:\n\tPGF: $consensus_min_PGF - $consensus_max_PGF\n\tGenomic (for VCF restriction): $consensus_min_genomic - $consensus_max_genomic \n\tChromotype coordinates (for chromotype restriction): $consensus_min_inChromotypes - $consensus_max_inChromotypes\n";

my $classical_VCF_restricted = $classical_VCF.'.validationRestricted';
restrictVCF($classical_VCF, $classical_VCF_restricted, '6', $consensus_min_PGF+$pgf_start, $consensus_max_PGF+$pgf_start);

if($cluster3)
{
	# copy over essential files - graph etc
	
	my $cluster3_expected_kMer_graph_file = $expected_kMer_graph_file;
	$cluster3_expected_kMer_graph_file =~ s/^\.\./\/gpfs1\/well\/gsk_hla\/MHC-PRG/;
	
	my $cluster3_kMer_count_sample = $kMer_count_sample;
	$cluster3_kMer_count_sample =~ s/^\.\./\/gpfs1\/well\/gsk_hla\/MHC-PRG/;
	
	my $copyOver = sub {
		my $source = shift;
		my $target = shift;
		
		my $target_dir = dirname($target);
		unless(-e $target_dir)
		{
			mkpath($target_dir) or die "Cannot mkdir $target_dir";
		}	
		
		copy($source, $target) or die "Cannot copy $source to $target";
	};
	
	$copyOver->($expected_kMer_graph_file, $cluster3_expected_kMer_graph_file);
	$copyOver->($kMer_count_sample, $cluster3_kMer_count_sample);
	
	foreach my $f (glob($kMer_count_sample.'*'))
	{
		my $f_cluster3 = $f;
		$f_cluster3 =~ s/^\.\./\/gpfs1\/well\/gsk_hla\/MHC-PRG/;
		# print $f, " => ", $f_cluster3, "\n";
		$copyOver->($f, $f_cluster3);		
	}
	
	die unless($classical_VCF_restricted =~ /^\/gpfs1/);
	die unless($referenceGenome =~ /^\/gpfs1/);
	die unless($sample_deBruijn_graph_file =~ /^\/gpfs1/);
		
	my $cluster3_contig_general_dir = qq(/gpfs1/well/gsk_hla/MHC-PRG/tmp/alignedContigs/);
	unless(-e $cluster3_contig_general_dir)
	{
		mkdir($cluster3_contig_general_dir) or die "Cannot mkdir $cluster3_contig_general_dir";
	}
	my $cluster3_contig_dir_graph = $cluster3_contig_general_dir . join('_', $graphForFileName, $sample, $kMer_size) . '_P';
	unless(-e $cluster3_contig_dir_graph)
	{
		mkdir($cluster3_contig_dir_graph) or die "Cannot mkdir $cluster3_contig_dir_graph";
	}
	
	my $contigFile_noSpecialCharacters = (fileparse($contigs_file))[0];
	$contigFile_noSpecialCharacters =~ s/\W/_/g;
	
	my $cluster3_contig_dir_for_contigsFile = $cluster3_contig_dir_graph . '/' . $contigFile_noSpecialCharacters;
	
	my $bases_per_job = 100000;
	
	if($cluster3 eq 'prepare')
	{		
		if(-e $cluster3_contig_dir_for_contigsFile)
		{
			warn "Directory " . $cluster3_contig_dir_for_contigsFile . " existing already, i.e. cluster 3 preparation has already taken place!\n";
		}
		mkdir($cluster3_contig_dir_for_contigsFile);
		unless(-e $cluster3_contig_dir_for_contigsFile)   
		{
			die "Cannot mkdir $cluster3_contig_dir_for_contigsFile";		
		}
		
		my $job_submission_file = $cluster3_contig_dir_for_contigsFile . '/submissions.txt';
		open(SUBMISSIONS, '>', $job_submission_file) or die "Cannot open $job_submission_file";
		
		my $current_fh_contig_printing;
		
		my $prepareJobDirectory = sub {
			my $threadN = shift;
			my $cluster3_contig_specific_dir = $cluster3_contig_dir_for_contigsFile . '/job_'.$threadN;
			unless(-e $cluster3_contig_specific_dir)
			{
				mkdir($cluster3_contig_specific_dir) or die "Cannot mkdir $cluster3_contig_specific_dir";
				mkdir($cluster3_contig_specific_dir.'/toReference') or die "Cannot mkdir a dir";
				mkdir($cluster3_contig_specific_dir.'/toViterbiChromotypes') or die "Cannot mkdir a dir";
				mkdir($cluster3_contig_specific_dir.'/toAmendedChromotypes') or die "Cannot mkdir a dir";
				mkdir($cluster3_contig_specific_dir.'/toVCF') or die "Cannot mkdir a dir";
			}
			
			if($current_fh_contig_printing)
			{
				close($current_fh_contig_printing);
			}
			
			my $contigs_for_job_fN = $cluster3_contig_specific_dir . '/' . 'contigs.txt';
			open($current_fh_contig_printing, '>', $contigs_for_job_fN) or die "Cannot open $contigs_for_job_fN";		
		
			my $command_align_chromotypes = qq(/gpfs1/well/gsk_hla/MHC-PRG/bin/MHC-PRG domode nextGenContigValidation $cluster3_expected_kMer_graph_file $cluster3_kMer_count_sample $consensus_min_inChromotypes $consensus_max_inChromotypes $classical_VCF_restricted $consensus_min_genomic $consensus_max_genomic $referenceGenome $sample_deBruijn_graph_file $cluster3_contig_specific_dir $contigs_for_job_fN);
	
			my $filename_for_qsub = $cluster3_contig_specific_dir . '/qsub.txt';
			open(QSUB, '>', $filename_for_qsub) or die "Cannot open $filename_for_qsub";
print QSUB qq(#!/bin/bash
#\$ -o $cluster3_contig_specific_dir
#\$ -P mcvean.prjb -q short.qb
#\$ -pe shmem 6
$command_align_chromotypes
);
			close(QSUB);
			
			print SUBMISSIONS qq(qsub $filename_for_qsub), "\n";
		};
		
		my $currentJob = 1;
		$prepareJobDirectory->($currentJob);
		
		my $currrent_job_bases = 0;
		open(CONTIGS, '<', $contigs_file) or die "Cannot open $contigs_file";
		while(<CONTIGS>)
		{
			my $contigID = $_;
			my $contigSequence = <CONTIGS>;
			
			print {$current_fh_contig_printing} $contigID, $contigSequence;
			
			$currrent_job_bases += length($contigSequence);
			
			if($currrent_job_bases > $bases_per_job)
			{
				$currentJob++;
				$prepareJobDirectory->($currentJob);				
				$currrent_job_bases = 0;  
			}
		}
		close(CONTIGS);
		close(SUBMISSIONS);	
		
		print "\n\nCluster3 splitting done, now execute commands in $job_submission_file \n\n";
	}
	elsif($cluster3 eq 'collect')
	{
		my $contig_count = 0;
		open(INPUT_CONTIGS, '<', $contigs_file) or die "Cannot open $contigs_file";
		while(<INPUT_CONTIGS>)
		{
			my $line = $_;
			if(substr($line, 0, 1) eq '>')
			{
				$contig_count++;
			}
		}
		close(INPUT_CONTIGS);
		
		my @methods = qw/toVCF toReference toViterbiChromotypes toAmendedChromotypes/;
	
		my @subdirectories = glob($cluster3_contig_dir_for_contigsFile.'/*');
		@subdirectories = grep {(-d $_) and ($_ =~ /job_/)} @subdirectories;
		
		my %have_contigs_per_method;
		my %resubmit;
		
		foreach my $method (@methods)
		{
			my $output_dir = $contig_specific_dir . '/' . $method;
			if(-e $output_dir)
			{
				warn "$output_dir existing, delete and re-create!\n";
				rmtree($output_dir);
				mkpath($output_dir) or die "Cannot mkpath $output_dir";
			}
		}
		
		foreach my $jobDirectory (@subdirectories)
		{		
			die unless($jobDirectory =~ /job_(\d+)/);
			my $jobNumber = $1;
					
			my $qsub_file = $jobDirectory.'/qsub.txt';
			die unless(-e $qsub_file);
			
			my %output_FHs = map {
				my $fN = $contig_specific_dir . '/' . $_ . '/' . $jobNumber . '.alignment';
				my $fH;
				open($fH, '>', $fN) or die "Cannot open $fN";
				$_ => $fH 
			} @methods;
			
			foreach my $method (@methods)
			{
				my $method_directory = $jobDirectory . '/'. $method;
				my $statusFile = $method_directory . '/all.status';
				my $statusOK = 1;
				if(-e $statusFile)
				{
					open(STATUS, '<', $statusFile) or die "Cannot open $statusFile";
					my $status = <STATUS>;
					$status = substr($status, 0, 1);
					unless($status eq '1')
					{
						$statusOK = 0;
					}
					close(STATUS);
				}
				else
				{
					$statusOK = 0;
				}
				
				unless($statusOK)
				{
					$resubmit{$qsub_file}++;
					warn "Job $jobNumber -- method $method not completed yet - ignore!\n";
					next;
				}
				
				my @files = glob($method_directory . '/*.alignment');
				@files = grep {$_ !~ /\.status/} @files;
				foreach my $file (@files)
				{
					open(F, '<', $file) or die "Cannot open $file";
					while(<F>)
					{
						my $line = $_;
						die unless($output_FHs{$method});
						print {$output_FHs{$method}} $line;
						if(substr($line, 0, 1) eq '>')
						{
							$have_contigs_per_method{$method}++;
						}
					}
					close(F);
				}
			}
			
			foreach my $method (keys %output_FHs)
			{
				close($output_FHs{$method});
			}			
		}
		
		print "Method collection statistics (first 2 x input contig count, then contig-count per method):\n";
		foreach my $method (@methods)
		{
			my $statusFile = $contig_specific_dir . '/' . $method . '/all.status';
			open(STATUS, '>', $statusFile) or die "Cannot open $statusFile";
			print STATUS 1;
			close(STATUS);
			
			my $expect_output_contigs = 5 * $contig_count;
			
			print "\t", join("\t", $method, $expect_output_contigs, $have_contigs_per_method{$method}), "\n";
			if($have_contigs_per_method{$method} != $expect_output_contigs)
			{
				warn "\n\nMethod ${method}: would like to see $contig_count x 5, but have $have_contigs_per_method{$method} ! \n\n";
			}
		}
		
		print "Completed cluster3 collection - now re-run, without --cluster3!\n";
		
		if(keys %resubmit)
		{
			my $resubmit_fn = '_resubmit.txt';
			open(RESUBMIT, '>', $resubmit_fn) or die "Cannot open $resubmit_fn";
			print RESUBMIT join ("\n", map {'qsub '.$_} keys %resubmit), "\n";   
			close(RESUBMIT);
			print "\n\nIf you want to resubmit, execute commands in file _resubmit.txt\n\n";
		}
	}
	else
	{
		die "Unknown value for argument --cluster3";
	}	
}
else
{
	my $command = qq(../bin/MHC-PRG domode nextGenValidation $expected_kMer_graph_file $kMer_count_sample $consensus_min_inChromotypes $consensus_max_inChromotypes $classical_VCF_restricted $consensus_min_genomic $consensus_max_genomic $referenceGenome $sample_deBruijn_graph_file $kMer_validation_output_dir $contigs_file $graph_dir);

	print "Execute:\n\n$command\n\n";

	my $command_align_chromotypes = qq(../bin/MHC-PRG domode nextGenContigValidation $expected_kMer_graph_file $kMer_count_sample $consensus_min_inChromotypes $consensus_max_inChromotypes $classical_VCF_restricted $consensus_min_genomic $consensus_max_genomic $referenceGenome $sample_deBruijn_graph_file $contig_specific_dir $contigs_file $graph_dir);

	print "And execute this for contig validation:\n\n";

	print $command_align_chromotypes;

	print "\n\n";
}

sub get_graph_loci
{
	my $graph_dir = shift;
	
	my %locus_origin;
	
	# die $graph_dir;
	
	my @loci_in_graph;
	my @locus_alleles;
	
	my $graph_segments_file = $graph_dir.'/segments.txt';
	die "File not there: $graph_segments_file" unless (-e $graph_segments_file);
	open(SEGMENTS, '<', $graph_segments_file) or die;
	while(<SEGMENTS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
		next unless($line);
		
		my $segment_file = $graph_dir.'/'.$line;

		open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
		my $firstLine = <SEGMENT>;
		chomp($firstLine);
		$firstLine =~ s/\n//g;
		$firstLine =~ s/\r//g;			
		my @line_fields = split(/ /, $firstLine);
		shift(@line_fields); # kick out individual ID
		
		push(@{$locus_origin{$segment_file}}, @line_fields);
		push(@loci_in_graph, @line_fields);
		
		my $beginning_of_locus_alleles = $#locus_alleles + 1;
		push(@locus_alleles, map {{}} @line_fields);
		
		while(<SEGMENT>)
		{
			my $line = $_;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;			
			my @this_line_fields = split(/ /, $line);
			shift(@this_line_fields); # kick out individual ID		
			for(my $j = 0; $j <= $#this_line_fields; $j++)
			{
				if($j == 0)
				{
					# print join(' ', $segment_file, $beginning_of_locus_alleles, $j, $beginning_of_locus_alleles + $j, $this_line_fields[$j]), "\n";
				}
				$locus_alleles[$beginning_of_locus_alleles + $j]{$this_line_fields[$j]}++;
			}
		}
		close(SEGMENT);
	}
	close(SEGMENTS);
	
	return [\@loci_in_graph, \%locus_origin, \@locus_alleles];
}


sub read_alignment
{
	# this is code duplicated from documents\analysis\04 Januar 2012\MHC alignment prepare for variation graph 2.pl , and should stay
	# consistent with that! (apart from the bit which determines alignment positions)
	
	my $file = shift;

	my %alignment;
	my %alignment_positions;
	unless (-e $file)
	{
		$file = 'all_aligned.fasta';
	}
	
	open(ALIGNMENT, "<", $file) or die "Cannot open $file";
	my $current_name;
	while(<ALIGNMENT>)
	{
		my $line = $_;
		chomp($line);

		if($line =~ /^>/)
		{
			if($line =~ /chr6_(\w+?)_/)
			{
				$current_name = $1;
			}
			else
			{
				$current_name = 'pgf';
			}
			die if(exists $alignment{$current_name});
		}
		else
		{
			$line =~ tr/acgtn/ACGTN/;
			die unless($current_name);
			$line =~ s/\./_/g;
			$line =~ s/\-/_/g;
			if($current_name eq 'pgf')
			{
				if($line =~ /n|N/)
				{
					die $line;
				}
			}
			
			$line =~ s/N/\*/g;			
			$alignment{$current_name} .= $line;
			die unless ($line =~ /^[ACGT\*\_]+$/);
		}
	}
	close(ALIGNMENT);
	
	foreach my $haplotype (keys %alignment)
	{
		my $sequence = $alignment{$haplotype};
		my $nonGapIndex = 0;
		for(my $i = 0; $i < length($sequence); $i++)
		{
			my $char = substr($sequence, $i, 1);
			die unless($char);
			if($char ne '_')
			{
				$nonGapIndex++;
			}
			push(@{$alignment_positions{$haplotype}}, $nonGapIndex);
		}
		
		die unless(scalar(@{$alignment_positions{$haplotype}}) == length($sequence));
	}
	
	
	# remove all gaps
	
	my @positions_gap_count;
	my $alignment_length = -1;	
	foreach my $haploID (keys %alignment)
	{
		my $al = $alignment{$haploID};
		
		if($alignment_length == -1)
		{
			$alignment_length = length($al);
		}
		
		die unless(length($al) == $alignment_length);
		
		for(my $i = 0; $i < length($al); $i++)
		{
			my $c = substr($al, $i, 1);
			if($c eq '_')
			{	
				$positions_gap_count[$i]++;
			}
		}
	}
	my @positions_all_gaps = grep {((defined $positions_gap_count[$_]) ? $positions_gap_count[$_] : 0) == scalar(keys %alignment)} (0 .. $#positions_gap_count);
	my %positions_all_gaps = map {$_ => 1} @positions_all_gaps;
	print "Remove ", scalar(@positions_all_gaps), " positions from alignment because all gaps\n";
	
	$alignment_length = -1;
	foreach my $haploID (keys %alignment)
	{
		my $al_old = $alignment{$haploID};
		my $al_new = '';
		
		for(my $i = 0; $i < length($al_old); $i++)
		{
			my $c = substr($al_old, $i, 1);
			if(! $positions_all_gaps{$i})
			{
				$al_new .= $c;
			}
		}
		
		if($alignment_length == -1)
		{
			$alignment_length = length($al_new);
		}
		
		die unless(length($al_new) == $alignment_length);
		die unless ($al_new =~ /^[ACGT\*\_]+$/);
		$alignment{$haploID} = $al_new;
	}
	
	my $remove_round_2 = 0;
	my %new_alignment;
	my %new_alignment_positions;
	
	my $delete_positions_in_row = 0;
	for(my $i = 0; $i < $alignment_length; $i++)
	{
		my $gappish_haplos = 0;
		foreach my $haploID (keys %alignment)
		{
			my $c = substr($alignment{$haploID}, $i, 1);
			if(($c eq '*') or ($c eq '_'))
			{
				$gappish_haplos++;
			}
		}
		my $c_pgf = substr($alignment{'pgf'}, $i, 1);
		die unless($c_pgf);
		
		my $potentially_delete_position = 0;
		if(($gappish_haplos == scalar(keys %alignment)) and ($c_pgf eq '_'))
		{
			$potentially_delete_position = 1;
		}
		
		my $do_delete_position = 0;
		if($potentially_delete_position)
		{
			$delete_positions_in_row++;
			if($delete_positions_in_row > $kMer_size)
			{
				$do_delete_position = 1;
			}
		}
		else
		{
			$delete_positions_in_row = 0;
		}
		
		unless($do_delete_position)
		{
			foreach my $haploID (keys %alignment)
			{
				my $c = substr($alignment{$haploID}, $i, 1);
				$new_alignment{$haploID} .= $c;
				my $alignment_pos = $alignment_positions{$haploID}[$i];
				die unless(defined $alignment_pos);
				if(not $new_alignment_positions{$haploID})
				{
					$new_alignment_positions{$haploID} = [];
				}
				
				push(@{$new_alignment_positions{$haploID}}, $alignment_pos);
			}		
		}
		
		$remove_round_2 += $do_delete_position;
	}	
	
	print "Remove ", $remove_round_2, " positions from alignment because in sequence > length $kMer_size of gaps or stars, and PGF is gap\n";
	
	$alignment_length = -1;
	foreach my $haploID (keys %new_alignment)
	{
		if($alignment_length == -1)
		{
			$alignment_length = length($new_alignment{$haploID});
		}
		
		die unless(length($new_alignment{$haploID}) == $alignment_length);
	}
	
	return {alignment => \%new_alignment, positions => \%new_alignment_positions};
}

sub restrictVCF
{
	my $file = shift;
	my $output_restricted_file = shift;
	my $extraction_chromosome = shift;
	my $extraction_min = shift;
	my $extraction_max = shift;
	
	open(VCF, '<', $file) or die "Cannot open $file";
	open(VCFOUT, '>', $output_restricted_file) or die "Cannot open $output_restricted_file";
	
	while(<VCF>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 2) eq '##')
		{
			print VCFOUT $line, "\n";
			next;
		}
		if(substr($line, 0, 1) eq '#')
		{
			my @header_fields = split(/\t/, $line);
			die unless($header_fields[0] eq '#CHROM');
			die unless($header_fields[1] eq 'POS');
			print VCFOUT $line, "\n";			
			next;
		}   
		
		my @fields = split(/\t/, $line);
		
		my $chromosome = $fields[0];
		my $position = $fields[1];
		
		next unless($chromosome eq $extraction_chromosome);
		next unless(($position >= $extraction_min) and ($position <= $extraction_max));
		
		print VCFOUT $line, "\n";
	}
	close(VCF);
	close(VCFOUT);
}

sub VCF_get_xMHC_minmax
{
	my $file = shift;
	
	my $p_min = -1;
	my $p_max = -1;
	
	open(VCF, '<', $file) or die "Cannot open $file";

	while(<VCF>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 2) eq '##')
		{
			next;   
		}
		if(substr($line, 0, 1) eq '#')
		{
			my @header_fields = split(/\t/, $line);         
			die Dumper("Problem with header fields", $file, $header_fields[0], @header_fields) unless($header_fields[0] eq '#CHROM');
			die unless($header_fields[1] eq 'POS');
			next;
		}
		
		my @fields = split(/\t/, $line);
		
		my $chromosome = $fields[0];
		my $position = $fields[1];
		
		next unless($chromosome eq '6');
		if(($p_min == -1) or ($position < $p_min))
		{
			$p_min = $position;
		}
		
		if(($p_max == -1) or ($position > $p_max))
		{
			$p_max = $position;
		}		
	}
	close(VCF);
	
	if(($p_min == -1) or ($p_max == -1))
	{
		die "Problems determining the min and max xMHC values from VCF $file -- expect chromosome number identifers, not prefix 'chr'";
	}
	
	return [$p_min - $pgf_start, $p_max - $pgf_start];
}


sub read_haplotypes
{
	my $haplotypes_file = shift;
	my $loci_in_graph_aref = shift;
	my $amendedHaplotypes = shift;
	
	open(HAPLOTYPES, '<', $haplotypes_file) or die "Cannot open $haplotypes_file";
	my $haplotype_1, my $haplotype_2;
	my @haplotype_1_alignment_fields;
	my @haplotype_2_alignment_fields;
	while(<HAPLOTYPES>)
	{
		my $line = $_;
		chomp($line);
		last if ($. > 2);
		my @fields = split(/ /, $line);
		my $lastIndex = $#fields - 2;
		if($amendedHaplotypes)
		{
			$lastIndex = $#fields;
		}
		my $haplotype = join('', @fields[2 .. $lastIndex]);
		if($. == 1)
		{
			$haplotype_1 = $haplotype;
			@haplotype_1_alignment_fields = @fields[2 .. $lastIndex];
		}
		elsif($. == 2)
		{
			$haplotype_2 = $haplotype;	
			@haplotype_2_alignment_fields = @fields[2 .. $lastIndex];
		}
	}
	close(HAPLOTYPES);

	die unless($haplotype_1 and $haplotype_2);
	die Dumper('$#haplotype_1_alignment_fields == $#{$loci_in_graph_aref}', $haplotypes_file, $#haplotype_1_alignment_fields, $#{$loci_in_graph_aref}) unless($#haplotype_1_alignment_fields == $#{$loci_in_graph_aref});
	die Dumper('$#haplotype_2_alignment_fields == $#{$loci_in_graph_aref}', $haplotypes_file, $#haplotype_2_alignment_fields, $#{$loci_in_graph_aref})  unless($#haplotype_2_alignment_fields == $#{$loci_in_graph_aref});
	
	return [[$haplotype_1, $haplotype_2], [\@haplotype_1_alignment_fields, \@haplotype_2_alignment_fields]];
}

