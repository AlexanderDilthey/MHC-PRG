#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;
  
my $kMer_size = 55;  

# input parameters

my $graph = 'hla';   
my $sampleIDs;
my $BAMs;
my $actions;

GetOptions ('graph:s' => \$graph,
 'sampleIDs:s' => \$sampleIDs, 
 'BAMs:s' => \$BAMs, 
 'actions:s' => \$actions, 
);         


my $genome_graph_file = qq(../tmp2/GS_nextGen/hla/derived/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph.ctx);
unless(-e $genome_graph_file)
{
	die "Please set variable \$genome_graph_file to an existing file - the current value $genome_graph_file is not accessible.";
}
my $expected_kMer_file = qq(../tmp2/GS_nextGen/${graph}/requiredkMers_graph.txt.kmers_25);
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
		
		my $command = qq($use_bin domode HLATypeInference --input_alignedReads $aligned_file --graphDir ../tmp2/GS_nextGen/${graph} --sampleID $sampleID &> $stdout_file);
	
		print "Now executing command:\n$command\n\n";
		
		system($command);		
	}
}



