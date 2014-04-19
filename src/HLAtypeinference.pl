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

if($actions =~ /n/)
{
	unless(@BAMs)
	{
		die "Please provide --BAMs for negative filtering";
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
		
		my $output_file = '../tmp2/hla/'.$sampleI.'/reads.p';
		unless(-e '../tmp2/hla/'.$sampleI)
		{
			mkdir('../tmp2/hla/'.$sampleI) or die "Cannot mkdir ".'../tmp2/hla/'.$sampleI;
		}
		
		my $command = qq($use_bin domode filterReads --input_BAM $BAM --positiveFilter $expected_kMer_file --output_FASTQ $output_file);
	}
}