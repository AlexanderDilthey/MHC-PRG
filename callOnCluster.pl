#!/usr/bin/perl -w

use strict;
use warnings;
use List::MoreUtils qw/all mesh any /; 
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;  
use Storable;
use FindBin;
use File::Spec;

# set so we get at least ~80G of RAM
my $scheduling_threads = 12;
my $scheduling_P = 'mcvean.prjb';
my $scheduling_q = 'short.qb';

my $action;
my $graph = 'hla';   
my $BAMDir;
my $HiSeq250bp = 0;
my $MiSeq250bp = 0;
my $referenceGenome = '/gpfs1/well/gsk_hla/GRCh37.60/fasta/combined/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.fa';
my $executeLocally = 0;
my $removeExistingFiles = 0;
my $aggregateOutput;
my $sampleIDPrefix = '';
GetOptions (
 'action:s' => \$action,
 'graph:s' => \$graph,
 'BAMDir:s' => \$BAMDir,
 'HiSeq250bp:s' => \$HiSeq250bp, 
 'referenceGenome:s' => \$referenceGenome, 
 'executeLocally:s' => \$executeLocally, 
 'aggregateOutput:s' => \$aggregateOutput, 
 'sampleIDPrefix:s' => \$sampleIDPrefix, 
);         

die "Please set argument --referenceGenome" unless(-e $referenceGenome);

unless($action)
{
	die "Please specify action - e.g. qsub, check, aggregate";
}

unless($BAMDir)
{
	die "Please specify --BAMDir";
}

my @BAMs = glob($BAMDir . '/*.bam');
my @sampleIDs;
foreach my $BAM (@BAMs)
{
	my $bai = $BAM . '.bai';
	unless(-e $bai)
	{
		die "BAM $bai not indexed";
	}
	
	die "Can't parse BAM file path $BAM" unless($BAM =~ /.+[\\\/](\S+?)\.bam$/);	
	my $sampleID = $sampleIDPrefix . $1;
	push(@sampleIDs, $sampleID);
}

unless((-e 'temp_qsub') and (-d 'temp_qsub'))
{
	mkdir('temp_qsub') or die "Cannot mkdir temp_qsub";
}
	
die unless($#BAMs == $#sampleIDs);
my %_BAMs = map {$_ => 1} @BAMs;
my %_sampleIDs = map {$_ => 1} @sampleIDs;
die Dumper("Non-unique BAMs", \@BAMs) unless(scalar(keys %_BAMs) == scalar(@BAMs));
die Dumper("Non-unique sample IDs", \@sampleIDs) unless(scalar(keys %_sampleIDs) == scalar(@sampleIDs));

my $aggregate_output_fh;
if($action eq 'aggregate')
{
	die "Please specify --aggregateOutput" unless(defined $aggregateOutput);
	open($aggregate_output_fh, '>', $aggregateOutput) or die "Cannot open $aggregateOutput";
}
my $this_bin_dir = $FindBin::RealBin;
my %qsub_files;
for(my $BAMi = 0; $BAMi < $#BAMs; $BAMi++)
{
	my $BAM = $BAMs[$BAMi];
	my $sampleID = $sampleIDs[$BAMi];
	my $working_dir = '../tmp/hla/'.$sampleID;
	my $results_file = $working_dir . '/R1_bestguess.txt';;
	if($action eq 'qsub')
	{
		if(-e $results_file)
		{
			if($removeExistingFiles)
			{
				unlink($results_file) or die "Cannot unlink $results_file";
			}
			else
			{
				warn "Results file $results_file exists, but don't delete because --removeExistingFiles is not set";
			}
		}
		
		my $qsub_file = $sampleID;
		$qsub_file =~ s/\W//g;
		$qsub_file = 'temp_qsub/' . $qsub_file;
		
		if(exists $qsub_files{$qsub_file})
		{
			die "Automatically generated qsub file name $qsub_file more than once, use different sample IDs";
		}
		$qsub_files{$qsub_file} = 1;
		
		die "File $BAM not present" unless(-e $BAM);
		my $BAM_abs = File::Spec->rel2abs($BAM);
		my $qsub_time = $qsub_file . ".output_and_timing";
		
		open(QSUB, ">", $qsub_file) or die "Cannot open $qsub_file";
		print QSUB qq(#!/bin/bash
#\$ -P $scheduling_P  -q $scheduling_q
#\$ -pe shmem $scheduling_threads
source ~/.profile
cd $this_bin_dir
/usr/bin/time -v ./HLAtypeinference.pl --actions pnai --sampleIDs $sampleID --BAMs $BAM_abs --referenceGenome $referenceGenome --HiSeq250bp $HiSeq250bp --MiSeq250bp $MiSeq250bp --threads $scheduling_threads --graph $graph &> $qsub_time
	);
	close(QSUB);
	}
	elsif($action eq 'check')
	{
		if(-e $results_file)
		{
			print join("\t", $sampleID, "OK"), "\n";
		}
		else
		{
			print join("\t", $sampleID, "NOT FOUND", "!"), "\n";
		
		}
	}
	elsif($action eq 'aggregate')
	{
		die unless(defined $aggregate_output_fh);
		if(not -e $results_file)
		{
			warn "File for $sampleID not present, skip aggregation";
			next;
		}		
	}
	else
	{
		die "Unknown action: $action";
	}
}

foreach my $f (keys %qsub_files)
{
	if($executeLocally)
	{
		my $cmd = "bash $f";
		print $cmd, "\n";		
		if(system($cmd))
		{
			die "Command $cmd failed";
		}
	}
	else
	{
		my $cmd = "qsub $f";
		print $cmd, "\n";
		if(system($cmd))
		{
			die "Command $cmd failed";
		}	
	}
}
