#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/all mesh any /; 
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;  
use Storable;
  
my $graph = 'hla';   
my $action = '';
my $iteration = 1;
my $sampleIDs;
my $BAMs; 
my $HiSeq250bp = 0;
my $threads = 1;

# my $referenceGenome = qq(/gpfs1/well/chimp/oa/ref/hs37d5.fasta);
my $referenceGenome = qq(/gpfs1/well/gsk_hla/GRCh37.60/fasta/combined/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.fa);

#$referenceGenome = '';

GetOptions ('graph:s' => \$graph,
 'action:s' => \$action, 
 'iteration:s' => \$iteration,  
 'sampleIDs:s' => \$sampleIDs,
 'BAMs:s' => \$BAMs,
 'HiSeq250bp:s' => \$HiSeq250bp, 
);         

my $target_directory_for_copying = qq(/Net/birch/data/dilthey/MHC-PRG/tmp/hla);

my @BAMs;
my @sampleIDs;
if($sampleIDs)
{
	@sampleIDs = split(/,/, $sampleIDs);
	if($BAMs)
	{
		@BAMs = split(/,/, $BAMs);
		die unless($#BAMs == $#sampleIDs);
	}
	else
	{
		@BAMs = map {'/gpfs1/well/gsk_hla/bam_output/'.$_.'.bam'} @sampleIDs;
	}
}
else
{
	@BAMs = glob('/gpfs1/well/gsk_hla/bam_output/*.bam');
	@sampleIDs = map {die unless($_ =~ /.+\/(.+?)\.bam/); $1} @BAMs;
}

if($action eq 'qsub')
{
	for(my $bI = 0; $bI <= $#BAMs; $bI++)
	{
		my $BAM = $BAMs[$bI];
		my $sampleID = $sampleIDs[$bI];
		
		unless(-e '../tmp/hla_qsub/')
		{
			mkdir('../tmp/hla_qsub/') or die;
		}
		
		die unless(defined $sampleID);
		
		my $dir_output = '../tmp/hla/'.$sampleID;
		unless(-e $dir_output)
		{
			mkdir($dir_output) or die "Cannot mkdir $dir_output";
		}
		
		my $time_file = $dir_output . '/time_p.txt';
		
		my $qsub_filename = '../tmp/hla_qsub/'.$sampleID.'.bash';
		
		
		my $command_positive_filtering = qq(perl HLAtypeinference.pl --graph $graph --sampleIDs $sampleID --BAMs $BAM --actions p --threads $threads);
		if($referenceGenome)
		{
			$command_positive_filtering .= qq( --referenceGenome $referenceGenome);
		}

		if($HiSeq250bp)
		{
			$command_positive_filtering .= qq( --HiSeq250bp 1);
		}
		
		print $command_positive_filtering, "\n";
		if($threads == 1)
		{
		open(QSUB, '>', $qsub_filename) or die "Cannot open $qsub_filename";
print QSUB qq(#!/bin/bash
#\$ -P mcvean.prja -q short.qa    
#\$ -pe shmem 3
export PERL5LIB=/users/mcvean/dilthey/perl5/lib/perl5:\$PERL5LIB
cd /gpfs1/well/gsk_hla/MHC-PRG/src
/usr/bin/time -v $command_positive_filtering &> $time_file
);
		close(QSUB);	
			
		}
		else
		{
			die "Unknown threads: $threads";
		}
		
		my $qsub_cmd = qq(qsub $qsub_filename);
		
		system($qsub_cmd);
	}
}
elsif($action eq 'copy')
{
	die "Can't access target directory $target_directory_for_copying" unless(-e $target_directory_for_copying);
	for(my $bI = 0; $bI <= $#BAMs; $bI++)
	{
		my $sampleID = $sampleIDs[$bI];
		# next if(($sampleID =~ /E2/) or ($sampleID =~ /E4/) or ($sampleID =~ /F4/)); # todo reactivate
		my $directory = '../tmp/hla/'.$sampleID;
		if(-e $directory)
		{	
			my $target = $target_directory_for_copying . '/I' . $iteration . '_' . $sampleID;
			if(-e $target)
			{
				warn "$target existing, skip";
				next;
			}
			my $cmd_copy = qq(cp -R $directory $target);
			
			my @files_input = glob($directory.'/*');
			@files_input = grep {$_ =~ /\.p_[12]$/} @files_input;
			die "Can't find two positive read files ".Dumper(\@files_input) unless(scalar(@files_input) == 2);
			
			my $md5_input_1 = `md5sum $files_input[0]`;			
			die "Can't parse wc1 output: $md5_input_1" unless($md5_input_1 =~ /^(\w+)\s+/);
			$md5_input_1 = $1;
			
			my $md5_input_2 = `md5sum $files_input[1]`;
			die "Can't parse wc2 output: $md5_input_2" unless($md5_input_2 =~ /^(\w+)\s+/);
			$md5_input_2 = $1;
			
			print "Executing command $cmd_copy\n";
			system($cmd_copy);
  
			my @files_output = glob($target.'/*');
			@files_output = grep {$_ =~ /\.p_[12]$/} @files_output;
			die "Can't find two positive read files ".Dumper(\@files_output) unless(scalar(@files_output) == 2);
			
			my $md5_output_1 = `md5sum $files_output[0]`;			
			die "Can't parse wc1 output: $md5_output_1" unless($md5_output_1 =~ /^(\w+)\s+/);
			$md5_output_1 = $1;
			
			my $md5_output_2 = `md5sum $files_output[1]`;
			die "Can't parse wc2 output: $md5_output_2" unless($md5_output_2 =~ /^(\w+)\s+/);
			$md5_output_2 = $1;
			
			unless($md5_input_1 eq $md5_output_1)
			{
				die "MD5 for $files_output[0] don't agree - $md5_input_1 / $md5_output_1";
			}
			
			unless($md5_input_2 eq $md5_output_2)
			{
				die "MD5 for $files_output[1] don't agree - $md5_input_2 / $md5_output_2";
			}
						

		}
	}
}
else
{
	die "Unknown --action!";
}

