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

my $referenceGenome = qq(/gpfs1/well/chimp/oa/ref/hs37d5.fasta);
#$referenceGenome = '';

GetOptions ('graph:s' => \$graph,
 'action:s' => \$action, 
 'iteration:s' => \$iteration,  
 'sampleIDs:s' => \$sampleIDs,
);         

my $target_directory_for_copying = qq(/Net/birch/data/dilthey/MHC-PRG/tmp/hla);

my @BAMs;
my @sampleIDs;
if($sampleIDs)
{
	@sampleIDs = split(/,/, $sampleIDs);
	@BAMs = map {'/gpfs1/well/gsk_hla/bam_output/'.$_.'.bam'} @sampleIDs;
}
else
{
	@BAMs = glob('/gpfs1/well/gsk_hla/bam_output/*.bam');
	@sampleIDs = map {die unless($_ =~ /.+\/(.+?)\.bam/); $1} @BAMs;
}

#die Dumper(\@BAMs);

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
		
		my $qsub_filename = '../tmp/hla_qsub/'.$sampleID.'.bash';
		
		my $command_negative_filtering = qq(perl HLAtypeinference.pl --graph $graph --sampleIDs $sampleID --BAMs $BAM --actions p);
		if($referenceGenome)
		{
			$command_negative_filtering .= qq( --referenceGenome $referenceGenome);
		}
		
		print $command_negative_filtering, "\n";
		open(QSUB, '>', $qsub_filename) or die "Cannot open $qsub_filename";
print QSUB qq(#!/bin/bash
#\$ -P mcvean.prjb -q long.qb
#\$ -pe shmem 2
export PERL5LIB=/users/dilthey/perl/lib/x86_64-linux-thread-multi/:/opt/perl/5.8.8/bioperl/1.6.901/Bio-SamTools-1.38/lib/perl5/x86_64-linux-thread-multi:/opt/perl/5.8.8/bioperl/1.6.901/bioperl-live:/users/dilthey/perl5/lib/perl5/x86_64-linux-thread-multi/:/users/dilthey/perl5/lib/perl5:\$PERL5LIB
cd /gpfs1/well/gsk_hla/MHC-PRG/src
$command_negative_filtering
);
		close(QSUB);	
		
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
			
			print "Executing command $cmd_copy\n";
			system($cmd_copy);

			my @files = glob($target.'/*');
			@files = grep {$_ =~ /\.p_[12]$/} @files;
			die "Can't find two positive read files ".Dumper(\@files) unless(scalar(@files) == 2);
			
			my $wc1_output = `wc $files[0]`;			
			die "Can't parse wc1 output: $wc1_output" unless($wc1_output =~ /^\s*(\d+)\s+(\d+)\s+(\d+)\s+/);
			my $lines1 = $1;
			
			my $wc2_output = `wc $files[1]`;
			die "Can't parse wc2 output: $wc2_output" unless($wc2_output =~ /^\s*(\d+)\s+(\d+)\s+(\d+)\s+/);
			my $lines2 = $1;
			
			unless($lines1 == $lines2)
			{
				die "Lines for $files[0] and $files[1] don't agree - $lines1 / $lines2 \n$wc1_output\n$wc2_output";
			}
						

		}
	}
}
else
{
	die "Unknown --action!";
}

