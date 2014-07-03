#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;

# input parameters

my $iteration = 1;
GetOptions ('iteration:s' => \$iteration,
);         


my @dirs = glob('../tmp/hla/simulations/*');

foreach my $dir (@dirs)
{
	die unless($dir =~ /.+\/(.+)/);
	my $baseName_outer = $1;
	
	my @subfiles = glob($dir.'/*_1.fastq');
	
	foreach my $subfile (@subfiles)
	{
		my $subfile_2 = $subfile;
		$subfile_2 =~ s/_1\.fastq/_2.fastq/;
		die unless(-e $subfile_2);
		
		die unless($subfile =~ /.+\/(.+)/);		
		my $basename = $1;
		
		die unless($basename =~ /S_(\d+)_/);
		my $simulationNumber = $1;
		#die Dumper($subfile_2, $basename, $simulationNumber);
		
		my $targetDir = '../tmp/hla/I'.$iteration.'_simulations_'.$baseName_outer.'_sample'.$simulationNumber;
		
		mkdir($targetDir);
		die "Cannot mkdir $targetDir" unless(-e $targetDir);

		my $targetfile_1 = $targetDir . '/reads.p_1';
		my $targetfile_2 = $targetDir . '/reads.p_2';
		
		print "Copying into $targetDir ...\n";
		copy($subfile, $targetfile_1) or die "Cannot copy file";
		copy($subfile_2, $targetfile_2) or die "Cannot copy file";

	}
}