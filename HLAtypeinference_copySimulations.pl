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
	
	# true HLA
	my $trueHLA = $dir . '/trueHLA.txt';
	my $trueHLA_out = '../tmp/hla/trueHLA_'.$iteration.'_'.$baseName_outer;
	
	open(TRUE, '<', $trueHLA) or die "Cannot open $trueHLA";
	open(TRUEOUT, '>', $trueHLA_out) or die "Cannot open $trueHLA_out";
	
	my $headerLine = <TRUE>;
	chomp($headerLine);
	my @header_fields = split(/ /, $headerLine);
	for(my $i = 1; $i <= $#header_fields; $i++)
	{
		$header_fields[$i] = 'HLA'.$header_fields[$i];
	}
	print TRUEOUT join("\t", @header_fields), "\n";
	while(<TRUE>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/ /, $line);
		
		die unless($fields[0] =~ /S_(\d+)/);
		my $sampleN = $1;
			
		$fields[0] = 'simulations_'.$baseName_outer.'_sample'.$sampleN;
		
		for(my $i = 1; $i <= $#fields; $i++)
		{
			my $v = $fields[$i];
			my @p = split(/\//, $v);
			die unless($#p == 1);
			for(my $j = 0; $j <= $#p; $j++)
			{
				my @p2 = split(/\:/, $p[$j]);
				if($#p2 >= 1)
				{
					$p[$j] = join(':', @p2[0, 1]);
				}
			}
			$v = join('/', @p);
			$fields[$i] = $v;
		}
		
		print TRUEOUT join("\t", @fields), "\n";
	}
	close(TRUE);
	close(TRUEOUT);
	
	print "Validation HLA file:\n$trueHLA_out\n\n";
	
	# true haplotypes
	my $trueHaplotypes = $dir . '/trueHaplotypes.txt';
	my $trueHaplotypes_out = '../tmp/hla/trueHaplotypes'.$iteration.'_'.$baseName_outer;
	
	open(TRUE, '<', $trueHaplotypes) or die "Cannot open $trueHaplotypes";
	open(TRUEOUT, '>', $trueHaplotypes_out) or die "Cannot open $trueHaplotypes_out";
	
	$headerLine = <TRUE>;
	chomp($headerLine);
	@header_fields = split(/ /, $headerLine);
	for(my $i = 1; $i <= $#header_fields; $i++)
	{
		$header_fields[$i] = 'HLA'.$header_fields[$i];
	}
	print TRUEOUT join("\t", @header_fields), "\n";
	while(<TRUE>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/ /, $line);
		
		die unless($fields[0] =~ /S_(\d+)/);
		my $sampleN = $1;
			
		$fields[0] = 'simulations_'.$baseName_outer.'_sample'.$sampleN;
		
		for(my $i = 1; $i <= $#fields; $i++)
		{
			my $v = $fields[$i];
			my @p = split(/\//, $v);
			die unless($#p == 1);
			for(my $j = 0; $j <= $#p; $j++)
			{
				my @p2 = split(/\:/, $p[$j]);
				if($#p2 >= 1)
				{
					$p[$j] = join('', @p2[0, 1]);
				}
			}
			$v = join('/', @p);
			$fields[$i] = $v;
		}
		
		print TRUEOUT join("\t", @fields), "\n";
	}
	close(TRUE);
	close(TRUEOUT);
	
	print "Validation haplotypes file:\n$trueHaplotypes_out\n\n";	
	
	# perturbed positions
	
	my $trueHaplotypesPerturbed = $dir . '/trueHaplotypes.txt.perturbed';
	my $trueHaplotypesPerturbed_out = '../tmp/hla/trueHaplotypes'.$iteration.'_'.$baseName_outer.".perturbed";
	
	open(TRUE, '<', $trueHaplotypesPerturbed) or die "Cannot open $trueHaplotypesPerturbed";
	open(TRUEOUT, '>', $trueHaplotypesPerturbed_out) or die "Cannot open $trueHaplotypesPerturbed_out";
	
	$headerLine = <TRUE>;
	chomp($headerLine);
	@header_fields = split(/ /, $headerLine);
	for(my $i = 1; $i <= $#header_fields; $i++)
	{
		$header_fields[$i] = 'HLA'.$header_fields[$i];
	}
	print TRUEOUT join("\t", @header_fields), "\n";
	while(<TRUE>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/ /, $line);
		
		die unless($fields[0] =~ /S_(\d+)/);
		my $sampleN = $1;
			
		$fields[0] = 'simulations_'.$baseName_outer.'_sample'.$sampleN;
		
		for(my $i = 1; $i <= $#fields; $i++)
		{
			my $v = $fields[$i];
			my @p = split(/\//, $v, -1);
			die "Problem with file $trueHaplotypesPerturbed: missing slash in line $. for field $header_fields[$i] - field value '$v' " unless($#p == 1);
			for(my $j = 0; $j <= $#p; $j++)
			{
				my @p2 = split(/\:/, $p[$j]);
				if($#p2 >= 1)
				{
					$p[$j] = join('', @p2[0, 1]);
				}
			}
			$v = join('/', @p);
			$fields[$i] = $v;
		}
		
		print TRUEOUT join("\t", @fields), "\n";
	}
	close(TRUE);
	close(TRUEOUT);
	
	print "File with perturbed positions: $trueHaplotypesPerturbed_out\n";	
		
	# files
	
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