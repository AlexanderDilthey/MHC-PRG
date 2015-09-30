#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/all/;
use File::Path;

my @coverage_targets = (40, 30, 20);
my $iterations = 3;
my %samples = (
	'I1_AA02O9Q_Z2' => 53.4,
#	'I3_NA12878' => 60.7,
);

open(CMD, '>', '_downsample_inference_commands.txt') or die;
foreach my $sample (keys %samples)
{
	my $current_coverage = $samples{$sample};
	die "Coverage problem" unless(all {$_ < $current_coverage} @coverage_targets);
	
	my $reads_file_1 = '../tmp/hla/' . $sample . '/reads.p.n_1';
	my $reads_file_2 = '../tmp/hla/' . $sample . '/reads.p.n_2';
	die unless(-e $reads_file_1);
	die unless(-e $reads_file_2);
	
	my @sampleIDs;
	my @sampleIDs_perTarget;
	
	foreach my $target (@coverage_targets)
	{
		my @sampleIDs_thisTarget;
		for(my $iteration = 1; $iteration <= $iterations; $iteration++)
		{
			my $downsampling_factor = $target / $current_coverage;
			print "Sample $sample target ${target}x iteration $iteration downsampling factor ", sprintf("%.2f", $downsampling_factor), "\n";
			
			my $new_sampleID = 'downsample_' . $sample . '_DSC' . $target . '_' . $iteration;
			my $output_dir = '../tmp/hla/' . $new_sampleID;
			
			if(-e $output_dir)
			{
				rmtree($output_dir);
			}
			mkdir($output_dir) or die "Cannot mkdir $output_dir";
			
			my $output_file_1 = $output_dir . '/reads.p.n_1';
			my $output_file_2 = $output_dir . '/reads.p.n_2';
			
			my $fh_in_1;
			my $fh_in_2;
			open($fh_in_1, '<', $reads_file_1) or die "Cannot open $reads_file_1";
			open($fh_in_2, '<', $reads_file_2) or die "Cannot open $reads_file_2";
			
			print "\t$reads_file_1 => $output_file_1\n";
			print "\t$reads_file_2 => $output_file_2\n";
			
			my $reads_read = 0;
			my $reads_printed = 0;
			open(READSOUT1, '>', $output_file_1) or die "Cannot open $output_file_1";			
			open(READSOUT2, '>', $output_file_2) or die "Cannot open $output_file_2";
			while(! eof($fh_in_1))
			{
				my @r = get_read_lines($fh_in_1, $fh_in_2);
				$reads_read++;
				if(rand() <= $downsampling_factor)
				{
					print READSOUT1 join('', @{$r[0]});
					print READSOUT2 join('', @{$r[1]});
					$reads_printed++;
					
				}
			}
			close(READSOUT2);			
			close(READSOUT1);	
			my $achieved_downsampling_factor = ($reads_read != 0) ? ($reads_printed / $reads_read) : -1;
			
			print "\t$reads_read reads read, $reads_printed reads printed, achieved factor: ", sprintf("%.2f", $achieved_downsampling_factor), "\n";
			close($fh_in_1);
			close($fh_in_2);
			
			push(@sampleIDs, $new_sampleID);
			push(@sampleIDs_thisTarget, $new_sampleID);
		}
		print "\n";
		push(@sampleIDs_perTarget, \@sampleIDs_thisTarget);
	}
	
	print CMD qq(./HLAtypeinference.pl --actions ai --sampleIDs ).join(',', @sampleIDs)."\n";
	foreach my $sampleList (@sampleIDs_perTarget)
	{
		print CMD qq(./HLAtypeinference.pl --actions v ---trueHLA /file/validation --sampleIDs ).join(',', @$sampleList)."\n";
	}
	print CMD "\n";
}

close(CMD);

sub get_read_lines
{
	my $fh_1 = shift;
	my $fh_2 = shift;
	
	my @l1, my @l2;
	for(my $i = 1; $i <= 4; $i++)
	{
		die if eof($fh_1);
		die if eof($fh_2);
		
		my $l1 = <$fh_1>;
		my $l2 = <$fh_2>;
		
		push(@l1, $l1);
		push(@l2, $l2);
		
	}
	
	die unless(substr($l1[0], 0, 1) eq '@');
	die unless(substr($l2[0], 0, 1) eq '@');
	die unless(substr($l1[2], 0, 1) eq '+');
	die unless(substr($l2[2], 0, 1) eq '+');
	
	return (\@l1, \@l2);
}