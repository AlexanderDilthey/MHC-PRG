#!/usr/bin/perl

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;
use Sys::Hostname;

my $graph = 'varigraph';
my $begin_pos;
my $end_pos;
GetOptions (
	'graph:s' => \$graph,
	'begin_pos:s' => \$begin_pos,
	'end_pos:s' => \$end_pos,
);

die "Please provide --begin_pos and --end_pos" unless((defined $begin_pos) and (defined $end_pos));
die unless($begin_pos <= $end_pos);

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);
print "Extracting from graph: $graph\n";

my @loci_in_graph;
my $graph_segments_file = $graph_dir.'/segments.txt';
die "File not there: $graph_segments_file" unless (-e $graph_segments_file);

my @print_matrix;
my $runningLocusOffset = 0;
my $runningMatrixOffset = 0;
open(SEGMENTS, '<', $graph_segments_file) or die;
while(<SEGMENTS>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/\n//g;
	$line =~ s/\r//g;	
	next unless($line);
	
	my $segment_file = $graph_dir.'/'.$line;

	my $thisSegmentLoci;
	open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
	
	my $lineNumber = -1;
	my $taken_fields = 0;	
	LINE: while(<SEGMENT>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;			
		my @line_fields = split(/ /, $line);
		
		$lineNumber++;		
		$taken_fields = 0;			
		
		if($lineNumber == 0)
		{
			$thisSegmentLoci = scalar(@line_fields) - 1;
			my $min_index = $runningLocusOffset;
			my $max_index = $runningLocusOffset + scalar(@line_fields) - 1;
			unless((($begin_pos >= $min_index) and ($begin_pos <= $max_index)) or (($end_pos >= $min_index) and ($end_pos <= $max_index)))
			{
				last LINE;
			}
		}
		
		for(my $i = 0; $i <= $#line_fields; $i++)
		{
			if($i == 0)
			{
				my $identifier = $line_fields[0];
				$print_matrix[$runningMatrixOffset][$lineNumber] = $identifier;
			}
			else
			{
				my $thisFeldTotalIndex = $runningLocusOffset + $i - 1;
				if(($thisFeldTotalIndex >= $begin_pos) and ($thisFeldTotalIndex <= $end_pos))
				{
					$print_matrix[$runningMatrixOffset + $taken_fields + 1][$lineNumber] = $line_fields[$i];
					$taken_fields++;
				}
			}
		}
	}
	close(SEGMENT);
	
	die unless($thisSegmentLoci);
	$runningLocusOffset += $thisSegmentLoci;
	$runningMatrixOffset += $taken_fields + 1 if($taken_fields > 0);
	
}
close(SEGMENTS);

my $max_print_y = -1;
foreach my $l (@print_matrix)
{ 
	if($#{$l} > $max_print_y)
	{
		$max_print_y = $#{$l};
	}
}

for(my $i = 0; $i <= $max_print_y; $i++)
{
	my @print_fields;
	for(my $x = 0; $x <= $#print_matrix; $x++)
	{ 
		if(defined $print_matrix[$x][$i])
		{
			push(@print_fields, $print_matrix[$x][$i]);
		}
		else
		{
			push(@print_fields, undef);		
		}
	}
	
	print join("\t", @print_fields), "\n";
}
