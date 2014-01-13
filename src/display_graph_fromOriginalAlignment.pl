#!/usr/bin/perl

# This script displays a range of positons from the original alignment, but in the haplotype
# coordinate system!

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use List::Util qw/min max/;
use Data::Dumper;
use Getopt::Long;
use Sys::Hostname;

my $graph = 'varigraph_new';
my $begin_pos;
my $end_pos;
my $original_alignment = '/gpfs1/well/gsk_hla/shared/mhc_ref_8_haplotypes/alignment/all_aligned.fasta';
my $haplotype = 'pgf';

my $kMer_size = 55; # This is curious and only required to reconstruct the original alignment layout - should stay at k = 55, I believe.

GetOptions (
	'graph:s' => \$graph,
	'original_alignment:s' => \$original_alignment,
	'haplotype:s' => \$haplotype,	
	'begin_pos:s' => \$begin_pos,
	'end_pos:s' => \$end_pos,
);

die "Please provide --begin_pos and --end_pos" unless((defined $begin_pos) and (defined $end_pos));
die unless($begin_pos <= $end_pos);
die "Alignment file $original_alignment not existing" unless (-e $original_alignment);

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);

print "Extracting from graph: $graph\n";
print "Original alignment: $original_alignment\n";

my $original_alignment_complete = read_alignment($original_alignment);
my $original_alignment = $original_alignment_complete->{alignment};
my $original_alignment_positions = $original_alignment_complete->{positions};

unless (exists $original_alignment->{$haplotype})
{
	die "Haplotype $haplotype not defined in original alignment\n";
}
unless (exists $original_alignment_positions->{$haplotype})
{
	die "Haplotype $haplotype not defined in original alignment positions\n";
}
unless(length($original_alignment->{$haplotype}) == scalar(@{$original_alignment_positions->{$haplotype}}))
{
	die "Alignment string and positions for haplotpye $haplotype not in agreement!";
}

# follow the selected haplotype and find the selected positions
my $pgf_start = 28702185;

my $begin_pos_inHaplotype = $begin_pos;
my $end_pos_inHaplotype = $end_pos;
if($haplotype eq 'pgf')
{
	$begin_pos_inHaplotype = $begin_pos_inHaplotype - $pgf_start;
	$end_pos_inHaplotype = $end_pos_inHaplotype - $pgf_start;
}
my $first_alignmentColumn_matchBegin = -1;
my $last_alignmentColumn_matchEnd = -1;
for(my $i = 0; $i < scalar(@{$original_alignment_positions->{$haplotype}}); $i++)
{
	my $position_in_haplotype = $original_alignment_positions->{$haplotype}[$i];
	if($position_in_haplotype < $begin_pos_inHaplotype)
	{
		$first_alignmentColumn_matchBegin = $i;
	}
	else
	{
		last;
	}
}
die if($first_alignmentColumn_matchBegin == -1);
$first_alignmentColumn_matchBegin++;

for(my $i = scalar(@{$original_alignment_positions->{$haplotype}}) - 1; $i >= 0 ; $i--)
{
	my $position_in_haplotype = $original_alignment_positions->{$haplotype}[$i];
	if($position_in_haplotype > $end_pos_inHaplotype)
	{
		$last_alignmentColumn_matchEnd = $i;
	}
	else
	{
		last;
	}
}
die if($last_alignmentColumn_matchEnd == -1);
$last_alignmentColumn_matchEnd--;

print "\n\nExtraction info:\n";
print "Haplotype name: $haplotype\n";
print "Extraction positions:\n";
print "\tBegin: $original_alignment_positions->{$haplotype}[$first_alignmentColumn_matchBegin] (achieved), $begin_pos_inHaplotype(specified), $begin_pos (original specified)\n";
 print "\tEnd: $original_alignment_positions->{$haplotype}[$last_alignmentColumn_matchEnd] (achieved), $end_pos_inHaplotype (specified), $end_pos (original specified)\n";
 
 my $first_alignmentColumn_matchBegin_PGF = $original_alignment_positions->{'pgf'}[$first_alignmentColumn_matchBegin];
 my $last_alignmentColumn_matchEnd_PGF = $original_alignment_positions->{'pgf'}[$last_alignmentColumn_matchEnd];

 
print "\n\nTransformation into PGF coordinates:\n";
print "\tMatched first alignment column: $first_alignmentColumn_matchBegin => PGF column $first_alignmentColumn_matchBegin_PGF\n";
print "\tMatched last alignment column: $last_alignmentColumn_matchEnd => PGF column $last_alignmentColumn_matchEnd_PGF\n";

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

	my @line_fields_PGF =(undef);
	
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
			my $last_used_PGFidx = -1;
			for(my $i = 1; $i <= $#line_fields; $i++)
			{
				my @parts = split(/_/, $line_fields[$i]);
				if($#parts == 2)
				{
					$line_fields_PGF[$i] = $parts[2];
					$last_used_PGFidx = $parts[2];
				}
				else
				{
					if($last_used_PGFidx != -1)
					{
						$line_fields_PGF[$i] = $last_used_PGFidx;
					}
					else
					{
						$line_fields_PGF[$i] = undef;					
					}
				}
			}
			
			
			my $min_index = min(@line_fields_PGF);
			my $max_index = max(@line_fields_PGF);
			
			unless((($first_alignmentColumn_matchBegin_PGF >= $min_index) and ($first_alignmentColumn_matchBegin_PGF <= $max_index)) or (($last_alignmentColumn_matchEnd_PGF >= $min_index) and ($last_alignmentColumn_matchEnd_PGF <= $max_index)))
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
				my $thisFeldTotalIndexPGF = $line_fields_PGF[$i];
				if(($thisFeldTotalIndexPGF >= $first_alignmentColumn_matchBegin_PGF) and ($thisFeldTotalIndexPGF <= $last_alignmentColumn_matchEnd_PGF))
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

