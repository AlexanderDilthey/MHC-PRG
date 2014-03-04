#!/usr/bin/perl

# This script identifies positions in the alignment which carry insertions relative to PGF; which samples contain calls to these
# these insertions; and where the re-mapping identified additional variation on top of these insertions.

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use List::Util qw/min max/;
use Data::Dumper;
use Getopt::Long;
use Sys::Hostname;

$| = 1;  

my $graph = 'varigraph_new';
my $begin_pos;
my $end_pos;
my $original_alignment = '/gpfs1/well/gsk_hla/shared/mhc_ref_8_haplotypes/alignment/all_aligned.fasta';
my $haplotype = 'pgf';

my $kMer_size = 55; # This is curious and only required to reconstruct the original alignment layout - should stay at k = 55, I believe.
my $kMer_size_inference = 31; # This influences filenames!
my $pgf_start = 28702185;
#my $printing_threshold_length = 50;
my $nearby_genes_distance = 10000;
my $temp_output_pgf_deletions = 'temp/_identified_PGF_deletions';

my $remapping_flank_length = 100;

GetOptions (
	'graph:s' => \$graph,
	'original_alignment:s' => \$original_alignment,
);

die "Alignment file $original_alignment not existing" unless (-e $original_alignment);

# Read graph loci and translate graph positions into PGF

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);
print "Extracting from graph: $graph\n";

my $genes_href = &get_chr6_genes();
my $graph_loci_data = get_graph_loci($graph_dir);
my @loci_in_graph = @{$graph_loci_data->[0]};
my @graph_locus_alleles = @{$graph_loci_data->[2]};

my %pgf_to_graphAlignment;
my $last_used_pgf = -1;
my $first_position_PGF = -1;
my $last_position_PGF = -1;
my %PGF_within_HLA;
foreach my $segment_file (keys %{$graph_loci_data->[1]})
{
	if($segment_file =~ /HLA\w+\.txt/)
	{
		foreach my $f (@{$graph_loci_data->[1]{$segment_file}})
		{
			my @parts = split(/_/, $f);
			die unless($#parts == 2);
			$PGF_within_HLA{$parts[2]} = 1;
		}
	}
}
  
for(my $i = 0; $i <= $#loci_in_graph; $i++)
{
	my @parts = split(/_/, $loci_in_graph[$i]);
	my $pgfPos;
	
	if($#parts == 2)
	{
		$pgfPos = $parts[2];	
		if($i == 0)
		{
			$first_position_PGF	= $pgfPos;
		}
		if($i == $#loci_in_graph)
		{
			$last_position_PGF	= $pgfPos;		
		}
	}
	else
	{
		die $loci_in_graph[$i];
		die unless($last_used_pgf > 0);
		$pgfPos = $last_used_pgf + 1;
	}
	
	if(not exists $pgf_to_graphAlignment{$pgfPos})
	{
		$pgf_to_graphAlignment{$pgfPos} = [];
	}
	push(@{$pgf_to_graphAlignment{$pgfPos}}, $i);
	$last_used_pgf = $pgfPos;	
}


die unless($first_position_PGF != -1);
die unless($last_position_PGF != -1);

# Read alignment and identify stretches of interest (where PGF has sequence, but other haplotypes have gaps

print "Original alignment: $original_alignment\n";

my $original_alignment_complete = read_alignment($original_alignment);
my $original_alignment = $original_alignment_complete->{alignment};
my $original_alignment_positions = $original_alignment_complete->{positions};
my $pgf_haplotype = $original_alignment->{'pgf'};
die unless($pgf_haplotype);

my %deletion_start_perHaplotype;
my %identified_deletions;
my $addGap = sub {
	my $p1 = shift;
	my $p2 = shift;
	die unless($p1 <= $p2);
	my $gap_str = join('--', $p1, $p2);
	$identified_deletions{$gap_str} = [$p1, $p2];
};

my $firstGapOver = 0;
for(my $posInAlignment = 0; $posInAlignment < length($pgf_haplotype); $posInAlignment++)
{
	my $pgf_character = substr($pgf_haplotype, $posInAlignment, 1);
	
	foreach my $haplotypeID (keys %{$original_alignment})
	{
		next if ($haplotypeID eq 'pgf');
		my $haplotype_character = substr($original_alignment->{$haplotypeID}, $posInAlignment, 1);
		if($pgf_character ne '_')
		{
			$firstGapOver = 1;
		}
		next unless($firstGapOver);
		
		if($pgf_character eq '_')
		{
			# pgf is gap, thus we will end all existing gaps and not start another one
			# from this position
			my $all_gap_end_position = $posInAlignment - 1;
			die unless($all_gap_end_position >= 0);
			foreach my $existingGapHaploID (keys %deletion_start_perHaplotype)
			{
				my $thisGapStartPos = $deletion_start_perHaplotype{$existingGapHaploID};
				die unless($thisGapStartPos <= $all_gap_end_position);
				$addGap->($thisGapStartPos, $all_gap_end_position);				
			}
			%deletion_start_perHaplotype = ();
		}
		else
		{
			if($haplotype_character eq '_')
			{
				if($deletion_start_perHaplotype{$haplotypeID})
				{
					# continue deletion
				}
				else
				{
					$deletion_start_perHaplotype{$haplotypeID} = $posInAlignment;
				}
			}
			else
			{
				if($deletion_start_perHaplotype{$haplotypeID})
				{
					# end deletion
					my $all_gap_end_position = $posInAlignment - 1;
					die unless($deletion_start_perHaplotype{$haplotypeID} <= $all_gap_end_position);					
					$addGap->($deletion_start_perHaplotype{$haplotypeID}, $all_gap_end_position);
					delete $deletion_start_perHaplotype{$haplotypeID};
				}
				else
				{
					# continue non-deletion
				}				
			}		
		}
	}
}

my @long_deletions_keys = grep {my $l = ($identified_deletions{$_}[1] - $identified_deletions{$_}[0]); (($l > 150))} keys %identified_deletions;
print "Identified ", scalar(keys %identified_deletions), " deletions in total and ", scalar(@long_deletions_keys), " after filtering for length\n";

# Create output directory for deletion data

unless(-e $temp_output_pgf_deletions)
{
	mkdir($temp_output_pgf_deletions) or die;
}
foreach my $existingFile (glob($temp_output_pgf_deletions.'/*'))
{
	unlink($existingFile) or die;
}

my $temp_output_summary = $temp_output_pgf_deletions.'/summary.txt';
open(SUMMARY, '>', $temp_output_summary) or die "Cannot open $temp_output_summary";
print SUMMARY join("\t", qw/i alignment_start alignment_stop pgf_start pgf_stop genomic_position inSamples inWhichSamples /, 'nearbyGenes'.$nearby_genes_distance), "\n";

# go through all Viterbi haplotypes and find out where there was a call on our deletions of interest

my @haplotype_files = glob('../tmp/*.viterbiHaplotypes');
die "No haplotype files to analyze?" unless(@haplotype_files);

my %sample_haplotypes;
my %sample_haplotypes_positions;
my %sample_haplotype_IGV_filenames;
my %sample_haplotype_IGV_filehandles;
foreach my $file (@haplotype_files)
{
	next unless($file =~ /Z1/);

	my $pattern = $graph.'_(\S+?)_'.$kMer_size_inference;
	next unless($file =~ /$pattern/);
	
	my $sampleID = $1;
	my $haplotypes_aref = read_haplotypes($file, \@loci_in_graph);
	
	$sample_haplotypes{$sampleID.'_H1'} = $haplotypes_aref->[0][0];
	$sample_haplotypes{$sampleID.'_H2'} = $haplotypes_aref->[0][1];	
	
	$sample_haplotypes_positions{$sampleID.'_H1'} = [];
	my $counter_nonGap = 0;
	for(split(//, $sample_haplotypes{$sampleID.'_H1'}))
	{
		if($_ ne '_')
		{
			$counter_nonGap++;
		}
		push(@{$sample_haplotypes_positions{$sampleID.'_H1'}}, $counter_nonGap); 
	}
	
	$sample_haplotypes_positions{$sampleID.'_H2'} = [];
	$counter_nonGap = 0;
	for(split(//, $sample_haplotypes{$sampleID.'_H2'}))
	{
		if($_ ne '_')
		{
			$counter_nonGap++;
		}
		push(@{$sample_haplotypes_positions{$sampleID.'_H2'}}, $counter_nonGap); 
	}
	
	$sample_haplotype_IGV_filenames{$sampleID.'_H1'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions/'.$sampleID.'_PnP2Viterbi_H1.txt';
	$sample_haplotype_IGV_filenames{$sampleID.'_H2'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions/'.$sampleID.'_PnP2Viterbi_H2.txt';
	$sample_haplotype_IGV_filenames{$sampleID.'_PGF'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions/'.$sampleID.'_PnP2Viterbi_PGF.txt';
	
	my $fh1, my $fh2, my $fh3;
	open($fh1, '>', $sample_haplotype_IGV_filenames{$sampleID.'_H1'}) or die "Cannot open ".$sample_haplotype_IGV_filenames{$sampleID.'_H1'};
	open($fh2, '>', $sample_haplotype_IGV_filenames{$sampleID.'_H2'}) or die "Cannot open ".$sample_haplotype_IGV_filenames{$sampleID.'_H2'};
	open($fh3, '>', $sample_haplotype_IGV_filenames{$sampleID.'_PGF'}) or die "Cannot open ".$sample_haplotype_IGV_filenames{$sampleID.'_PGF'};
	
	
	$sample_haplotype_IGV_filehandles{$sampleID.'_H1'} = $fh1;
	$sample_haplotype_IGV_filehandles{$sampleID.'_H2'} = $fh2;
	$sample_haplotype_IGV_filehandles{$sampleID.'_PGF'} = $fh3;
	
	print { $fh1 } qq(new
genome ${sampleID}_1
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_1.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions
);

	print { $fh2 } qq(new
genome ${sampleID}_2
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_2.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions
);

	print { $fh3 } qq(new
genome hg19
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\samtools_merged_xMHC.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions
);
}

my $deletion_number = 0;
foreach my $deletionKey (@long_deletions_keys)
{
	$deletion_number++;
	
	my $deletionStart_Alignment = $identified_deletions{$deletionKey}[0];
	my $deletionStop_Alignment = $identified_deletions{$deletionKey}[1];
	
	die unless($original_alignment_positions->{pgf});
	my $deletionStart_PGF = $original_alignment_positions->{pgf}[$deletionStart_Alignment] - 1;
	my $deletionStop_PGF = $original_alignment_positions->{pgf}[$deletionStop_Alignment] - 1;
	my $deletionStart_Genomic = $pgf_start + $deletionStart_PGF;

	next if($PGF_within_HLA{$deletionStart_PGF} or $PGF_within_HLA{$deletionStop_PGF});
	
	die unless($pgf_to_graphAlignment{$deletionStart_PGF} and $pgf_to_graphAlignment{$deletionStop_PGF});
	my $deletionStart_Graph = $pgf_to_graphAlignment{$deletionStart_PGF}[0];
	my $deletionStop_Graph = $pgf_to_graphAlignment{$deletionStop_PGF}[0];
	
	die unless($deletionStop_Graph >= $deletionStart_Graph);
	
	
	my %deletion_in_haploIDs;
	foreach my $haplotypeID (keys %sample_haplotypes)
	{
		my $haplotype = $sample_haplotypes{$haplotypeID};
		
		my $extraction_area = substr($haplotype, $deletionStart_Graph, $deletionStop_Graph - $deletionStart_Graph + 1);
		if($extraction_area =~ /^_+$/)
		{
			my $inHaplotype_deletion_start = $sample_haplotypes_positions{$haplotypeID}[$deletionStart_Graph] + $remapping_flank_length + 1;
			my $inHaplotype_deletion_stop = $sample_haplotypes_positions{$haplotypeID}[$deletionStop_Graph] + $remapping_flank_length + 1;
			
			my $key_with_position = $haplotypeID.' ['.$inHaplotype_deletion_start.'-'.$inHaplotype_deletion_start.']';
			$deletion_in_haploIDs{$key_with_position} = 1;
		
			my $display_rightCoordinate_Genomic = $deletionStart_Genomic - 400;		
			my $display_leftCoordinate_Genomic = $deletionStart_Genomic + 400;
			my $filename_for_IGV_snapshot_PGF = "${deletion_number}_${deletionStart_Genomic}_PGF.png";
			my $pgf_haploID = $haplotypeID;
			$pgf_haploID =~ s/_H[12]/_PGF/;
			print { $sample_haplotype_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${display_rightCoordinate_Genomic}-${display_leftCoordinate_Genomic}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
			
			if($haplotypeID =~ /H1/)
			{
				my $filename_for_IGV_snapshot_H1 = "${deletion_number}_${deletionStart_Genomic}_H1.png";
				print { $sample_haplotype_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${inHaplotype_deletion_start}\nsnapshot ${filename_for_IGV_snapshot_H1}\n);										
			}
			else
			{
				die unless($haplotypeID =~ /H2/);
				my $filename_for_IGV_snapshot_H2 = "${deletion_number}_${deletionStart_Genomic}_H2.png";
				print { $sample_haplotype_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${inHaplotype_deletion_start}\nsnapshot ${filename_for_IGV_snapshot_H2}\n);								
			}
		}
	}
	
	if(scalar(keys %deletion_in_haploIDs) > 0)
	{
		# open deletion-specific file
		open(F, '>', $temp_output_pgf_deletions.'/deletion'.$deletion_number.'.txt') or die;
		print F "Deletion $deletion_number\n";
		print F "Alignment: $deletionStart_Alignment - $deletionStop_Alignment \n";
		print F "PGF: $deletionStart_PGF - $deletionStop_PGF (genomic start $deletionStart_Genomic -- IGV plus 100 or so!!) \n";
		print F "Graph: $deletionStart_Graph - $deletionStop_Graph\n\n\n";
		
		foreach my $haplotype (sort keys %{$original_alignment})
		{
			my $haplotype_string = $original_alignment->{$haplotype};

			my $insertion_sequence = substr($haplotype_string, $deletionStart_Alignment, $deletionStop_Alignment - $deletionStart_Alignment + 1);
			my $pre_extraction = substr($haplotype_string, $deletionStart_Alignment-15, 15);
			my $post_extraction = substr($haplotype_string, $deletionStop_Alignment+1, 15);
					
			print F sprintf("%20s", $haplotype), "\t", $pre_extraction, '  ', $insertion_sequence, '  ', $post_extraction, "\n";
		}
		
		print F "\n\n";
	
		foreach my $haplotypeID (keys %sample_haplotypes)
		{
			my $haplotype = $sample_haplotypes{$haplotypeID};
			
			my $inHaplotype_deletion_start = $sample_haplotypes_positions{$haplotypeID}[$deletionStart_Graph] + $remapping_flank_length + 1;
			my $inHaplotype_deletion_stop = $sample_haplotypes_positions{$haplotypeID}[$deletionStop_Graph] + $remapping_flank_length + 1;
				
			my $extraction_area = substr($haplotype, $deletionStart_Graph, $deletionStop_Graph - $deletionStart_Graph + 1);
			my $pre_extraction = substr($haplotype, $deletionStart_Graph-15, 15);
			my $post_extraction = substr($haplotype, $deletionStop_Graph+1, 15);

			print F sprintf("%30s", $haplotypeID), "\t", $inHaplotype_deletion_start, "  ",
$pre_extraction, "  ", $extraction_area, "  ", $post_extraction, "  ", $inHaplotype_deletion_stop, "\n";
		}
	
		close(F);
	}
	
	my @nearby_genes;
	foreach my $gene (keys %$genes_href)
	{
		my $start = $genes_href->{$gene}[0];
		my $stop = $genes_href->{$gene}[0];
		my $start_distance = abs($deletionStart_PGF - $start);
		my $stop_distance = abs($deletionStart_PGF - $stop);
		my $start_II_distance = abs($deletionStop_PGF - $start);
		my $stop_II_distance = abs($deletionStop_PGF - $stop);
		
		my $min_I_distance = ($start_distance < $stop_distance) ? $start_distance : $stop_distance;
		my $min_II_distance = ($start_II_distance < $stop_II_distance) ? $start_II_distance : $stop_II_distance;
		
		my $min_distance = ($min_I_distance < $min_II_distance) ? $min_I_distance : $min_II_distance;
		
		die unless($min_distance >= 0);
		if($min_distance < $nearby_genes_distance)
		{
			push(@nearby_genes, $gene.'-'.$min_distance);
		}
	}	
		
	print SUMMARY join("\t",
		$deletion_number,
		$deletionStart_Alignment,
		$deletionStop_Alignment,
		$deletionStart_PGF,
		$deletionStop_PGF,
		$deletionStart_Genomic,
		scalar(keys %deletion_in_haploIDs),
		join(';', keys %deletion_in_haploIDs),
		join(';', @nearby_genes),
	), "\n";
	
}

close(SUMMARY);

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
	close(ALIGNMENT);
	
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

sub get_graph_loci
{
	my $graph_dir = shift;
	
	my %locus_origin;
	
	# die $graph_dir;
	
	my @loci_in_graph;
	my @locus_alleles;
	
	my $graph_segments_file = $graph_dir.'/segments.txt';
	die "File not there: $graph_segments_file" unless (-e $graph_segments_file);
	open(SEGMENTS, '<', $graph_segments_file) or die;
	while(<SEGMENTS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
		next unless($line);
		
		my $segment_file = $graph_dir.'/'.$line;

		open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
		my $firstLine = <SEGMENT>;
		chomp($firstLine);
		$firstLine =~ s/\n//g;
		$firstLine =~ s/\r//g;			
		my @line_fields = split(/ /, $firstLine);
		shift(@line_fields); # kick out individual ID
		
		push(@{$locus_origin{$segment_file}}, @line_fields);
		push(@loci_in_graph, @line_fields);
		
		my $beginning_of_locus_alleles = $#locus_alleles + 1;
		push(@locus_alleles, map {{}} @line_fields);
		
		while(<SEGMENT>)
		{
			my $line = <SEGMENT>;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;			
			my @this_line_fields = split(/ /, $line);
			shift(@this_line_fields); # kick out individual ID		
			for(my $j = 0; $j <= $#this_line_fields; $j++)
			{
				$locus_alleles[$beginning_of_locus_alleles + $j]{$this_line_fields[$j]}++;
			}
		}
		close(SEGMENT);
	}
	close(SEGMENTS);
	
	return [\@loci_in_graph, \%locus_origin, \@locus_alleles];
}

sub read_haplotypes
{
	my $haplotypes_file = shift;
	my $loci_in_graph_aref = shift;
	
	open(HAPLOTYPES, '<', $haplotypes_file) or die "Cannot open $haplotypes_file";
	my $haplotype_1, my $haplotype_2;
	my @haplotype_1_alignment_fields;
	my @haplotype_2_alignment_fields;
	while(<HAPLOTYPES>)
	{
		my $line = $_;
		chomp($line);
		last if ($. > 2);
		my @fields = split(/ /, $line);
		my $haplotype = join('', @fields[2 .. ($#fields-2)]);
		if($. == 1)
		{
			$haplotype_1 = $haplotype;
			@haplotype_1_alignment_fields = @fields[2 .. ($#fields-2)];
		}
		elsif($. == 2)
		{
			$haplotype_2 = $haplotype;	
			@haplotype_2_alignment_fields = @fields[2 .. ($#fields-2)];
		}
	}
	close(HAPLOTYPES);

	die unless($haplotype_1 and $haplotype_2);
	die Dumper('$#haplotype_1_alignment_fields == $#{$loci_in_graph_aref}', $haplotypes_file, $#haplotype_1_alignment_fields, $#{$loci_in_graph_aref}) unless($#haplotype_1_alignment_fields == $#{$loci_in_graph_aref});
	die Dumper('$#haplotype_2_alignment_fields == $#{$loci_in_graph_aref}', $haplotypes_file, $#haplotype_2_alignment_fields, $#{$loci_in_graph_aref})  unless($#haplotype_2_alignment_fields == $#{$loci_in_graph_aref});
	
	return [[$haplotype_1, $haplotype_2], [\@haplotype_1_alignment_fields, \@haplotype_2_alignment_fields]];
}

sub get_chr6_genes
{
	my %genes;
	my $input = '/Net/birch/data/dilthey/MHC-PRG/src/temp/UCSC_chr6_RefSeq_genes.txt';
	open(F, '<', $input) or die "Cannot open $input";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next if(substr($line, 0, 1 ) eq '#');
		my @fields = split(/\t/, $line);
		my $gene_name = $fields[5];
		my $gene_start = $fields[3];
		my $gene_stop = $fields[4];
		die unless($fields[1] eq 'chr6');
		$genes{$gene_name} = [$gene_start, $gene_stop];
	}
	close(F);
	
	return \%genes;
}

