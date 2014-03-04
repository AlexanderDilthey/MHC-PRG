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

my $remapping_flank_length = 100;


GetOptions (
	'graph:s' => \$graph,
	'original_alignment:s' => \$original_alignment,
);



die "Alignment file $original_alignment not existing" unless (-e $original_alignment);

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);

my $genes_href = &get_chr6_genes();

my $graph_loci_data = get_graph_loci($graph_dir);
my @loci_in_graph = @{$graph_loci_data->[0]};

my @graph_locus_alleles = @{$graph_loci_data->[2]};

# create index into graph PGF positions
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

# die Dumper($pgf_to_graphAlignment{1877}, $pgf_to_graphAlignment{1878}, $pgf_to_graphAlignment{1879});

print "Extracting from graph: $graph\n";
print "Original alignment: $original_alignment\n";

my $original_alignment_complete = read_alignment($original_alignment);
my $original_alignment = $original_alignment_complete->{alignment};
my $original_alignment_positions = $original_alignment_complete->{positions};

my $pgf_haplotype = $original_alignment->{'pgf'};
die unless($pgf_haplotype);
my $pgf_noGaps = $pgf_haplotype;
$pgf_noGaps =~ s/[^ACGTacgt]//g;
die if($pgf_noGaps =~ /\*/);
die if($pgf_noGaps =~ /\_/);

my @pgf_insertions;
my $in_insertion = 0;
my $insertion_begin;
my $last_scanned_position;

my $insertion_has_non_StarGap_haplotype = sub {
	my $insertion_spec = shift;
	foreach my $haplotypeID (sort keys %{$original_alignment})
	{
		my $insertion_sequence = substr($original_alignment->{$haplotypeID}, $insertion_spec->[0],  $insertion_spec->[1] -  $insertion_spec->[0] + 1);
		if($insertion_sequence =~ /[^\*_]/)
		{
			return 1;
		}
	}
	return 0;
};

for(my $i = 0; $i < length($pgf_haplotype); $i++)
{
	my $pgf_unaligned_position = $original_alignment_positions->{pgf}[$i];
	die unless(defined $pgf_unaligned_position);
	if($pgf_unaligned_position < $first_position_PGF)
	{
		next;
	}	
	if($pgf_unaligned_position > $last_position_PGF)
	{
		last;
	}	
	
	$last_scanned_position = $i;
	my $char_pgf = substr($pgf_haplotype, $i, 1);
	if($in_insertion)
	{
		if($char_pgf eq '_')
		{
			# continue with insertion
		}
		else
		{
			# insertion end
			my $insertion_last_position = $i - 1;
			my $ins_spec = [$insertion_begin, $insertion_last_position];
			if($insertion_has_non_StarGap_haplotype->($ins_spec))
			{
				push(@pgf_insertions, $ins_spec);
			}
			$in_insertion = 0;
		}	
	}
	else
	{
		if($char_pgf eq '_')
		{
			$in_insertion = 1;
			$insertion_begin = $i;
		}
		else
		{
			# continue normal!
		}
	}
}
if($in_insertion)
{
	my $insertion_last_position = $last_scanned_position;
	my $ins_spec = [$insertion_begin, $insertion_last_position];	
	if($insertion_has_non_StarGap_haplotype->($ins_spec))
	{
		push(@pgf_insertions, $ins_spec);
	}
}

print "Identified ", scalar(@pgf_insertions), " insertions relative to PGF";

my $temp_output_pgf_insertions = 'temp/_identified_PGF_insertions';
unless(-e $temp_output_pgf_insertions)
{
	mkdir($temp_output_pgf_insertions) or die;
}
foreach my $existingFile (glob($temp_output_pgf_insertions.'/*'))
{
	unlink($existingFile) or die;
}

my %extraction_targets;
my %long_insertions;

my @haplotype_files = glob('../tmp/*.viterbiHaplotypes');
die "No haplotype files to analyze?" unless(@haplotype_files);

my %sample_haplotypes;
my %sample_haplotype_positions;
my %sample_haplotype_insertion_IGV_filenames;
my %sample_haplotype_insertion_IGV_filehandles;
my %sample_haplotype_novelSeq_IGV_filenames;
my %sample_haplotype_novelSeq_IGV_filehandles;
foreach my $file (@haplotype_files)
{
	my $pattern = $graph.'_(\S+?)_'.$kMer_size_inference;
	next unless($file =~ /Z1/);
	
	next unless($file =~ /$pattern/);
		
	my $sampleID = $1;
	my $haplotypes_aref = read_haplotypes($file, \@loci_in_graph);
	$sample_haplotypes{$sampleID.'_H1'} = $haplotypes_aref->[0][0];
	$sample_haplotypes{$sampleID.'_H2'} = $haplotypes_aref->[0][1];
	
	$sample_haplotype_positions{$sampleID.'_H1'} = [];
	$sample_haplotype_positions{$sampleID.'_H2'} = [];
	
	my $non_gap_symbols_H1 = 0;
	for(my $i = 0; $i < length($sample_haplotypes{$sampleID.'_H1'}); $i++)
	{
		my $c = substr($sample_haplotypes{$sampleID.'_H1'}, $i, 1);
		if($c ne '_')
		{
			$non_gap_symbols_H1++;
		}
		push(@{$sample_haplotype_positions{$sampleID.'_H1'}}, $non_gap_symbols_H1);
	}
	
	die unless(scalar(@loci_in_graph) == scalar(@{$sample_haplotype_positions{$sampleID.'_H1'}}));
	
	my $non_gap_symbols_H2 = 0;
	for(my $i = 0; $i < length($sample_haplotypes{$sampleID.'_H2'}); $i++)
	{
		my $c = substr($sample_haplotypes{$sampleID.'_H2'}, $i, 1);
		if($c ne '_')
		{
			$non_gap_symbols_H2++;
		}
		push(@{$sample_haplotype_positions{$sampleID.'_H2'}}, $non_gap_symbols_H2);
	}	
	
	die unless(scalar(@loci_in_graph) == scalar(@{$sample_haplotype_positions{$sampleID.'_H2'}}));
	

	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H1'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions/'.$sampleID.'_PnP2Viterbi_H1.txt';
	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H2'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions/'.$sampleID.'_PnP2Viterbi_H2.txt';
	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_PGF'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions/'.$sampleID.'_PnP2Viterbi_PGF.txt';
	
	my $fh1, my $fh2, my $fh3;
	open($fh1, '>', $sample_haplotype_insertion_IGV_filenames{$sampleID.'_H1'}) or die "Cannot open ".$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H1'};
	open($fh2, '>', $sample_haplotype_insertion_IGV_filenames{$sampleID.'_H2'}) or die "Cannot open ".$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H2'};
	open($fh3, '>', $sample_haplotype_insertion_IGV_filenames{$sampleID.'_PGF'}) or die "Cannot open ".$sample_haplotype_insertion_IGV_filenames{$sampleID.'_PGF'};
	
	
	$sample_haplotype_insertion_IGV_filehandles{$sampleID.'_H1'} = $fh1;
	$sample_haplotype_insertion_IGV_filehandles{$sampleID.'_H2'} = $fh2;
	$sample_haplotype_insertion_IGV_filehandles{$sampleID.'_PGF'} = $fh3;
	
	print { $fh1 } qq(new
genome ${sampleID}_1
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_1.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions
);

	print { $fh2 } qq(new
genome ${sampleID}_2
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_2.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions
);

	print { $fh3 } qq(new
genome hg19
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\samtools_merged_xMHC.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions
);

	$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H1'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/NewVariants/'.$sampleID.'_PnP2Viterbi_H1.txt';
	$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H2'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/NewVariants/'.$sampleID.'_PnP2Viterbi_H2.txt';
	$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_PGF'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/NewVariants/'.$sampleID.'_PnP2Viterbi_PGF.txt';
	
	my $fh4, my $fh5, my $fh6;
	open($fh4, '>', $sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H1'}) or die "Cannot open ".$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H1'};
	open($fh5, '>', $sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H2'}) or die "Cannot open ".$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_H2'};
	open($fh6, '>', $sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_PGF'}) or die "Cannot open ".$sample_haplotype_novelSeq_IGV_filenames{$sampleID.'_PGF'};
	
	
	$sample_haplotype_novelSeq_IGV_filehandles{$sampleID.'_H1'} = $fh4;
	$sample_haplotype_novelSeq_IGV_filehandles{$sampleID.'_H2'} = $fh5;
	# $sample_haplotype_novelSeq_IGV_filehandles{$sampleID.'_PGF'} = $fh6;

	print { $fh4 } qq(new
genome ${sampleID}_1
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_1.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\NewVariants
);

	print { $fh5 } qq(new
genome ${sampleID}_2
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_2.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\NewVariants
);

	# print { $fh6 } qq(new
# genome hg19
# load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\samtools_merged_xMHC.bam
# snapshotDirectory C:\\Users\AlexanderDilthey\\Desktop\temp\\${sampleID}\\snapshots\\NewVariants
# );		
}

my @amended_haplotype_files = glob('../tmp/*.amendedHaplotypes');
warn "No haplotype files to analyze?" unless(@amended_haplotype_files);

my %amended_sample_haplotypes;
foreach my $file (@amended_haplotype_files)
{
	my $pattern = $graph.'_(\S+?)_'.$kMer_size_inference;
	next unless($file =~ /Z1/);
	
	next unless($file =~ /$pattern/);
	my $sampleID = $1;	
	
	my $haplotypes_aref = read_haplotypes($file, \@loci_in_graph, 1);
	$amended_sample_haplotypes{$sampleID.'_H1'} = $haplotypes_aref->[1][0];
	$amended_sample_haplotypes{$sampleID.'_H2'} = $haplotypes_aref->[1][1];
}

my $printing_threshold_length = 30;
my $nearby_genes_distance = 30000;

my $temp_output_summary = $temp_output_pgf_insertions.'/summary.txt';
open(SUMMARY, '>', $temp_output_summary) or die "Cannot open $temp_output_summary";
print SUMMARY join("\t", qw/i genomic_position length pgf_before pgf_after inSamples inWhichSamples affectedByModification affectedByModification_WhichSamples properlyModified_Known properlyModified_Known_WhichSamples properlyModified_Novel properlyModified_Novel_WhichSamples uniqueUnderlyingSequence_properlyModified_Novel uniqueUnderlyingSequence_properlyModified_Novel_WhichSamples/, 'nearbyGenes'.$nearby_genes_distance), "\n";

for(my $i = 0; $i <= $#pgf_insertions; $i++)
{
	my $insertID = $i;
	my $insertion_data = $pgf_insertions[$i];
	
	die unless($original_alignment_positions->{pgf});
	
	die unless($pgf_insertions[$i][0] > 0);
	my $pgf_position_before = $original_alignment_positions->{pgf}[$pgf_insertions[$i][0]-1] - 1;
	my $pgf_position_after = $original_alignment_positions->{pgf}[$pgf_insertions[$i][1]+1] - 1;

	my $pgf_position_before_genomic = $pgf_position_before + $pgf_start; 
	my $pgf_position_after_genomic =  $pgf_position_after + $pgf_start; 

	next if($PGF_within_HLA{$pgf_position_before} or $PGF_within_HLA{$pgf_position_after});
	my $insertion_length = $pgf_insertions[$i][1] - $pgf_insertions[$i][0];

	next if($insertion_length < $printing_threshold_length);
	
	print "\rPrint insertion $i";
	
	open(F, '>', $temp_output_pgf_insertions.'/insertion_'.$i.'.txt') or die;
	print F "Insertion #${i}\n";
	print F "Start in alignment: $pgf_insertions[$i][0] (genomic ", $pgf_start+$pgf_position_before, ")\n";
	print F "Stop in alignment: $pgf_insertions[$i][1]\n";	
	
	print F "PGF position before: $pgf_position_before\n";
	print F "PGF position after: $pgf_position_after \n";	

	foreach my $haplotype (sort keys %{$original_alignment})
	{
		my $haplotype_string = $original_alignment->{$haplotype};

		my $insertion_sequence = substr($haplotype_string, $pgf_insertions[$i][0], $pgf_insertions[$i][1] - $pgf_insertions[$i][0] + 1);
		my $pre_extraction = substr($haplotype_string, $pgf_insertions[$i][0]-15, 15);
		my $post_extraction = substr($haplotype_string, $pgf_insertions[$i][1]+1, 15);
				
		print F sprintf("%20s", $haplotype), "\t", $pre_extraction, '  ', $insertion_sequence, '  ', $post_extraction, "\n";
	}
	
	die if($pgf_insertions[$i][0] == 0);
	die unless($pgf_position_before and $pgf_position_after);
	die Dumper('$pgf_position_after == ($pgf_position_before + 1)', $pgf_position_after, $pgf_position_before) unless($pgf_position_after == ($pgf_position_before + 1));	
	die if(substr($original_alignment->{pgf}, $pgf_insertions[$i][0]-1, 1) eq '_');
	die if(substr($original_alignment->{pgf}, $pgf_insertions[$i][1]+1, 1) eq '_');
	
	if($insertion_length >= $printing_threshold_length)
	{
		$long_insertions{$i} = $insertion_length;
	}
	
	# translate this to position in the graph alignment (which might differ from haplotype alignment due to HLA allele sequences)
	
	die Dumper($pgf_to_graphAlignment{$pgf_position_before}, $pgf_to_graphAlignment{$pgf_position_before-1}, $pgf_to_graphAlignment{$pgf_position_before+1}) unless(scalar(@{$pgf_to_graphAlignment{$pgf_position_before}}) > 1);
	my $first_insertion_position_graphAlignment = $pgf_to_graphAlignment{$pgf_position_before}[1];
	my $last_insertion_position_graphAlignment = $pgf_to_graphAlignment{$pgf_position_before}[$#{$pgf_to_graphAlignment{$pgf_position_before}}];
	
	die unless($first_insertion_position_graphAlignment and $last_insertion_position_graphAlignment);
	die unless($last_insertion_position_graphAlignment >= $first_insertion_position_graphAlignment);
	
	my $this_extraction_target = {
		i => $i,
		start_8hap_alignment => $pgf_insertions[$i][0], 
		stop_8hap_alignment => $pgf_insertions[$i][1],
		start_graph_alignment => $first_insertion_position_graphAlignment,
		stop_graph_alignment => $last_insertion_position_graphAlignment,
		pgf_position_before => $pgf_position_before, 
		pgf_position_after => $pgf_position_after,
	};
	
	unless(($last_insertion_position_graphAlignment - $first_insertion_position_graphAlignment) == ($pgf_insertions[$i][1] - $pgf_insertions[$i][0]))
	{
		warn "LENGTH DIFFERENCE in inferred positions from original alignment vs graph alignment!\n".Dumper($this_extraction_target);
	}
	
	print F "\n\nSample haplotypes: \n\n";
	
	my $insertion_appears = 0;
	my $insertion_area_modified = 0;	
	my $insertion_appears_modified = 0;	
	my $insertion_appears_nonRef_modified = 0;
	my $insertion_uniqueUnderlyingSequence_nonRefModification = 0;
	
	my @haplotypeIDs_insertion_appears;
	my @haplotypeIDs_insertion_area_modified;		
	my @haplotypeIDs_insertion_appears_modified;	
	my @haplotypeIDs_insertion_appears_nonRef_modified;
	my @haplotypeIDs_uniqueUnderlyingSequence_nonRefModification;

	foreach my $haplotypeID (keys %sample_haplotypes)
	{
		my $haplotype = $sample_haplotypes{$haplotypeID};
		my $extraction_area = substr($haplotype, $first_insertion_position_graphAlignment, $last_insertion_position_graphAlignment - $first_insertion_position_graphAlignment + 1);
		my $pre_extraction = substr($haplotype, $first_insertion_position_graphAlignment-15, 15);
		my $post_extraction = substr($haplotype, $last_insertion_position_graphAlignment+1, 15);
		
		die unless(exists $sample_haplotype_positions{$haplotypeID}[$first_insertion_position_graphAlignment]);
		die unless(exists $sample_haplotype_positions{$haplotypeID}[$last_insertion_position_graphAlignment]);
		
		my $start_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haplotypeID}[$first_insertion_position_graphAlignment] + $remapping_flank_length + 1;
		my $stop_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haplotypeID}[$last_insertion_position_graphAlignment] + $remapping_flank_length + 1;
		
		
		
		my $this_sample_counts_as_appearing = 0;
		if($extraction_area =~ /[^\*_]/)
		{
			$insertion_appears++;
			push(@haplotypeIDs_insertion_appears, $haplotypeID.'['.$start_insertion_in_originalHaplotype_IGV.'-'.$stop_insertion_in_originalHaplotype_IGV.']');
			$this_sample_counts_as_appearing = 1;
			
			if($insertion_length < 2500)
			{
				my $filename_for_IGV_snapshot_PGF = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_PGF.png";
				my $pgf_haploID = $haplotypeID;
				$pgf_haploID =~ s/_H[12]/_PGF/;
				
				my $genomic_left_position_IGV= $pgf_position_before_genomic - 200;
				my $genomic_right_position_IGV= $pgf_position_after_genomic + 200;
				
				print { $sample_haplotype_insertion_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${genomic_left_position_IGV}-${genomic_right_position_IGV}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
				
				my $startPos_for_display = $start_insertion_in_originalHaplotype_IGV - 100;
				my $stopPos_for_display = $stop_insertion_in_originalHaplotype_IGV + 100;
				if($haplotypeID =~ /H1/)
				{
					my $filename_for_IGV_snapshot_H1 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_for_display}-${stopPos_for_display}\nsnapshot ${filename_for_IGV_snapshot_H1}\n);										
				}		
				else
				{
					die unless($haplotypeID =~ /H2/);
					my $filename_for_IGV_snapshot_H2 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_for_display}-${stopPos_for_display}\nsnapshot ${filename_for_IGV_snapshot_H2}\n);										
				}
			}
			else
			{
			
				my $genomic_left_position_IGV= $pgf_position_before_genomic - 200;
				my $genomic_right_position_IGV= $pgf_position_after_genomic + 200;
							
				my $filename_for_IGV_snapshot_PGF = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_PGF.png";
				my $pgf_haploID = $haplotypeID;
				$pgf_haploID =~ s/_H[12]/_PGF/;
				print { $sample_haplotype_insertion_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${genomic_left_position_IGV}-${genomic_right_position_IGV}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
				
				my $startPos_1_for_display = $start_insertion_in_originalHaplotype_IGV - 400;
				my $stopPos_1_for_display = $start_insertion_in_originalHaplotype_IGV + 400;
				
				my $startPos_2_for_display = $stop_insertion_in_originalHaplotype_IGV - 400;
				my $stopPos_2_for_display = $stop_insertion_in_originalHaplotype_IGV + 400;


				if($haplotypeID =~ /H1/)
				{
					my $filename_1_for_IGV_snapshot_H1 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1_front.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_1_for_display}-${stopPos_1_for_display}\nsnapshot ${filename_1_for_IGV_snapshot_H1}\n);										
					
					my $filename_2_for_IGV_snapshot_H1 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1_back.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_2_for_display}-${stopPos_2_for_display}\nsnapshot ${filename_2_for_IGV_snapshot_H1}\n);															
				}		
				else
				{
					die unless($haplotypeID =~ /H2/);
					
					my $filename_1_for_IGV_snapshot_H2 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2_front.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_1_for_display}-${stopPos_1_for_display}\nsnapshot ${filename_1_for_IGV_snapshot_H2}\n);										
					
					my $filename_2_for_IGV_snapshot_H2 = "${insertID}_${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2_back.png";
					print { $sample_haplotype_insertion_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${startPos_2_for_display}-${stopPos_2_for_display}\nsnapshot ${filename_2_for_IGV_snapshot_H2}\n);															
				}			
			}			
		}
		
		my $have_insertion_ontop_of_nongap = 0;
		
		if($amended_sample_haplotypes{$haplotypeID})
		{
			my @amended_haplotype = @{$amended_sample_haplotypes{$haplotypeID}};
			die unless(($#amended_haplotype+1) == length($haplotype));
			my @original_extracted_split = split(//, $extraction_area);
			my @amended_extracted = @amended_haplotype[$first_insertion_position_graphAlignment .. $last_insertion_position_graphAlignment];
			die unless($#original_extracted_split == $#amended_extracted);
			
			my $have_modification = 0;
			my $properModification_changedCharacters = 0;
			my $properModification_noPGFkMer = 0;
			
			my @properModification_positions_originalHaplotype;
			my @uniqueModifiedkMers;
		
			my @have_new_modifications;
			my @print_amended_extracted;
			
			for(my $i = 0; $i <= $#amended_extracted; $i++)
			{
				my $position_in_graph_alignment = $first_insertion_position_graphAlignment + $i;
				my $position_in_graph_alignment_2_IGV = $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment] + $remapping_flank_length;
				
				if($original_extracted_split[$i] ne $amended_extracted[$i])
				{				
					$have_modification = 1;
					my $thisI_have_insertion_ontop_of_nongap = 0;
					if(($original_extracted_split[$i] !~ /_/) and ($amended_extracted[$i] !~ /\*/))
					{
						$have_insertion_ontop_of_nongap = 1;
						$thisI_have_insertion_ontop_of_nongap = 1;
						$properModification_changedCharacters++;
						
						die Dumper("Undefined original position in haplotype $haplotypeID for alignment position $position_in_graph_alignment", $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment-1], $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment], $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment+1]) unless(exists $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment]);
						my $position_in_original_haplotype = $sample_haplotype_positions{$haplotypeID}[$position_in_graph_alignment];
						push(@properModification_positions_originalHaplotype, $position_in_original_haplotype.'-'.$i.'-'.$original_extracted_split[$i].'-'.$amended_extracted[$i]);
						

					}
					
					if(not exists $graph_locus_alleles[$position_in_graph_alignment]{$amended_extracted[$i]})
					{	
						push(@have_new_modifications, $i);
						
						if($insertID == 42)
						{
							if($haplotypeID eq 'AA02O9Q_A2_H1')
							{
								print Dumper($haplotypeID, $i, $first_insertion_position_graphAlignment, $position_in_graph_alignment, $graph_locus_alleles[$position_in_graph_alignment], $graph_locus_alleles[$position_in_graph_alignment-1], $graph_locus_alleles[$position_in_graph_alignment+1]), "\n\n";
							}
						}		
					
						if($thisI_have_insertion_ontop_of_nongap)
						{
							my $surrounding_kMer;
							my $kMer_left_position = ($i >= 7) ? ($i - 7) : 0;
							my $kMer_right_position = (($#amended_extracted - $i - 7) > 0) ? ($i + 7) : ($#amended_extracted);
							die unless($kMer_right_position >= $kMer_left_position);
							die unless($kMer_left_position >= 0);
							die unless($kMer_right_position <= $#amended_extracted);
							$surrounding_kMer = join('', @original_extracted_split[$kMer_left_position .. $kMer_right_position]);

							if($surrounding_kMer =~ /^[ACGT]+$/)
							{
								my $surrounding_kMer_2 = join('', @original_extracted_split[$kMer_left_position .. $i - 1]).$amended_extracted[$i].join('', @original_extracted_split[$i + 1 .. $kMer_right_position]);
								die unless($surrounding_kMer_2 =~ /^[ACGT]+$/);
								
								die Dumper($surrounding_kMer, $surrounding_kMer_2) unless(length($surrounding_kMer_2) == length($surrounding_kMer));
								
								if((index($pgf_noGaps, $surrounding_kMer) == -1) and (index($pgf_noGaps, $surrounding_kMer_2) == -1))
								{
									$properModification_noPGFkMer = 1;
									push(@uniqueModifiedkMers, $surrounding_kMer.','.join('-', $kMer_left_position, $i, $kMer_right_position));			

									my $IGV_left_position = $position_in_graph_alignment_2_IGV - 100;
									my $IGV_right_position = $position_in_graph_alignment_2_IGV + 100;
								   
									if($haplotypeID =~ /H1/)
									{
										my $filename_for_IGV_snapshot_H1 = "${i}_${pgf_position_before_genomic}_H1.png";
										print { $sample_haplotype_novelSeq_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${IGV_left_position}-${IGV_right_position}\nsnapshot ${filename_for_IGV_snapshot_H1}\n);										
									}
									else
									{
										die unless($haplotypeID =~ /H2/);
										my $filename_for_IGV_snapshot_H2 = "${i}_${pgf_position_before_genomic}_H2.png";
										print { $sample_haplotype_novelSeq_IGV_filehandles{$haplotypeID} } qq(goto xMHCpseudoChromHaplotype:${IGV_left_position}-${IGV_right_position}\nsnapshot ${filename_for_IGV_snapshot_H2}\n);								
									}									
								}
								else
								{
									warn "Dimiss variant\n\tkMer 1: $surrounding_kMer index: ".index($pgf_noGaps, $surrounding_kMer)."\n\tkMer 2: $surrounding_kMer_2 index: ".index($pgf_noGaps, $surrounding_kMer_2)."\n";
								}
							}										
						}
					}
				}
				
				my $max_length = (length($original_extracted_split[$i]) > length($amended_extracted[$i])) ? length($original_extracted_split[$i]) : length($amended_extracted[$i]);
				$original_extracted_split[$i] = sprintf("%".$max_length."s", $original_extracted_split[$i]);
				$amended_extracted[$i] = sprintf("%".$max_length."s", $amended_extracted[$i]);
				
				if($original_extracted_split[$i] ne $amended_extracted[$i])
				{
					push(@print_amended_extracted, $amended_extracted[$i]);
					# $insertion_appears_nonRef_modified = 1;
					# warn "Weird - insertion counts as not-appearing, but still modified???" unless($this_sample_counts_as_appearing);
				}
				else
				{
					push(@print_amended_extracted, sprintf("%".$max_length."s", ''));				
				}
			}
			
			if($have_modification)
			{
				if(@have_new_modifications)
				{
					print F "\t[NON-REF (i.e. POP VAR GRAPH) MODIFICATIONS at positions ".join(', ', @have_new_modifications)."]\n";
					print F sprintf("%30s", $haplotypeID), "\t", $pre_extraction, "  ", $extraction_area, "  ", $post_extraction, "\n\n";					
					print F sprintf("%30s", ''), "\t", "\t", join(' ', @original_extracted_split), "\n";
					print F sprintf("%30s", ''), "\t", "\t", join(' ', @print_amended_extracted), "\n";
					
					$insertion_area_modified++;
					push(@haplotypeIDs_insertion_area_modified, $haplotypeID);
					
					if($have_insertion_ontop_of_nongap)
					{
						$insertion_appears_nonRef_modified++;
						my $haplotypeID_withInfo = $haplotypeID.' ['.$properModification_changedCharacters.': '.join(', ', @properModification_positions_originalHaplotype).']';
						push(@haplotypeIDs_insertion_appears_nonRef_modified, $haplotypeID_withInfo);
						
						if($properModification_noPGFkMer)
						{
							$insertion_uniqueUnderlyingSequence_nonRefModification++;							
							push(@haplotypeIDs_uniqueUnderlyingSequence_nonRefModification, $haplotypeID.' ['.join(', ', @uniqueModifiedkMers).']');						
						}
					}

				}
				else
				{
					print F "\t[modified]\n";
					print F sprintf("%30s", $haplotypeID), "\t", $pre_extraction, "  ", $extraction_area, "  ", $post_extraction, "\n\n";										
					print F sprintf("%30s", ''), "\t",  join(' ', @original_extracted_split), "\n";
					print F sprintf("%30s", ''), "\t",  join(' ', @print_amended_extracted), "\n";
					
					$insertion_area_modified++;
					push(@haplotypeIDs_insertion_area_modified, $haplotypeID);					

					if($have_insertion_ontop_of_nongap)
					{					
						$insertion_appears_modified++;
						my $haplotypeID_withInfo = $haplotypeID.' ['.$properModification_changedCharacters.': '.join(', ', @properModification_positions_originalHaplotype).']';
						push(@haplotypeIDs_insertion_appears_modified, $haplotypeID_withInfo);
					}
				}

			}
			else
			{
				print F "\t[not modified]\n";		
				print F sprintf("%30s", $haplotypeID), "\t", $pre_extraction, "  ", $extraction_area, "  ", $post_extraction, "\n";					
			}
		}
		else
		{
			print F "\t[no amended haplotype file present]\n";				
			print F sprintf("%30s", $haplotypeID), "\t", $pre_extraction, "  ", $extraction_area, "  ", $post_extraction, "\n";								
		}
		print F "\t\t\t(IGV ", $start_insertion_in_originalHaplotype_IGV, " - ", $stop_insertion_in_originalHaplotype_IGV, "\n\n";
	}
	
	push(@{$extraction_targets{$first_insertion_position_graphAlignment}}, $this_extraction_target);
	
	
	my @nearby_genes;
	foreach my $gene (keys %$genes_href)
	{
		my $start = $genes_href->{$gene}[0];
		my $stop = $genes_href->{$gene}[0];
		my $start_distance = abs($pgf_position_before - $start);
		my $stop_distance = abs($pgf_position_before - $stop);
		my $min_distance = ($start_distance < $stop_distance) ? $start_distance : $stop_distance;
		die unless($min_distance >= 0);
		if($min_distance < $nearby_genes_distance)
		{
			push(@nearby_genes, $gene.'-'.$min_distance);
		}
	}

	my @fields_for_summary = ($i, $pgf_start+$pgf_position_before, $insertion_length, $pgf_position_before, $pgf_position_after, $insertion_appears, join('; ', @haplotypeIDs_insertion_appears), $insertion_area_modified, join('; ', @haplotypeIDs_insertion_area_modified), $insertion_appears_modified, join('; ', @haplotypeIDs_insertion_appears_modified), $insertion_appears_nonRef_modified, join('; ', @haplotypeIDs_insertion_appears_nonRef_modified), $insertion_uniqueUnderlyingSequence_nonRefModification, join('; ', @haplotypeIDs_uniqueUnderlyingSequence_nonRefModification), join('; ', @nearby_genes));
	print SUMMARY join("\t", @fields_for_summary), "\n";
	
	print F "\n\nSUMMARY line:\n";	
	print F join("\t", @fields_for_summary), "\n";
	 
	close(F);

}     

print "\n\nIdentified ", scalar(keys %long_insertions), " long insertions:\n";
foreach my $i (sort {$long_insertions{$b} <=> $long_insertions{$a}} keys %long_insertions)
{
	print "$i: $long_insertions{$i}\n";
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
			my $line = $_;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;			
			my @this_line_fields = split(/ /, $line);
			shift(@this_line_fields); # kick out individual ID		
			for(my $j = 0; $j <= $#this_line_fields; $j++)
			{
				if($j == 0)
				{
					# print join(' ', $segment_file, $beginning_of_locus_alleles, $j, $beginning_of_locus_alleles + $j, $this_line_fields[$j]), "\n";
				}
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
	my $amendedHaplotypes = shift;
	
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
		my $lastIndex = $#fields - 2;
		if($amendedHaplotypes)
		{
			$lastIndex = $#fields;
		}
		my $haplotype = join('', @fields[2 .. $lastIndex]);
		if($. == 1)
		{
			$haplotype_1 = $haplotype;
			@haplotype_1_alignment_fields = @fields[2 .. $lastIndex];
		}
		elsif($. == 2)
		{
			$haplotype_2 = $haplotype;	
			@haplotype_2_alignment_fields = @fields[2 .. $lastIndex];
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

