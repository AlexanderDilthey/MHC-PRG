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
use Storable;

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
my $printing_threshold_length = 50;

GetOptions (
	'graph:s' => \$graph,
	'original_alignment:s' => \$original_alignment,
);


die "Alignment file $original_alignment not existing" unless (-e $original_alignment);

my $graph_dir = '../tmp2/GS_nextGen/'.$graph;
die "Not existing: $graph_dir" unless(-e $graph_dir);

# my $genes_href = &get_chr6_genes();

my $graph_loci_data;
my @loci_in_graph;
my %pgf_to_graphAlignment;
my %PGF_within_HLA;
my $original_alignment_complete;

if(-e 'temp/data_graph_loci.data')
{
	warn "Use cached data! If you don't want this, delete files from temp/!";
	
	$graph_loci_data = retrieve 'temp/data_graph_loci.data';
	@loci_in_graph = @{retrieve 'temp/data_loci_in_graph.data'};
	%pgf_to_graphAlignment = %{retrieve 'temp/data_pgf_to_graphAlignment.data'};
	%PGF_within_HLA = %{retrieve 'temp/data_PGF_within_HLA.data'};
	$original_alignment_complete = retrieve 'temp/data_original_alignment_complete.data';
	
	print "\tloading data done\n";
}
else
{
	$graph_loci_data = get_graph_loci($graph_dir);
	@loci_in_graph = @{$graph_loci_data->[0]};

	# create index into graph PGF positions

	my $last_used_pgf = -1;
	my $first_position_PGF = -1;
	my $last_position_PGF = -1;

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

	$original_alignment_complete = read_alignment($original_alignment);

	store $graph_loci_data, 'temp/data_graph_loci.data';
	store \@loci_in_graph, 'temp/data_loci_in_graph.data';
	store \%pgf_to_graphAlignment, 'temp/data_pgf_to_graphAlignment.data';
	store \%PGF_within_HLA, 'temp/data_PGF_within_HLA.data';
	store $original_alignment_complete, 'temp/data_original_alignment_complete.data';	
}

my @graph_locus_alleles = @{$graph_loci_data->[2]};

my @graph_pos_2_pgf;
for(my $i = 0; $i <= $#loci_in_graph; $i++)
{
	my @parts = split(/_/, $loci_in_graph[$i]);
	die unless ($#parts == 2);
	my $pgfPos = $parts[2];	
	push(@graph_pos_2_pgf, $pgfPos);
}
	
my $original_alignment = $original_alignment_complete->{alignment};
my $original_alignment_positions = $original_alignment_complete->{positions};
my $original_alignment_positions_pgf = $original_alignment_positions->{'pgf'};
die unless(@$original_alignment_positions_pgf);
my @pgf_2_firstAlignmentColumn = (undef);

my $last_pgf_pos = 0;
for(my $i = 0; $i <= $#{$original_alignment_positions_pgf}; $i++)
{
	my $pgfPos = $original_alignment_positions_pgf->[$i];
	if($pgfPos != $last_pgf_pos)
	{
		die unless(($pgfPos - $last_pgf_pos) == 1);
		
		# print $pgfPos, "\t", $i, "\n";
		# die if ($#pgf_2_firstAlignmentColumn > 100);
		
		push(@pgf_2_firstAlignmentColumn, $i);
		$last_pgf_pos = $pgfPos;
	}
}

my $pgf_haplotype = $original_alignment->{'pgf'};
die unless($pgf_haplotype);
my $pgf_noGaps = $pgf_haplotype;
$pgf_noGaps =~ s/[^ACGTacgt]//g;
die if($pgf_noGaps =~ /\*/);
die if($pgf_noGaps =~ /\_/);


die Dumper(scalar(@pgf_2_firstAlignmentColumn), (length($pgf_noGaps)+1)) unless(scalar(@pgf_2_firstAlignmentColumn) == (length($pgf_noGaps)+1));



my $temp_output_pgf_insertions = 'temp/_identified_PGF_insertions_perSample';
unless(-e $temp_output_pgf_insertions)
{
	mkdir($temp_output_pgf_insertions) or die;
}
foreach my $existingFile (glob($temp_output_pgf_insertions.'/*'))
{
	unlink($existingFile) or die;
}

my $temp_output_pgf_deletions = 'temp/_identified_PGF_deletions_perSample';
unless(-e $temp_output_pgf_deletions)
{
	mkdir($temp_output_pgf_deletions) or die;
}
foreach my $existingFile (glob($temp_output_pgf_deletions.'/*'))
{
	unlink($existingFile) or die;
}

my @haplotype_files = glob('../tmp/*.viterbiHaplotypes');
die "No haplotype files to analyze?" unless(@haplotype_files);

my %sample_haplotypes;
my %sample_haplotype_positions;
my %sample_haplotype_insertion_IGV_filenames;
my %sample_haplotype_insertion_IGV_filehandles;
my %sample_haplotype_deletion_IGV_filenames;
my %sample_haplotype_deletion_IGV_filehandles;

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
	

	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H1'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions_perSample/'.$sampleID.'_PnP2Viterbi_H1.txt';
	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_H2'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions_perSample/'.$sampleID.'_PnP2Viterbi_H2.txt';
	$sample_haplotype_insertion_IGV_filenames{$sampleID.'_PGF'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Insertions_perSample/'.$sampleID.'_PnP2Viterbi_PGF.txt';
	
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
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions_perSample
);

	print { $fh2 } qq(new
genome ${sampleID}_2
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_2.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions_perSample
);
  
	print { $fh3 } qq(new
genome hg19
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\samtools_merged_xMHC.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Insertions_perSample
);

	$sample_haplotype_deletion_IGV_filenames{$sampleID.'_H1'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions_perSample/'.$sampleID.'_PnP2Viterbi_H1.txt';
	$sample_haplotype_deletion_IGV_filenames{$sampleID.'_H2'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions_perSample/'.$sampleID.'_PnP2Viterbi_H2.txt';
	$sample_haplotype_deletion_IGV_filenames{$sampleID.'_PGF'} = '/gpfs1/well/gsk_hla/Platypus/temp/IGVscripts/Deletions_perSample/'.$sampleID.'_PnP2Viterbi_PGF.txt';
	
	my $fh4, my $fh5, my $fh6;
	open($fh4, '>', $sample_haplotype_deletion_IGV_filenames{$sampleID.'_H1'}) or die "Cannot open ".$sample_haplotype_deletion_IGV_filenames{$sampleID.'_H1'};
	open($fh5, '>', $sample_haplotype_deletion_IGV_filenames{$sampleID.'_H2'}) or die "Cannot open ".$sample_haplotype_deletion_IGV_filenames{$sampleID.'_H2'};
	open($fh6, '>', $sample_haplotype_deletion_IGV_filenames{$sampleID.'_PGF'}) or die "Cannot open ".$sample_haplotype_deletion_IGV_filenames{$sampleID.'_PGF'};
	
	
	$sample_haplotype_deletion_IGV_filehandles{$sampleID.'_H1'} = $fh4;
	$sample_haplotype_deletion_IGV_filehandles{$sampleID.'_H2'} = $fh5;
	$sample_haplotype_deletion_IGV_filehandles{$sampleID.'_PGF'} = $fh6;
	
	print { $fh4 } qq(new
genome ${sampleID}_1
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_1.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions_perSample
);

	print { $fh5 } qq(new
genome ${sampleID}_2
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\mapped_xMHC_2.bam.exclusiveReads.sorted.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions_perSample
);

	print { $fh6 } qq(new
genome hg19
load C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\samtools_merged_xMHC.bam
snapshotDirectory C:\\Users\\AlexanderDilthey\\Desktop\\temp\\${sampleID}\\snapshots\\Deletions_perSample
);
}

my %insertions_perSample;
my %deletions_perSample;

foreach my $haploKey (keys %sample_haplotypes)
{
	my $haplotype = $sample_haplotypes{$haploKey};
	
	my $inDeletion = 0;
	my $deletionStart = 0;
	my $deletionStop = 0;
	
	my $inInsertion = 0;
	my $insertionStart = 0;
	my $insertionStop = 0;
	
	my @identified_insertions;
	my @identified_deletions;
	
	die unless(scalar(@graph_pos_2_pgf) == length($haplotype));	
	my $last_i_pgf;
	for(my $i = 0; $i < length($haplotype); $i++)
	{
		my $character = substr($haplotype, $i, 1);
		
		my $character_pgf_is_gap = ($last_i_pgf == $graph_pos_2_pgf[$i]);
		die if($character_pgf_is_gap and ($i == 0));
		
		die unless($character);

		if(($character_pgf_is_gap) and ($character ne '_'))
		{
			if((! $inInsertion) and ($character ne '_'))
			{			
				$inInsertion = 1;
				$insertionStart = $i;
			}
		}
		else
		{
			if($inInsertion)
			{
				$inInsertion = 0;
				$insertionStop = $i - 1;
				die unless($insertionStop >= $insertionStart);				
				push(@identified_insertions, [$insertionStart, $insertionStop]);
			}
		}
		
		if(($character eq '_') and (! $character_pgf_is_gap))
		{
			if((! $inDeletion) and (! $character_pgf_is_gap))
			{			
				$inDeletion = 1;
				$deletionStart = $i;
			}		
		}
		else
		{
			if($inDeletion)
			{
				$inDeletion = 0;
				$deletionStop = $i - 1;
				die unless($deletionStop >= $deletionStart);
				push(@identified_deletions, [$deletionStart, $deletionStop]);
			}		
		}
		
		$last_i_pgf = $graph_pos_2_pgf[$i];
	}
	
	if($inInsertion)
	{
		die unless($insertionStop >= $insertionStart);					
		push(@identified_insertions, [$insertionStart, length($haplotype)-1]);
	}
	
	if($inDeletion)
	{
		die unless($deletionStop >= $deletionStart);	
		push(@identified_deletions, [$deletionStart, length($haplotype)-1]);
	}
	
	@identified_insertions = grep {($_->[1] - $_->[0]) > $printing_threshold_length} @identified_insertions;
	@identified_deletions = grep {($_->[1] - $_->[0]) > $printing_threshold_length} @identified_deletions;

	print "${haploKey}: identified ", scalar(@identified_insertions), " insertions and ", scalar(@identified_deletions), " deletions\n";
	
	$insertions_perSample{$haploKey} = \@identified_insertions;
	$deletions_perSample{$haploKey} = \@identified_deletions;
}

foreach my $haploKey (keys %sample_haplotypes)
{
	open(INSERTIONS, '>', $temp_output_pgf_insertions.'/'.$haploKey.'.txt') or die;
	open(DELETIONS, '>', $temp_output_pgf_deletions.'/'.$haploKey.'.txt') or die;

	foreach my $insertionSpec (@{$insertions_perSample{$haploKey}})
	{
		my $pgf_position_before;
		my $pgf_position_after;
		
		my $pgf_position_before_genomic;
		my $pgf_position_after_genomic;
		
		my $first_insertion_position_graph = $insertionSpec->[0];
		my $last_insertion_position_graph = $insertionSpec->[1];

		my $first_insertion_position_alignment;
		my $last_insertion_position_alignment;
	
		my $start_insertion_in_originalHaplotype_IGV;
		my $stop_insertion_in_originalHaplotype_IGV;
		my $insertion_length;
		
		if($first_insertion_position_graph > 0)
		{
			# $pgf_position_before = $original_alignment_positions->{pgf}[$insertionSpec->[0]-1] - 1;			
			# $pgf_position_after = $original_alignment_positions->{pgf}[$insertionSpec->[1]+1] - 1;
			
			$pgf_position_before = $graph_pos_2_pgf[$first_insertion_position_graph-1];
			$pgf_position_after = $graph_pos_2_pgf[$last_insertion_position_graph+1];
			
			next if($PGF_within_HLA{$pgf_position_before} or $PGF_within_HLA{$pgf_position_after});		

			$pgf_position_before_genomic = $pgf_position_before + $pgf_start; 
			$pgf_position_after_genomic =  $pgf_position_after + $pgf_start; 			
			
			my $pgf_first_position = $graph_pos_2_pgf[$first_insertion_position_graph];
			my $pgf_last_position = $graph_pos_2_pgf[$last_insertion_position_graph];
			die unless($pgf_first_position == $pgf_last_position);
			
			$first_insertion_position_alignment = $pgf_2_firstAlignmentColumn[$pgf_first_position+1]+1;
			$last_insertion_position_alignment = $pgf_2_firstAlignmentColumn[$pgf_last_position+2]-1;
			die Dumper($insertionSpec, $pgf_first_position, [$first_insertion_position_alignment, $last_insertion_position_alignment], [@pgf_2_firstAlignmentColumn[$pgf_first_position - 5 .. $pgf_first_position + 5]]) unless($first_insertion_position_alignment <= $last_insertion_position_alignment);
			
			die unless(defined $sample_haplotype_positions{$haploKey}[$first_insertion_position_graph]);
			die unless(defined $sample_haplotype_positions{$haploKey}[$last_insertion_position_graph]);	
			
			$start_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haploKey}[$first_insertion_position_graph] + $remapping_flank_length + 1;
			$stop_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haploKey}[$last_insertion_position_graph] + $remapping_flank_length + 1;			
			
			$insertion_length = $insertionSpec->[1] - $insertionSpec->[0] + 1;
		}
		else
		{
		
			# $pgf_position_before = -1;		
			# $pgf_position_after = $original_alignment_positions->{pgf}[$insertionSpec->[1]+1] - 1;
			
			$pgf_position_before = $graph_pos_2_pgf[0];
			$pgf_position_after = $graph_pos_2_pgf[$last_insertion_position_graph+1];
						
			next if($PGF_within_HLA{$pgf_position_before} or $PGF_within_HLA{$pgf_position_after});		
						
			$pgf_position_before_genomic = $pgf_position_before + $pgf_start; 
			$pgf_position_after_genomic =  $pgf_position_after + $pgf_start; 			
		
			my $pgf_first_position = $graph_pos_2_pgf[$first_insertion_position_graph+1];
			my $pgf_last_position = $graph_pos_2_pgf[$last_insertion_position_graph];
			die unless($pgf_first_position == $pgf_last_position);
			
			$first_insertion_position_alignment = $pgf_2_firstAlignmentColumn[$pgf_first_position+1]+1;
			$last_insertion_position_alignment = $pgf_2_firstAlignmentColumn[$pgf_last_position+2]-1;
			die Dumper($insertionSpec, $pgf_first_position, [$first_insertion_position_alignment, $last_insertion_position_alignment], [@pgf_2_firstAlignmentColumn[$pgf_first_position - 5 .. $pgf_first_position + 5]]) unless($first_insertion_position_alignment <= $last_insertion_position_alignment);
			
			die unless(defined $sample_haplotype_positions{$haploKey}[$first_insertion_position_graph]);
			die unless(defined $sample_haplotype_positions{$haploKey}[$last_insertion_position_graph]);		
			
			$start_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haploKey}[$first_insertion_position_graph] + $remapping_flank_length + 1;
			$stop_insertion_in_originalHaplotype_IGV = $sample_haplotype_positions{$haploKey}[$last_insertion_position_graph] + $remapping_flank_length + 1;		
			
			$insertion_length = $insertionSpec->[1] - $insertionSpec->[0] + 1;		
		}
		next if($PGF_within_HLA{$pgf_position_before} or $PGF_within_HLA{$pgf_position_after});		
		
		die Dumper($insertionSpec, [$pgf_position_before, $pgf_position_after], [$pgf_position_before_genomic, $pgf_position_after_genomic], [$stop_insertion_in_originalHaplotype_IGV, $start_insertion_in_originalHaplotype_IGV]) unless($pgf_position_after_genomic >= $pgf_position_before_genomic);
		die Dumper($insertionSpec, [$pgf_position_before, $pgf_position_after], [$pgf_position_before_genomic, $pgf_position_after_genomic], [$stop_insertion_in_originalHaplotype_IGV, $start_insertion_in_originalHaplotype_IGV]) unless($stop_insertion_in_originalHaplotype_IGV >= $start_insertion_in_originalHaplotype_IGV);

		print INSERTIONS "INSERTION FROM $insertionSpec->[0] TO $insertionSpec->[1]\n====================================================\n";
		print INSERTIONS "Length: $insertion_length\n";
		print INSERTIONS "On PGF: $pgf_position_before_genomic - $pgf_position_after_genomic\n";
		print INSERTIONS "For IGV: $start_insertion_in_originalHaplotype_IGV - $stop_insertion_in_originalHaplotype_IGV\n";
		print INSERTIONS "In alignment: $first_insertion_position_alignment - $last_insertion_position_alignment\n\n";
		
		
		foreach my $haplotype (sort keys %{$original_alignment})
		{
			my $haplotype_string = $original_alignment->{$haplotype};

			my $insertion_sequence = substr($haplotype_string, $first_insertion_position_alignment, $last_insertion_position_alignment - $first_insertion_position_alignment + 1);
			my $pre_extraction = substr($haplotype_string, $first_insertion_position_alignment-15, 15);
			my $post_extraction = substr($haplotype_string, $last_insertion_position_alignment+1, 15);
					
			print INSERTIONS sprintf("%20s", $haplotype), "\t", $pre_extraction, '  ', $insertion_sequence, '  ', $post_extraction, "\n";
		}
	
		my $thisHaplotype_extraction_area = substr($sample_haplotypes{$haploKey}, $first_insertion_position_graph, $last_insertion_position_graph - $first_insertion_position_graph + 1);
		my $thisHaplotype_extraction_areapre_extraction = substr($sample_haplotypes{$haploKey}, $first_insertion_position_graph-15, 15);
		my $thisHaplotype_extraction_areapost_extraction = substr($sample_haplotypes{$haploKey}, $last_insertion_position_graph+1, 15);
		
		die Dumper($insertionSpec, $thisHaplotype_extraction_area) if($thisHaplotype_extraction_area =~ /_/);
		
		print INSERTIONS "\nSample haplotype:\n\n";		
		
		print INSERTIONS sprintf("%20s", $haploKey), "\t", $thisHaplotype_extraction_areapre_extraction, '  ', $thisHaplotype_extraction_area, '  ', $thisHaplotype_extraction_areapost_extraction, "\n\n";
			
		if($insertion_length < 2500)
		{
			my $filename_for_IGV_snapshot_PGF = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_PGF.png";
			my $pgf_haploID = $haploKey;
			$pgf_haploID =~ s/_H[12]/_PGF/;
			
			my $genomic_left_position_IGV = $pgf_position_before_genomic - 200;
			my $genomic_right_position_IGV = $pgf_position_after_genomic + 200;
			
			print { $sample_haplotype_insertion_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${genomic_left_position_IGV}-${genomic_right_position_IGV}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
			
			my $startPos_for_display = $start_insertion_in_originalHaplotype_IGV - 100;
			my $stopPos_for_display = $stop_insertion_in_originalHaplotype_IGV + 100;
			if($haploKey =~ /H1/)
			{
				my $filename_for_IGV_snapshot_H1 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_for_display}-${stopPos_for_display}\nsnapshot ${filename_for_IGV_snapshot_H1}\n);										
			}		
			else
			{
				die unless($haploKey =~ /H2/);
				my $filename_for_IGV_snapshot_H2 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_for_display}-${stopPos_for_display}\nsnapshot ${filename_for_IGV_snapshot_H2}\n);										
			}
		}
		else
		{
		
			my $genomic_left_position_IGV = $pgf_position_before_genomic - 200;
			my $genomic_right_position_IGV = $pgf_position_after_genomic + 200;
						
			my $filename_for_IGV_snapshot_PGF = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_PGF.png";
			my $pgf_haploID = $haploKey;
			$pgf_haploID =~ s/_H[12]/_PGF/;
			print { $sample_haplotype_insertion_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${genomic_left_position_IGV}-${genomic_right_position_IGV}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
			
			my $startPos_1_for_display = $start_insertion_in_originalHaplotype_IGV - 100;
			my $stopPos_1_for_display = $start_insertion_in_originalHaplotype_IGV + 100;
			
			my $startPos_2_for_display = $stop_insertion_in_originalHaplotype_IGV - 100;
			my $stopPos_2_for_display = $stop_insertion_in_originalHaplotype_IGV + 100;

			my $startPos_3_for_display = $start_insertion_in_originalHaplotype_IGV - 100;
			my $stopPos_3_for_display = $stop_insertion_in_originalHaplotype_IGV + 100;
			
			if($haploKey =~ /H1/)
			{
				my $filename_1_for_IGV_snapshot_H1 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1_front.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_1_for_display}-${stopPos_1_for_display}\nsnapshot ${filename_1_for_IGV_snapshot_H1}\n);										
				
				my $filename_2_for_IGV_snapshot_H1 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1_back.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_2_for_display}-${stopPos_2_for_display}\nsnapshot ${filename_2_for_IGV_snapshot_H1}\n);															
				
				my $filename_3_for_IGV_snapshot_H1 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H1_all.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_3_for_display}-${stopPos_3_for_display}\nsnapshot ${filename_3_for_IGV_snapshot_H1}\n);															
				
			}		   
			else
			{
				die unless($haploKey =~ /H2/);
				
				my $filename_1_for_IGV_snapshot_H2 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2_front.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_1_for_display}-${stopPos_1_for_display}\nsnapshot ${filename_1_for_IGV_snapshot_H2}\n);										
				
				my $filename_2_for_IGV_snapshot_H2 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2_back.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_2_for_display}-${stopPos_2_for_display}\nsnapshot ${filename_2_for_IGV_snapshot_H2}\n);															
				
				my $filename_3_for_IGV_snapshot_H2 = "${pgf_position_before_genomic}_${pgf_position_after_genomic}_H2_all.png";
				print { $sample_haplotype_insertion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${startPos_3_for_display}-${stopPos_3_for_display}\nsnapshot ${filename_3_for_IGV_snapshot_H2}\n);															
				
			}			
		}				
	}
	

	foreach my $deletionSpec (@{$deletions_perSample{$haploKey}})
	{
	

		my $deletionStart_Graph = $deletionSpec->[0];
		my $deletionStop_Graph = $deletionSpec->[1];
		
		
		die unless($original_alignment_positions->{pgf});
		
		my $deletionStart_PGF = $graph_pos_2_pgf[$deletionStart_Graph];
		my $deletionStop_PGF = $graph_pos_2_pgf[$deletionStop_Graph];
		
		my $deletionStart_Alignment = $pgf_2_firstAlignmentColumn[$deletionStart_PGF+1];
		my $deletionStop_Alignment = $pgf_2_firstAlignmentColumn[$deletionStop_PGF+1];
		die unless($deletionStop_Alignment >= $deletionStart_Alignment);

		
		my $deletionStart_Genomic = $pgf_start + $deletionStart_PGF;
		my $deletionStop_Genomic = $pgf_start + $deletionStop_PGF;
		die unless($deletionStart_Genomic <= $deletionStop_Genomic);
		
		next if($PGF_within_HLA{$deletionStart_PGF} or $PGF_within_HLA{$deletionStop_PGF});
		
		die unless($pgf_to_graphAlignment{$deletionStart_PGF} and $pgf_to_graphAlignment{$deletionStop_PGF});
		
		# my $deletionStart_Graph = $pgf_to_graphAlignment{$deletionStart_PGF}[0];
		# my $deletionStop_Graph = $pgf_to_graphAlignment{$deletionStop_PGF}[0];
		
		my $deletion_length = $deletionStop_Graph - $deletionStart_Graph;
		
		die unless($deletionStop_Graph >= $deletionStart_Graph);
		
		my $IGV_deletion_start = $sample_haplotype_positions{$haploKey}[$deletionStart_Graph] + $remapping_flank_length + 1;
		my $IGV_deletion_stop = $sample_haplotype_positions{$haploKey}[$deletionStop_Graph] + $remapping_flank_length + 1;
		
		die Dumper($deletionSpec, [$deletionStart_PGF, $deletionStop_PGF], [$deletionStart_Graph, $deletionStop_Graph], [$IGV_deletion_stop, $IGV_deletion_start]) unless($IGV_deletion_stop >= $IGV_deletion_start);
			
		print DELETIONS "DELETION FROM $deletionStart_Alignment TO $deletionStop_Alignment\n====================================================\n";
		print DELETIONS "Length: $deletion_length\n";
		print DELETIONS "On PGF: $deletionStart_Genomic - $deletionStop_Genomic\n";
		print DELETIONS "For IGV: $IGV_deletion_start - $IGV_deletion_stop\n\n";
					
		foreach my $haplotype (sort keys %{$original_alignment})
		{
			my $haplotype_string = $original_alignment->{$haplotype};

			my $insertion_sequence = substr($haplotype_string, $deletionStart_Alignment, $deletionStop_Alignment - $deletionStart_Alignment + 1);
			my $pre_extraction = substr($haplotype_string, $deletionStart_Alignment-15, 15);
			my $post_extraction = substr($haplotype_string, $deletionStop_Alignment+1, 15);
					
			print DELETIONS sprintf("%20s", $haplotype), "\t", $pre_extraction, '  ', $insertion_sequence, '  ', $post_extraction, "\n";
		}
	
		my $thisHaplotype_extraction_area = substr($sample_haplotypes{$haploKey}, $deletionStart_Graph, $deletionStop_Graph - $deletionStart_Graph + 1);
		my $thisHaplotype_extraction_areapre_extraction = substr($sample_haplotypes{$haploKey}, $deletionStart_Graph-15, 15);
		my $thisHaplotype_extraction_areapost_extraction = substr($sample_haplotypes{$haploKey}, $deletionStop_Graph+1, 15);
		
		die unless($thisHaplotype_extraction_area =~ /^_+$/);
		
		print DELETIONS "\nSample haplotype:\n\n";		
		print DELETIONS sprintf("%20s", $haploKey), "\t", $thisHaplotype_extraction_areapre_extraction, '  ', $thisHaplotype_extraction_area, '  ', $thisHaplotype_extraction_areapost_extraction, "\n\n";
	

		my $display_rightCoordinate_Genomic = $deletionStart_Genomic - 400;		
		my $display_leftCoordinate_Genomic = $deletionStart_Genomic + 400;
		my $filename_for_IGV_snapshot_PGF = "${deletionStart_Genomic}_PGF.png";
		my $pgf_haploID = $haploKey;
		$pgf_haploID =~ s/_H[12]/_PGF/;
		print { $sample_haplotype_deletion_IGV_filehandles{$pgf_haploID}  } qq(goto 6:${display_rightCoordinate_Genomic}-${display_leftCoordinate_Genomic}\nsnapshot ${filename_for_IGV_snapshot_PGF}\n);			
		
		if($haploKey =~ /H1/)
		{
			my $filename_for_IGV_snapshot_H1 = "${deletionStart_Genomic}_H1.png";
			print { $sample_haplotype_deletion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${IGV_deletion_start}\nsnapshot ${filename_for_IGV_snapshot_H1}\n);										
		}
		else
		{
			die unless($haploKey =~ /H2/);
			my $filename_for_IGV_snapshot_H2 = "${deletionStart_Genomic}_H2.png";
			print { $sample_haplotype_deletion_IGV_filehandles{$haploKey} } qq(goto xMHCpseudoChromHaplotype:${IGV_deletion_start}\nsnapshot ${filename_for_IGV_snapshot_H2}\n);								
		}
	}	
	
	close(INSERTIONS);
	close(DELETIONS);
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
	my $input = '/Net/birch/data/dilthey/PnPHaploGraph2/src/temp/UCSC_chr6_RefSeq_genes.txt';
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

