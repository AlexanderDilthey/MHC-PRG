#!/usr/bin/perl

use strict;
use Modern::Perl '2012';
use List::MoreUtils qw/mesh/;
use Data::Dumper;

$| = 1;

my $method = 'toViterbiChromotypes';
my $posInMiddleOfTarget = 32460000;

# my $level_start = 0;
# my $level_stop = 5000000;
# my $chromotype1 = ('A' x 5000000);  
my ($level_start, $level_stop, $reference_start, $reference_stop, $chromotype1, $chromotype2) = findChromotypesAndLevels($method, $posInMiddleOfTarget, 0);
die unless(length($chromotype1) == length($chromotype2));
print "Identified segment from $level_start to $level_stop (", ($level_stop - $level_start + 1), " levels, in reference space: $reference_start to $reference_stop)\n";

my $which_chromotype_onlyGaps = -1;   
my $align_to_chromotype;
die if((gapProportion($chromotype1) > 0.9) and (gapProportion($chromotype2) > 0.9));
if(gapProportion($chromotype1) > 0.9)
{
	$which_chromotype_onlyGaps = 1;
	$align_to_chromotype = 2;
}
if(gapProportion($chromotype2) > 0.9)
{
	$which_chromotype_onlyGaps = 2;
	$align_to_chromotype = 1;
}
die "No chromotype consists exclusively of gaps, which we would have expected!\n\nC1 = $chromotype1\n\nC2 = $chromotype2\n\n" unless($which_chromotype_onlyGaps != -1);

print "Chromotype $which_chromotype_onlyGaps consists of gaps entirely\n";

my $chromotype_sequence_to_align_to = ($which_chromotype_onlyGaps == 1) ? $chromotype2 : $chromotype1;
$chromotype_sequence_to_align_to =~ s/_//g;

# my $chromotype_sequence_to_align_to = $chromotype1;

my $output_filename = '../tmp/dotplots/'.$method.'.R';
my $output_fh;
open($output_fh, '>', $output_filename) or die;

my @contig_files = findAlignmentFiles($method);
my $found_contigs = 0;

foreach my $file (@contig_files)
{
	print $file, "\n";
	
	my %fasta;
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $firstLine = $_;
		chomp($firstLine);
		next unless($firstLine);
		
		my $secondLine = <F>;
		chomp($secondLine);
		
		die unless(substr($firstLine, 0, 1) eq '>');
		substr($firstLine, 0, 1) = '';
		
		$firstLine =~ s/ - score .+//;
		$fasta{$firstLine} = $secondLine;	
	}
	close(F);
	
	my %contigs = map {$_ =~ s/ - .+//; $_ => 1} keys %fasta;
	
	foreach my $contigID (keys %contigs)
	{
		my $levelsKey = $contigID.' - graph - levels';
		my $sequenceKey = $contigID.' - sequence';
		my $graphKey = $contigID.' - graph';
		
		die unless($fasta{$levelsKey});
		die unless($fasta{$sequenceKey});
		die Dumper("Cannot find key $graphKey in fasta hash", $contigID, $file, [keys %fasta]) unless($fasta{$graphKey});
		
		my $levels_string = $fasta{$levelsKey};   
		my @levels = split(/ /, $levels_string);
		
		if(scalar((grep{$_ ne '-1'} @levels)) == 0)
		{
			next;
		}
		
		my $alignmentFirstLevel = (grep{$_ ne '-1'} @levels)[0];
		my $alignmentLastLevel = (grep{$_ ne '-1'} @levels)[-1];
		die unless($alignmentLastLevel > $alignmentFirstLevel);
		die unless($alignmentFirstLevel >= 0);
		# die unless($alignmentFirstLevel < length($chromotype1));
		die unless($alignmentLastLevel >= 0);
		# die unless($alignmentLastLevel < length($chromotype1));
		
		# print $contigID, "\n";
		# print "\t", "First level: ", $alignmentFirstLevel, "\n";
		# print "\t", "Last level: ", $alignmentLastLevel, "\n";
		
		# is this contig interesting?
		unless((($alignmentFirstLevel <= $level_start) and ($alignmentLastLevel >= $level_stop)) or
		       (($alignmentFirstLevel >= $level_start) and ($alignmentFirstLevel <= $level_stop)) or
			   (($alignmentLastLevel >= $level_start) and ($alignmentLastLevel <= $level_stop))
		)
		{
			# print "\tDISMISS!\n";
			next;
			
		}
		
		my $alignmentLength = $alignmentLastLevel - $alignmentFirstLevel + 1;
		next if($alignmentLength > 200000);		
		
		# print $contigID, "\n";
		# print "\t", "First level: ", $alignmentFirstLevel, "\n";
		# print "\t", "Last level: ", $alignmentLastLevel, "\n";
		
		
		my @alignmentIndices = (0 .. $#levels);
		my @alignmentIndices_overlap = grep {($levels[$_] ne '-1') and ($levels[$_] >= $level_start) and ($levels[$_] <= $level_stop)} @alignmentIndices;		
		die unless($#alignmentIndices_overlap > -1);
		my $firstIndexOverlap = $alignmentIndices_overlap[0];
		my $lastIndexOverlap = $alignmentIndices_overlap[-1];
		
		# print "\tImmediately overlapping levels ", scalar(@alignmentIndices_overlap), ": from contig-alignment-index $firstIndexOverlap to $lastIndexOverlap \n";
		
		@alignmentIndices_overlap = ($firstIndexOverlap .. $lastIndexOverlap);
		# print "\tIncluding all gaps, overlapping levels: ", scalar(@alignmentIndices_overlap), "\n";
		
		my @graph_characters = split(//, $fasta{$graphKey});
		my @sequence_characters = split(//, $fasta{$sequenceKey});
		
		my @alignmentIndices_overlap_noGraphGaps = grep {$levels[$_] ne '-1'} @alignmentIndices_overlap;
		
		my $reconstructedGraphSequence = join('', map {die unless($graph_characters[$_]); $graph_characters[$_]} @alignmentIndices_overlap_noGraphGaps);
		my $extractedContigSequence = join('', map {die unless($sequence_characters[$_]); $sequence_characters[$_]} @alignmentIndices_overlap_noGraphGaps);
		next if($reconstructedGraphSequence =~ /^_+$/);
		
		die Dumper("Problem with first overlapping alignment index", $alignmentIndices_overlap_noGraphGaps[0], $level_start, $level_stop) unless(($levels[$alignmentIndices_overlap_noGraphGaps[0]] >= $level_start) and ($levels[$alignmentIndices_overlap_noGraphGaps[0]] <= $level_stop));
		die unless(($levels[$alignmentIndices_overlap_noGraphGaps[-1]] >= $level_start) and ($levels[$alignmentIndices_overlap_noGraphGaps[-1]] <= $level_stop));		
		
		my $reconstructedSequence_fromExtractedChromotypeLevel = $levels[$alignmentIndices_overlap_noGraphGaps[0]] - $level_start;
		my $reconstructedSequence_toExtractedChromotypeLevel = $levels[$alignmentIndices_overlap_noGraphGaps[-1]] - $level_start;
		die unless($reconstructedSequence_toExtractedChromotypeLevel >= $reconstructedSequence_fromExtractedChromotypeLevel);
		
		my $equivalentChromotypeSequence_1 = substr($chromotype1, $reconstructedSequence_fromExtractedChromotypeLevel, $reconstructedSequence_toExtractedChromotypeLevel - $reconstructedSequence_fromExtractedChromotypeLevel + 1);
		my $equivalentChromotypeSequence_2 = substr($chromotype2, $reconstructedSequence_fromExtractedChromotypeLevel, $reconstructedSequence_toExtractedChromotypeLevel - $reconstructedSequence_fromExtractedChromotypeLevel + 1);
		
		die unless(length($equivalentChromotypeSequence_1) == ($reconstructedSequence_toExtractedChromotypeLevel - $reconstructedSequence_fromExtractedChromotypeLevel + 1));
		die unless(length($equivalentChromotypeSequence_2) == ($reconstructedSequence_toExtractedChromotypeLevel - $reconstructedSequence_fromExtractedChromotypeLevel + 1));

		
		unless(($equivalentChromotypeSequence_1 eq $reconstructedGraphSequence) or ($equivalentChromotypeSequence_2 eq $reconstructedGraphSequence))
		{
			die Dumper(
				"Incongruence between contig-specified chromotype sequence and directly extracted chromotype sequence", 
				$contigID,
				$equivalentChromotypeSequence_1,
				$equivalentChromotypeSequence_2,
				$reconstructedGraphSequence
			);
		}
		
		if($align_to_chromotype == 1)
		{
			unless($equivalentChromotypeSequence_1 eq $reconstructedGraphSequence)
			{
				next;
			}
		}
		else
		{
			die unless($align_to_chromotype == 2);
			unless($equivalentChromotypeSequence_2 eq $reconstructedGraphSequence)
			{
				next;
			}			
		}
				
		my $image_fn = "../tmp/dotplots/".$contigID.".png";
		
		my $reconstructedGraphSequence_noGaps = $reconstructedGraphSequence;
		$reconstructedGraphSequence_noGaps =~ s/_//g;
		
		my $extractedContigSequence_noGaps = $extractedContigSequence;
		$extractedContigSequence_noGaps =~ s/_//g;
		
		my $plotTitle = qq($contigID [$method graph:$levels[$alignmentIndices_overlap_noGraphGaps[0]] - $levels[$alignmentIndices_overlap_noGraphGaps[-1]] contig(alignment):$alignmentIndices_overlap_noGraphGaps[0]-$@alignmentIndices_overlap_noGraphGaps[-1]]);

		# plot complete contig sequence
		my $completeContigSequence_noGaps = join('', @sequence_characters);
		$completeContigSequence_noGaps =~ s/_//g;		
		produce_dotPlot($chromotype_sequence_to_align_to, $completeContigSequence_noGaps, 10, $plotTitle, $image_fn, $output_fh);
		
		
		# plot only extracted contig sequence
		# produce_dotPlot($chromotype_sequence_to_align_to, $extractedContigSequence_noGaps, 5, $plotTitle, $image_fn, $output_fh);
		
		$found_contigs++;
	}
	
	print "\tPrinted data for $found_contigs contigs.\n\n";
}

close($output_fh);

sub produce_dotPlot
{   
	my $x_string = shift;
	my $y_string = shift;
	my $k = shift;
	my $plotTitle = shift;
	my $image_fn = shift;
	my $output_fh = shift;
	
	die if($x_string =~ /_/);
	die if($y_string =~ /_/);
	
	my %kMers_x;
	my %kMers_y;
	
	for(my $kMerI_x = 0; $kMerI_x < (length($x_string) - $k + 1); $kMerI_x++)
	{
		my $kMer = substr($x_string, $kMerI_x, $k);
		die unless(length($kMer) == $k);
		push(@{$kMers_x{$kMer}}, $kMerI_x);		
	}
	
	for(my $kMerI_y = 0; $kMerI_y < (length($y_string) - $k + 1); $kMerI_y++)
	{
		my $kMer = substr($y_string, $kMerI_y, $k);
		die unless(length($kMer) == $k);
		push(@{$kMers_y{$kMer}}, $kMerI_y);		
	}	
	
	my %kMers_unified = ((map {$_ => 1} keys %kMers_x), (map {$_ => 1} keys %kMers_y));
	
	print "produce_dotPlot(..): ", scalar(keys %kMers_unified), " combined kMers\n";
	
	my @x_dots;
	my @y_dots;
	
	foreach my $kMer (keys %kMers_unified)
	{
		if($kMers_x{$kMer} and $kMers_y{$kMer})
		{
			my @coordinates_x = @{$kMers_x{$kMer}};
			my @coordinates_y = @{$kMers_y{$kMer}};
			
			for(my $i = 0; $i <= $#coordinates_x; $i++)
			{
				for(my $j = 0; $j <= $#coordinates_y; $j++)
				{
					my $X = $coordinates_x[$i];
					my $Y = $coordinates_y[$j];
					push(@x_dots, $X);
					push(@y_dots, $Y);
				}
			}	
		}
	}
	
	my $length_x_string = length($x_string);
	print {$output_fh} 'xV <- c(', join(', ', @x_dots) ,')', "\n";
	print {$output_fh} 'yV <- c(', join(', ', @y_dots) ,')', "\n";
	print {$output_fh} qq(filename <- "$image_fn"), "\n";	
	print {$output_fh} qq(png(filename = filename, width = 1280, height = 1024)), "\n";
	print {$output_fh} qq(plot(xV, yV, type="p", xlim=c(0, $length_x_string), pch = 15, xlab = "Chromotype", ylab="Contig")), "\n";
	print {$output_fh} qq(dev.off()), "\n";
	
	print {$output_fh} "\n";
}

sub findAlignmentFiles
{
	my $method = shift;
	
	my $directory = '/Net/birch/data/dilthey/MHC-PRG/tmp/alignedContigs/_GS_nextGen_varigraph3_AA02O9Q_Z2_31/contigs_xMHC_fasta/'.$method;
	my @files = glob($directory.'/*.alignment');
	@files = grep {$_ !~ /\.status/} @files;

	print "findAlignmentFiles(..): Found ", scalar(@files), " files.\n";
	die unless(scalar(@files));
	
	return @files;
	
	# local testing
	return ('C:\Users\AlexanderDilthey\Desktop\temp\AA02O9Q_Z1\959.alignment');
}
sub findChromotypesAndLevels
{
	my $method = shift;
	my $referencePosInMiddle = shift;
	my $atLeastOneGap = shift;
	
	my $onlyDiploid = 1;
	
	print "\n";
	
	my $evaluationFile = &findChromotypeEvaluationFile($method);

	my $startSearchPos = -1;
	
	my @supportFile_level;	
	my @supportFile_referenceCoordinate;
	my @supportFile_diploid;
	my @supportFile_losePhase;
	my @supportFile_character_h1;
	my @supportFile_character_h2;
	open(EVALUATION, '<', $evaluationFile) or die "Cannot open $evaluationFile";
	my $headerLine = <EVALUATION>;
	chomp($headerLine);
	my @evaluation_header_fields = split(/\t/, $headerLine);
	die unless($evaluation_header_fields[0] eq 'Level');		
	die unless($evaluation_header_fields[1] eq 'ReferenceCoordinate');
	die unless($evaluation_header_fields[2] eq 'DiploidChromotype');
	die unless($evaluation_header_fields[3] eq 'ChromotypeLostPhase');
	die unless($evaluation_header_fields[4] eq 'ChromotypeCharacter_1');
	die unless($evaluation_header_fields[5] eq 'ChromotypeCharacter_2');
		
	while(<EVALUATION>)
	{
		if(($. % 100000) == 0)
		{
			print "\r$evaluationFile line $.";
		}
		
		my $line = $_;
		chomp($line);
		
		my @line_fields = split(/\t/, $line, 7);

		push(@supportFile_level, $line_fields[0]);
		push(@supportFile_referenceCoordinate, $line_fields[1]);
		push(@supportFile_diploid, $line_fields[2]);
		push(@supportFile_losePhase, $line_fields[3]);
		push(@supportFile_character_h1, $line_fields[4]);
		push(@supportFile_character_h2, $line_fields[5]);
		
		if($line_fields[1] == $referencePosInMiddle)
		{
			$startSearchPos = $#supportFile_referenceCoordinate;
		}
	}
	close(EVALUATION);
	
	if($startSearchPos == -1)
	{
		die "Cannot find middle-reference position $referencePosInMiddle"
	}
	
	die unless($supportFile_referenceCoordinate[$startSearchPos] == $referencePosInMiddle);
	
	if($onlyDiploid)
	{
		die "Position $referencePosInMiddle is not diploid" unless($supportFile_diploid[$startSearchPos] == 1);
	}
	
	print "\n";
	
	# now expand both sides
	my $start_i = $startSearchPos;
	my $stop_i = $startSearchPos;
	
	if($onlyDiploid)
	{
		while((($stop_i + 1) <= $#supportFile_level) and ($supportFile_diploid[$stop_i + 1] == 1) and ($supportFile_losePhase[$stop_i + 1] == 0) and ((not $atLeastOneGap) or (($supportFile_character_h1[$stop_i + 1] eq '_') or ($supportFile_character_h2[$stop_i + 1] eq '_'))))
		{
			$stop_i++;			
		}
		while((($start_i - 1) >= 0) and ($supportFile_diploid[$start_i - 1] == 1) and ($supportFile_losePhase[$start_i] == 0) and ((not $atLeastOneGap) or (($supportFile_character_h1[$start_i - 1] eq '_') or ($supportFile_character_h2[$start_i - 1] eq '_'))))
		{
			$start_i--;
		}		
	}
	else
	{
		die "Non-'only diploid' not implemented yet.";
	}
	
	my $start_level = $supportFile_level[$start_i];
	my $stop_level = $supportFile_level[$stop_i];
	my $reference_start = $supportFile_referenceCoordinate[$start_i];
	my $reference_stop = $supportFile_referenceCoordinate[$stop_i];	
	my $chromotype_1 = join('', @supportFile_character_h1[$start_i .. $stop_i]);
	my $chromotype_2 = join('', @supportFile_character_h2[$start_i .. $stop_i]);
	
	return ($start_level, $stop_level, $reference_start, $reference_stop, $chromotype_1, $chromotype_2);
}

sub findChromotypeEvaluationFile
{
	my $method = shift;
	
	return qq(/Net/birch/data/dilthey/MHC-PRG/tmp/alignedContigs/_GS_nextGen_varigraph3_AA02O9Q_Z2_31/contigs_xMHC_fasta/chromotypeSupport_by_${method}.txt);

	# local testing	
	return qq(C:\\Users\\AlexanderDilthey\\Documents\\OxfordSVN\\documents\\analysis\\22 Oktober 2013\\contigSupportSpatially\\chromotypeSupport_by_toViterbiChromotypes.txt);
}

sub gapProportion
{
	my $string = shift;
	if(length($string) == 0)
	{
		return 0;
	}
	
	my $string_noGaps = $string;
	$string_noGaps =~ s/_//g;
	
	my $number_gaps = length($string) - length($string_noGaps);
	return ($number_gaps / length($string));
}