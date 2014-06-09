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
my @_levels_chromotype_sequence_to_align_to_withGaps = ($level_start .. $level_stop);
die unless(length($chromotype_sequence_to_align_to) == scalar(@_levels_chromotype_sequence_to_align_to_withGaps));
my @levels_chromotype_sequence_to_align;
die if(substr($chromotype_sequence_to_align_to, 0, 1) eq '_');
die if(substr($chromotype_sequence_to_align_to, length($chromotype_sequence_to_align_to)-1, 1) eq '_');
for(my $i = 0; $i < length($chromotype_sequence_to_align_to); $i++)
{
	my $C = substr($chromotype_sequence_to_align_to, $i, 1);
	if($C ne '_')
	{
		push(@levels_chromotype_sequence_to_align, $_levels_chromotype_sequence_to_align_to_withGaps[$i]);
	}
}
$chromotype_sequence_to_align_to =~ s/_//g;
die unless(length($chromotype_sequence_to_align_to) == scalar(@levels_chromotype_sequence_to_align));

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
		
		my $totalSequenceGaps = 0;
		my $totalGraphGaps = 0;
		die unless(length($fasta{$sequenceKey}) == length($fasta{$graphKey}));
		for(my $i = 0; $i < length($fasta{$sequenceKey}); $i++)
		{
			my $C_G = substr($fasta{$graphKey}, $i, 1);
			my $C_S = substr($fasta{$sequenceKey}, $i, 1);
			my $level = $levels[$i];
			
			if($C_G eq '_')
			{
				if($level == -1)
				{
					die if ($C_S eq '_');
					$totalGraphGaps++;
				}
			}
			elsif($C_S eq '_')
			{
				if($C_G ne '_')
				{
					$totalSequenceGaps++;
				}
			}
		}
		my $sequence_noGaps = $fasta{$sequenceKey};
		$sequence_noGaps =~ s/_//g;
		my $alignedSequenceLength = length($sequence_noGaps);

		
		my $dismiss = 0;
		if($totalGraphGaps >= (0.5* $alignedSequenceLength))
		{
			$dismiss = 1;
		}
		if($totalSequenceGaps >= 150000)
		{
			$dismiss = 1;
		}
		next if($dismiss > 200000);		
		
		
		my $alignmentLength = $alignmentLastLevel - $alignmentFirstLevel + 1;
		
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
		next if($reconstructedGraphSequence =~ /^_+$/);
		die unless(length($reconstructedGraphSequence) == scalar(@alignmentIndices_overlap_noGraphGaps));
		
		# my $extractedContigSequence = join('', map {die unless($sequence_characters[$_]); $sequence_characters[$_]} @alignmentIndices_overlap_noGraphGaps);
		
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
		
		# if($align_to_chromotype == 1)
		# {
			# die Dumper($equivalentChromotypeSequence_1, $reconstructedGraphSequence)
		# }
		# else
		# {
			# die Dumper($equivalentChromotypeSequence_2, $reconstructedGraphSequence)	
		# }		
				
		my $image_fn = "../tmp/dotplots/".$contigID.".png";
		
		my $reconstructedGraphSequence_noGaps = $reconstructedGraphSequence;
		$reconstructedGraphSequence_noGaps =~ s/_//g;
		
		#my $extractedContigSequence_noGaps = $extractedContigSequence;
		#$extractedContigSequence_noGaps =~ s/_//g;
		  
		my $plotTitle = qq($contigID [$method graph:$levels[$alignmentIndices_overlap_noGraphGaps[0]] - $levels[$alignmentIndices_overlap_noGraphGaps[-1]] contig(alignment):$alignmentIndices_overlap_noGraphGaps[0]-$@alignmentIndices_overlap_noGraphGaps[-1]]);

		# plot complete contig sequence
		my $completeContigSequence_noGaps = join('', @sequence_characters);
		die unless(length($completeContigSequence_noGaps) == scalar(@levels));
		my @levels_completeContigSequence_noGaps;
		for(my $i = 0; $i < length($completeContigSequence_noGaps); $i++)
		{
			my $C = substr($completeContigSequence_noGaps, $i, 1);
			if($C ne '_')
			{
				push(@levels_completeContigSequence_noGaps, $levels[$i]);
			}		
		}
		$completeContigSequence_noGaps =~ s/_//g;		
		die unless(length($completeContigSequence_noGaps) == scalar(@levels_completeContigSequence_noGaps));
		produce_dotPlot($chromotype_sequence_to_align_to, \@levels_chromotype_sequence_to_align, $completeContigSequence_noGaps, \@levels_completeContigSequence_noGaps, 10, $output_fh, $plotTitle, $image_fn, );
		
		$found_contigs++;
	}
	
	print "\tPrinted data for $found_contigs contigs.\n\n";
}

close($output_fh);

sub produce_dotPlot
{   
	my $x_string = shift;
	my $x_string_levels_aref = shift;
	my $y_string = shift;
	my $y_string_levels_aref = shift;
	my $k = shift;
	my $output_fh = shift;	
	my $plotTitle = shift;
	my $image_fn = shift;
	
	die if($x_string =~ /_/);
	die if($y_string =~ /_/);
	die unless(scalar(@{$x_string_levels_aref}) == length($x_string));
	die unless(scalar(@{$y_string_levels_aref}) == length($y_string));
	
	my %kMers_x;
	my %kMers_y;
	
	for(my $kMerI_x = 0; $kMerI_x < (length($x_string) - $k + 1); $kMerI_x++)
	{
		my $kMer = &uniqueMer(substr($x_string, $kMerI_x, $k));
		die unless(length($kMer) == $k);
		my @levels = @{$x_string_levels_aref}[$kMerI_x .. ($kMerI_x + $k - 1)];
		die unless(scalar(@levels) == $k);		
		push(@{$kMers_x{$kMer}}, [$kMerI_x, \@levels]);		
	}
	
	my @y_levels_in_alignment_andNoGap = grep {
			my $level = $y_string_levels_aref->[$_];
			die unless(defined $level);
			my $in = (($level != -1) and ($level >= $level_start) and ($level <= $level_stop));
			$in
		}
		(0 .. (length($y_string) - 1));
		
	my $first_y_index_in = $y_levels_in_alignment_andNoGap[0];
	my $last_y_index_in = $y_levels_in_alignment_andNoGap[-1];
	
	for(my $kMerI_y = 0; $kMerI_y < (length($y_string) - $k + 1); $kMerI_y++)
	{
		my $kMer = &uniqueMer(substr($y_string, $kMerI_y, $k));
		die unless(length($kMer) == $k);
		my @levels = (@{$y_string_levels_aref})[$kMerI_y .. ($kMerI_y + $k - 1)];
		die unless(scalar(@levels) == $k);		
		
		my $inAlignment = 0;
		my $firstYindex = $kMerI_y;
		my $lastYindex = $kMerI_y + $k - 1;
		if(($firstYindex >= $first_y_index_in) and ($firstYindex <= $last_y_index_in) and ($lastYindex >= $first_y_index_in) and ($lastYindex <= $last_y_index_in) )
		{
			$inAlignment = 1;
		}
		push(@{$kMers_y{$kMer}}, [$kMerI_y, \@levels, $inAlignment]);		
	}	
	
	my %kMers_unified = ((map {$_ => 1} keys %kMers_x), (map {$_ => 1} keys %kMers_y));
	
	print "produce_dotPlot(..): ", scalar(keys %kMers_unified), " combined kMers\n";
	
	my @x_dots;
	my @y_dots;
	my @dot_colours;
	
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
					my $X = $coordinates_x[$i][0];
					my $Y = $coordinates_y[$j][0];
					

					my $colour = 'black';
					if(defined $y_string_levels_aref)
					{
						my @Y_levels = @{$coordinates_y[$j]};
						
						my $Y_firstLevel = $coordinates_y[$j][1][0];
						my $Y_lastLevel = $coordinates_y[$j][1][-1];					
											
						if($coordinates_y[$j][2])
						{
							$colour = 'green';
						}
						
						if(defined $x_string_levels_aref)
						{
							if(join('|', @{$coordinates_x[$i][1]}) eq join('|', @{$coordinates_y[$j][1]}))
							{
								$colour = 'blue';
							}
						}
					}
					push(@x_dots, $X);
					push(@y_dots, $Y);
					push(@dot_colours, $colour);
				}
			}	
		}
	}
	
	my @contig_x_dots;
	my @contig_y_dots;	
	foreach my $kMer (keys %kMers_y)
	{
		my @coordinates_y = @{$kMers_y{$kMer}};
		
		for(my $i = 0; $i <= $#coordinates_y; $i++)
		{
			for(my $j = 0; $j <= $#coordinates_y; $j++)
			{
				my $X = $coordinates_y[$i][0];
				my $Y = $coordinates_y[$j][0];

				push(@contig_x_dots, $X);
				push(@contig_y_dots, $Y);
			}
		}	
	}
		
		
		
	my @chromotype_x_dots;
	my @chromotype_y_dots;	
	foreach my $kMer (keys %kMers_x)
	{
		my @coordinates_x = @{$kMers_x{$kMer}};
		
		for(my $i = 0; $i <= $#coordinates_x; $i++)
		{
			for(my $j = 0; $j <= $#coordinates_x; $j++)
			{
				my $X = $coordinates_x[$i][0];
				my $Y = $coordinates_x[$j][0];

				push(@chromotype_x_dots, $X);
				push(@chromotype_y_dots, $Y);
			}
		}	
	}
	
	my $length_x_string = length($x_string);
	
	print {$output_fh} 'xV <- c(', join(', ', @x_dots) ,')', "\n";
	print {$output_fh} 'yV <- c(', join(', ', @y_dots) ,')', "\n";
	print {$output_fh} 'contig_xV <- c(', join(', ', @contig_x_dots) ,')', "\n";
	print {$output_fh} 'contig_yV <- c(', join(', ', @contig_y_dots) ,')', "\n";
	print {$output_fh} 'chromotype_xV <- c(', join(', ', @chromotype_x_dots) ,')', "\n";
	print {$output_fh} 'chromotype_yV <- c(', join(', ', @chromotype_y_dots) ,')', "\n";
		
	print {$output_fh} 'dotC <- c(', join(', ', map {qq("$_")} @dot_colours) ,')', "\n";
	print {$output_fh} qq(filename <- "$image_fn"), "\n";	
	print {$output_fh} qq(png(filename = filename, width = 2250, height = 750)), "\n";	

	print {$output_fh} qq(layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), widths=c(1,1,1), heights=c(1))), "\n";
	print {$output_fh} qq(par(mar = c(5, 4, 4, 2)+3.1)), "\n";

	print {$output_fh} qq(plot(xV, yV, type="p", col=dotC, xlim=c(0, $length_x_string), pch = 15, xlab = "Chromotype", ylab="Contig", main = "Contig / Chromotype", cex.main = 3, cex.axis = 2.0, cex.lab = 2.5)), "\n";
	print {$output_fh} qq(legend("bottomright", legend = c("Actual alignment", "Alignment area", "Outside alignment area"), fill = c("blue", "green", "black"), cex = 2.5)), "\n";
	
	print {$output_fh} qq(plot(contig_xV, contig_yV, type="p", pch = 15, xlab = "Contig", ylab="Contig", main = "Contig / Contig", cex.main = 3, cex.axis = 2.0, cex.lab = 2.5)), "\n";	
	print {$output_fh} qq(plot(chromotype_xV, chromotype_yV, type="p", pch = 15, xlab = "Chromotype", ylab="Chromotype", main = "Chromotype / Chromotype", cex.main = 3, cex.axis = 2.0, cex.lab = 2.5)), "\n";		
	print {$output_fh} qq(dev.off()), "\n";
	
	print {$output_fh} "\n";
}

sub findAlignmentFiles
{
	my $method = shift;
	
	my $directory = '/Net/birch/data/dilthey/PnPHaploGraph2/tmp/alignedContigs/_GS_nextGen_varigraph3_AA02O9Q_Z2_31/contigs_xMHC_fasta/'.$method;
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
			print "\r$evaluationFile line $. ";
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
	
	return qq(/Net/birch/data/dilthey/PnPHaploGraph2/tmp/alignedContigs/_GS_nextGen_varigraph3_AA02O9Q_Z2_31/contigs_xMHC_fasta/chromotypeSupport_by_${method}.txt);

	# local testing	
	return qq(C:\\Users\\AlexanderDilthey\\Documents\\OxfordSVN\\documents\\analysis\\22 Oktober 2013\\contigSupportSpatially\\chromotypeSupport_by_toViterbiChromotypes.txt);
}


sub uniqueMer
{
	my $kMer = shift;
	my $kMer_complement = reverseComplement($kMer);
	if($kMer_complement le $kMer)
	{
		return $kMer_complement;
	}
	else
	{
		return $kMer;
	}	
}

sub reverseComplement
{
	my $in = shift;
	$in =~ tr/ACGTacgt/TGCAtgca/;
	return (reverse $in);
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