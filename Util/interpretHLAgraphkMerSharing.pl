#!/usr/bin/perl
use Modern::Perl;
use List::MoreUtils qw/mesh all/;
use Data::Dumper;

my $kMer_sharing_file = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/graph_kMersPerLevel_HLAGRAPH_25.txt);
my $output_dir = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/kMerSharing);
mkdir($output_dir) unless (-e $output_dir);
my $plot_output_dir = $output_dir.'/plots';
mkdir($plot_output_dir) unless (-e $plot_output_dir);
	
my %data_per_locus;
my %kMers_2_levels;
my %kMers_in_reference;

# first round

print "Reading $kMer_sharing_file .. \n";
open(F, '<', $kMer_sharing_file) or die "Cannot open $kMer_sharing_file";

my $header_line = <F>;  
chomp($header_line);
my @header_fields = split(/\t/, $header_line);

die unless('Level' ~~ \@header_fields);
die unless('NucleotideGraph_Locus' ~~ \@header_fields);
die unless('kMerMultiplicityAtLevel' ~~ \@header_fields);
die unless('kMerMultiplicityInReferenceGenome' ~~ \@header_fields);

while(<F>)
{
	my $line = $_;
	chomp($line);
	my @line_fields = split(/\t/, $line);
	my %line_hash = (mesh @header_fields, @line_fields);
	
	next if($line_hash{'NucleotideGraph_Locus'} =~ /^before/);
	my $level = $line_hash{'Level'};
	die unless(defined $level);
	
	die "Can't parse $line_hash{'NucleotideGraph_Locus'}" unless($line_hash{'NucleotideGraph_Locus'} =~ /^(\w+)\//);
	my $locus = $1;
	
	
	die "Can't parse $line_hash{'NucleotideGraph_Locus'}" unless($line_hash{'NucleotideGraph_Locus'} =~ /^(\S+?)_/);
	my $locusWithInfo = $1;
	my @locusWithInfo_split = split(/\//, $locusWithInfo);
	die unless($#locusWithInfo_split >= 1);
	my $locus_part = join('_', @locusWithInfo_split);
	
	# normal kMers
	
	my %thisLevel_kMers;
	my $kMers_str = $line_hash{'kMerMultiplicityAtLevel'};
	my @kMers_str = split(/,/, $kMers_str);
	foreach my $kMer_str (@kMers_str)
	{
		my @parts = split(/\:/, $kMer_str);
		die unless($#parts == 1);
		my $kMer = $parts[0];
		my $count = $parts[1];
		
		push(@{$kMers_2_levels{$kMer}}, $level);
		push(@{$kMers_2_levels{reverseComplement($kMer)}}, $level);	

		$thisLevel_kMers{$kMer}++;
	}
	
	my $thisLine_locusData = [$level, $locus_part, \%thisLevel_kMers];
	push(@{$data_per_locus{$locus}}, $thisLine_locusData);
	
	# reference genome
	
	my $ref_kMers_str = $line_hash{'kMerMultiplicityInReferenceGenome'};
	if(defined $ref_kMers_str)
	{
		my @ref_kMers_str = split(/,/, $ref_kMers_str);
		foreach my $ref_kMer_str (@ref_kMers_str)
		{
			my @parts = split(/\:/, $ref_kMer_str);
			die unless($#parts == 1);
			my $ref_kMer = $parts[0];
			my $ref_count = $parts[1];
			if($ref_count > 0)
			{
				$kMers_in_reference{$ref_kMer} = 1;
				$kMers_in_reference{reverseComplement($ref_kMer)} = 1;
			}
		}
	}
}
close(F);

# second round

foreach my $locus (keys %data_per_locus)
{	
	my $output_file = $output_dir.'/'.$locus.'.txt';
	my $output_file_R = $output_dir.'/'.$locus.'.R';
	my $output_file_plot = $plot_output_dir . '/' .$locus.'.png';
	
	if(-e $output_file_plot)
	{
		unlink($output_file_plot) or die;
	}
	
	open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
	open(R, '>', $output_file_R) or die "Cannot open $output_file_R";
	
	my @forR_levels;
	my @forR_parts;
	
	my @forR_kMers_unique;
	my @forR_kMers_unique_inOtherLevels;
	my @forR_kMers_unique_inReference;
	my @forR_kMers_unique_OK;

	my @forR_kMers_weighted;
	my @forR_kMers_weighted_inOtherLevels;
	my @forR_kMers_weighted_inReference;
	my @forR_kMers_weighted_OK;
	
	print OUTPUT join("\t", qw/Level Part kMers_unweighted kMers_unweighted_inOtherLevels kMers_unweighted_inReference kMers_unweighted_OK kMers_weighted kMers_weighted_inOtherLevels kMers_weighted_inReference kMers_weighted_OK/), "\n";
	
	for(my $locusLevelI = 0; $locusLevelI <= $#{$data_per_locus{$locus}}; $locusLevelI++)
	{
		my $graphLevel = $data_per_locus{$locus}[$locusLevelI][0];
		my $part = $data_per_locus{$locus}[$locusLevelI][1];
		
		die unless(defined $graphLevel);
		
		my $kMers_unique = 0;
		my $kMers_unique_inOtherLevels = 0;
		my $kMers_unique_inReference = 0;
		my $kMers_unique_OK = 0;
		
		my $kMers_weighted = 0;
		my $kMers_weighted_inOtherLevels = 0;
		my $kMers_weighted_inReference = 0;		
		my $kMers_weighted_OK = 0;
		
		foreach my $kMer (keys %{$data_per_locus{$locus}[$locusLevelI][2]})
		{
			if(($kMer eq '_') or ($kMer eq '*'))
			{
				next;
			}
			
			my $count = $data_per_locus{$locus}[$locusLevelI][2]{$kMer};
			
			$kMers_unique++;			
			$kMers_weighted += $count;
			
			if($kMers_in_reference{$kMer})
			{
				$kMers_unique_inReference++;
				$kMers_weighted_inReference += $count;
			}
			
			die unless($kMers_2_levels{$kMer});
			my @levels_this_kMer = @{$kMers_2_levels{$kMer}};
			my @other_levels_this_kMer = grep {$_ != $graphLevel} @levels_this_kMer;
			
			if(scalar(@other_levels_this_kMer))
			{
				$kMers_unique_inOtherLevels++;
				$kMers_weighted_inOtherLevels += $count			
			}	
			
			if((not scalar(@other_levels_this_kMer)) and not($kMers_in_reference{$kMer}))
			{
				$kMers_unique_OK++;
				$kMers_weighted_OK++;
			}			
		}
				
		print OUTPUT join("\t",
			$graphLevel,
			$part,
			$kMers_unique,
			$kMers_unique_inOtherLevels,
			$kMers_unique_inReference,
			$kMers_unique_OK,
			$kMers_weighted,
			$kMers_weighted_inOtherLevels,
			$kMers_weighted_inReference,
			$kMers_weighted_OK,
		), "\n";	

		push(@forR_levels, $graphLevel);
		push(@forR_parts, $part);
		
		push(@forR_kMers_unique, $kMers_unique);
		push(@forR_kMers_unique_inOtherLevels, $kMers_unique_inOtherLevels);
		push(@forR_kMers_unique_inReference, $kMers_unique_inReference);
		push(@forR_kMers_unique_OK, $kMers_unique_OK);
		push(@forR_kMers_weighted, $kMers_weighted);
		push(@forR_kMers_weighted_inOtherLevels, $kMers_weighted_inOtherLevels);
		push(@forR_kMers_weighted_inReference, $kMers_weighted_inReference);
		push(@forR_kMers_weighted_OK, $kMers_weighted_OK);
	}

	my @forR_kMers_unique_inOtherLevels_prop;
	my @forR_kMers_unique_inReference_prop;
	my @forR_kMers_unique_OK_prop;

	my @forR_kMers_weighted_inOtherLevels_prop;
	my @forR_kMers_weighted_inReference_prop;
	my @forR_kMers_weighted_OK_prop;

	
	for(my $i = 0; $i <= $#forR_kMers_unique; $i++)
	{
		if($forR_kMers_unique[$i] == 0)
		{
			push(@forR_kMers_unique_inOtherLevels_prop, 0);
			push(@forR_kMers_unique_inReference_prop, 0);
			push(@forR_kMers_unique_OK_prop, 0);
		}
		else
		{
			push(@forR_kMers_unique_inOtherLevels_prop, $forR_kMers_unique_inOtherLevels[$i] / $forR_kMers_unique[$i]);
			push(@forR_kMers_unique_inReference_prop, $forR_kMers_unique_inReference[$i] / $forR_kMers_unique[$i]);
			push(@forR_kMers_unique_OK_prop, $forR_kMers_unique_OK[$i] / $forR_kMers_unique[$i]);
		}
		
		if($forR_kMers_weighted[$i] == 0)
		{
			push(@forR_kMers_weighted_inOtherLevels_prop, 0);
			push(@forR_kMers_weighted_inReference_prop, 0);
			push(@forR_kMers_weighted_OK_prop, 0);
		}
		else
		{
			push(@forR_kMers_weighted_inOtherLevels_prop, $forR_kMers_weighted_inOtherLevels[$i] / $forR_kMers_weighted[$i]);
			push(@forR_kMers_weighted_inReference_prop, $forR_kMers_weighted_inReference[$i] / $forR_kMers_weighted[$i]);
			push(@forR_kMers_weighted_OK_prop, $forR_kMers_weighted_OK[$i] / $forR_kMers_weighted[$i]);
		}		
	}
	
	
	print R 'locusName <- "', $locus, '"', "\n";

	print R 'graphLevels <- c(', join(', ', @forR_levels), ')', "\n";
	print R 'parts <- c(', join(', ', map {qq("$_")} @forR_parts), ')', "\n";
	
	print R 'kMers_unique <- c(', join(', ', @forR_kMers_unique), ')', "\n";
	print R 'kMers_unique_inOtherLevels <- c(', join(', ', @forR_kMers_unique_inOtherLevels), ')', "\n";
	print R 'kMers_unique_inReference <- c(', join(', ', @forR_kMers_unique_inReference), ')', "\n";
	print R 'kMers_unique_OK <- c(', join(', ', @forR_kMers_unique_OK), ')', "\n";

	print R 'kMers_weighted <- c(', join(', ', @forR_kMers_weighted), ')', "\n";
	print R 'kMers_weighted_inOtherLevels <- c(', join(', ', @forR_kMers_weighted_inOtherLevels), ')', "\n";
	print R 'kMers_weighted_inReference <- c(', join(', ', @forR_kMers_weighted_inReference), ')', "\n";
	print R 'kMers_weighted_OK <- c(', join(', ', @forR_kMers_weighted_OK), ')', "\n";	
	
	print R 'kMers_unique_inOtherLevels_prop <- c(', join(', ', @forR_kMers_unique_inOtherLevels_prop), ')', "\n";
	print R 'kMers_unique_inReference_prop <- c(', join(', ', @forR_kMers_unique_inReference_prop), ')', "\n";
	print R 'kMers_unique_OK_prop <- c(', join(', ', @forR_kMers_unique_OK_prop), ')', "\n";

	print R 'kMers_weighted_inOtherLevels_prop <- c(', join(', ', @forR_kMers_weighted_inOtherLevels_prop), ')', "\n";
	print R 'kMers_weighted_inReference_prop <- c(', join(', ', @forR_kMers_weighted_inReference_prop), ')', "\n";
	print R 'kMers_weighted_OK_prop <- c(', join(', ', @forR_kMers_weighted_OK_prop), ')', "\n";	
	
	print R qq(
	
png(filename = "$output_file_plot", width = 1280, height = 640)

par(mar = c(2, 2, 2, 1))
layout(matrix(1:2, nrow = 2, ncol = 1, byrow = TRUE))

pl <- function(kMers_unique, kMers_unique_OK, kMers_weighted_inOtherLevels, kMers_weighted_inReference, plotUnique, absOrRel) {
	cexText <- 0.7
	partLineMargin <- 20

	unique_parts <- unique(parts)
	part_coordinates <- list()
	for(part in unique_parts)
	{
		min_ind <- min(which(parts == part))
		max_ind <- max(which(parts == part))
		min_level <- graphLevels[[min_ind]]
		max_level <- graphLevels[[max_ind]]
		
		part_coordinates[[part]] <- c(min_level, max_level)
	}

	kMers_unique_max <- max(kMers_unique)
	yLim_plot <- kMers_unique_max * 1.30
	linePartPos <- kMers_unique_max * 1.15
	linePartPosUpI <- kMers_unique_max * 1.20
	linePartPosUpII <- kMers_unique_max * 1.25
	linePartPosDownI <- kMers_unique_max * 1.10
	linePartPosDownII <- kMers_unique_max * 1.05

	plotTitle <- paste(locusName, " ", absOrRel, sep = "")
	plot(x = graphLevels, y = kMers_unique, type = "l", col = "black", lwd = 1, ylim = c(0, yLim_plot), main = plotTitle, xlab = "")
	points(x = graphLevels, y = kMers_weighted_inOtherLevels, col = "orange", lwd = 1, lty = "solid", pch = 15, cex = 0.8)
	points(x = graphLevels, y = kMers_weighted_inReference, col = "red", lwd = 1, lty = "solid", pch = 15, cex = 1.0)
	# legend("topright", legend = c("Total kMers", "Fully usable", "In other levels of graph", "In non-covered reference parts"), fill = c("black", "green", "orange", "red"), cex = cexText)
	if(plotUnique == T)
	{
		points(x = graphLevels, y = kMers_unique_OK, type = "l", col = "green", lwd = 1.2)
	}

	pI <- 0
	for(part in names(part_coordinates))
	{
		pI <- pI + 1
		
		pData <- part_coordinates[[part]]
		
		line_left <- pData[[1]] + partLineMargin
		line_right <- pData[[2]] - partLineMargin
		xAlign <- 0
		
		if((pI %% 2) == 0)
		{
			yCoord <- linePartPosUpI
			if((pI %% 4) == 0)
			{
				yCoord <- linePartPosUpII
			}
			lines(c(line_left, line_right), c(linePartPos, linePartPos), col = "black")
			lines(c(line_left, line_left), c(linePartPos, yCoord), col = "black")
			text(x = line_left, y = yCoord, adj = c(xAlign, 0), labels = part, cex = cexText)
		} else
		{
			yCoord <- linePartPosDownI
			if(((pI + 1) %% 4) == 0)
			{
				yCoord <- linePartPosDownII
			}	
			lines(c(line_left, line_right), c(linePartPos, linePartPos), col = "black")
			lines(c(line_left, line_left), c(linePartPos, yCoord), col = "black")
			text(x = line_left, y = yCoord, adj = c(xAlign, 1), labels = part, cex = cexText)	
		}
	}
}


pl(kMers_unique, kMers_unique_OK, kMers_weighted_inOtherLevels, kMers_weighted_inReference, T, "absolute")

pl(rep(1, length(kMers_unique_OK_prop)), kMers_unique_OK_prop, kMers_weighted_inOtherLevels_prop, kMers_weighted_inReference_prop, F, "relative")

dev.off()

);
	
	
	close(OUTPUT);
	close(R);
	
	my $Rcmd = qq(Rscript $output_file_R);
	system($Rcmd);
	unless(-e $output_file_plot)
	{
		warn "Problem with plotting for $locus -- command $Rcmd failed.";
	}
}

print "\n\nOutput written to $output_dir\n\n";




sub reverseComplement
{
	my $in = shift;
	$in =~ tr/ACGTacgt/TGCAtgca/;
	return scalar(reverse $in);
}
