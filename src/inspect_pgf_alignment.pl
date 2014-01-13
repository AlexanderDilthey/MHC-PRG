#!/usr/bin/perl
use strict;
use Data::Dumper;

my %contigs;

my $alignment = qq(pgfContigs_for_testing.txt.aligned);
open(ALIGNMENT, '<', $alignment) or die "Cannot open $alignment";
my $currentContig;
while(<ALIGNMENT>)
{
	my $line = $_;
	chomp($line);
	if(substr($line, 0, 1) eq '>')
	{
		$currentContig = $line;
	}
	else
	{
		die unless($currentContig);
		$contigs{$currentContig} .= $line;
	}
}
close(ALIGNMENT);

my @contig_IDs = sort keys %contigs;

print "Have ", scalar(@contig_IDs), " contigs, \n";#, join("\n", map {' - -'.$_.'-'} @contig_IDs), "\n";

my @graph_contigs = grep {$_ =~ /\- graph \-/} @contig_IDs;

my $contigI = 0;
foreach my $graph_contig (@graph_contigs)
{
	print $graph_contig, "\n";
	$contigI++;
	
	my $sequence_contig = $graph_contig;
	$sequence_contig =~ s/ graph \-.+$/ sequence/;  
	#die Dumper($graph_contig, $sequence_contig);
	
	die "No sequence contig for -${graph_contig}- , would expect -${sequence_contig}- \n" unless(exists $contigs{$sequence_contig});
	
	my $sequence = $contigs{$sequence_contig};
	my $graph = $contigs{$graph_contig};
	
	die unless($sequence and $graph);
	
	my $sequence_startGap;
	if($sequence =~ /^(_+)/)
	{
		$sequence_startGap = $1;
	}
	
	my $sequence_stopGap;
	if($sequence =~ /(_+)$/)
	{
		$sequence_stopGap = $1;
	}
	
	my $sequence2 = $sequence;
	$sequence2 =~ s/_//g;
	die "Problem with sequence for contig $graph_contig!" unless(length($sequence2) > 0);
	
	print "Contig $sequence_contig, start gaps ", length($sequence_startGap), ", stop gaps ", length($sequence_stopGap), "\n";
	
	substr($sequence, 0, length($sequence_startGap)) = '';
	substr($graph, 0, length($sequence_startGap)) = '';

	substr($sequence, length($sequence) - length($sequence_stopGap), length($sequence_stopGap)) = '';
	substr($graph, length($graph) - length($sequence_stopGap), length($sequence_stopGap)) = '';
	
	if(length($sequence) == 0)
	{
		die "Problem with contig $graph_contig\n";
	}	
	
	die unless(length($sequence) == length($graph));
	
	my $disagreementNumber = 0;
	my $continuousOK = 0;
	my $inDisagreement = 0;
	my $disagreementStart = 0;
	my $disagreementStop = 0;
	my $totalOK = 0;
	
	my $nucleotideOK = 0;
	my $nucleotideTotal = 0;	
	my $nucleotideInIsolation = 0;
	
	for(my $i = 0; $i < length($sequence); $i++)
	{
		if(substr($sequence, $i, 1) ne '_')
		{
			$nucleotideTotal++;
		}
		
		if(($i > 0) and ($i < (length($sequence) - 1)))
		{
			if(substr($sequence, $i, 1) ne '_')
			{
				if((substr($sequence, $i-1, 1) eq '_') and (substr($sequence, $i+1, 1) eq '_'))
				{
					$nucleotideInIsolation++;
				}
			}
		}
		
		if(substr($sequence, $i, 1) eq substr($graph, $i, 1))
		{
			if(substr($sequence, $i, 1) ne '_')
			{
				$nucleotideOK++;
			}
			
			$continuousOK++;
			$totalOK++;
			if($inDisagreement and ($continuousOK > 20))
			{
				$disagreementStop = $i - $continuousOK;
				$inDisagreement = 0;
				
				
				my $flank_left_graph = substr($graph, $disagreementStart - 15, 15);
				my $flank_left_sequence = substr($sequence, $disagreementStart - 15, 15);

				my $flank_right_graph = substr($graph, $disagreementStop + 1, 15);
				my $flank_right_sequence = substr($sequence, $disagreementStop + 1, 15);				
				
				my $disagreement_graph = substr($graph, $disagreementStart, $disagreementStop - $disagreementStart + 1);
				my $disagreement_sequence = substr($sequence, $disagreementStart, $disagreementStop - $disagreementStart + 1);
				
				if(length($disagreement_graph) > 5)
				{
					print "\tDisagreement $disagreementNumber $disagreementStart - $disagreementStop \n";
				
					print "\t\tG: [", $flank_left_graph, "] ", $disagreement_graph, " [", $flank_right_graph, "]", "\n";
					print "\t\tS: [", $flank_left_sequence, "] ", $disagreement_sequence, " [", $flank_right_sequence, "]", "\n";
					print "\n";
				}
			}
		}
		else
		{
			if(not $inDisagreement)
			{
				$inDisagreement = 1;
				$disagreementStart = $i;
				$disagreementNumber++;
			}
			$continuousOK = 0;
		}
	}
	
	print "\tContig $contigI $totalOK / ", length($sequence), " positions identical -- $nucleotideOK / $nucleotideTotal of sequence characters, $nucleotideInIsolation in isolation!\n";
		
	# my @agreement_scores;
	
	# my $iPlus = 10;
	# for(my $i = 0; $i < length($sequence); $i += $iPlus)
	# {
		# my $ok = 0;
		# for(my $j = $i; $j < $iPlus; $j++)
		# {
			# if(substr($sequence, $j, 1) eq substr($graph, $sequence, 1))
			# {
				# $ok++;
			# }
		# }
		# push(@agreement_scores, $ok/$iPlus);
	# }
}

