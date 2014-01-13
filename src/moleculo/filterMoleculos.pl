#!/usr/bin/perl
use strict;
use Data::Dumper;
use File::Copy;
use Getopt::Long;

$| = 1;
# die Dumper(&string2Kmers(4, "ACGTAC"));
die Dumper("reverse complement problem!", "ACGTacgtNn", reverseComplement("ACGTacgtNn"), reverseComplement(reverseComplement("ACGTacgtNn"))) unless(reverseComplement(reverseComplement("ACGTacgtNn")) eq "ACGTacgtNn");

my $thisThread = 0;
my $totalThreads = 1;
my $mode;

GetOptions(
	'tI:s' => \$thisThread,
	'T:s' => \$totalThreads,
	'm:s' => \$mode,
);

my $contigs_file = '/data1/projects/gsk/moleculo/contigs.fasta';
my $contigs_filtered_file = $contigs_file.'.filtered';
my $contigs_stats_file = $contigs_file.'.stats';

my $path_to_PGF_haplotype = qq(/Net/birch/data/dilthey/forDownload/mhc_ref_8_haplotypes/pgf_ref.fasta);

	
if($mode eq '')
{
	for(my $tI = 0; $tI < $totalThreads; $tI++)
	{
		print &cmd_for_tI($tI), "\n";
	}
}
elsif($mode eq 'execute')
{
	my $kMer_size = 31;
	my $graph = '../../tmp2/GS_nextGen/varigraph3';

		
	my $expected_normal_graph = $graph.'/graph.txt';
	die "Normal graph $expected_normal_graph not there!" unless(-e $expected_normal_graph);

	my $expected_kmer_graph_file = $expected_normal_graph.".kmers_".$kMer_size;
	die "Augmented kMer graph $expected_kmer_graph_file not there!" unless(-e $expected_kmer_graph_file);

	my $path_to_HumanGenome_graph = qq(/Net/birch/data/dilthey/forDownload/mhc_ref_8_haplotypes/alex_grc37_auto_X_Y.k${kMer_size}.ctx);
	my $interesting_kMers_file = $expected_kmer_graph_file.'.requiredKMers';

	die "Cannot access $path_to_HumanGenome_graph" unless (-e $path_to_HumanGenome_graph);
	die "Cannot access $interesting_kMers_file" unless (-e $interesting_kMers_file);

	my $expected_normal_graph = $graph.'/graph.txt';
	die "Normal graph $expected_normal_graph not there!" unless(-e $expected_normal_graph);

	print "Reading $interesting_kMers_file ... \n";
	my %interesting_kMers;
	open(KMERS, '<', $interesting_kMers_file) or die "Cannot open $interesting_kMers_file";
	while(<KMERS>)
	{
		my $line = $_;
		chomp($line);
		$interesting_kMers{$line}++;
	}
	close(KMERS);

	my $kMers_reference_counts_file = $interesting_kMers_file.'.binaryCount';
	my $kMers_reference_counts_corrected = 'tmp/referenceCounts_corrected.txt';

	unless(-e $kMers_reference_counts_file)
	{
		my $cmd = qq(../readCortexCoverage.pl --cortex_bin $path_to_HumanGenome_graph  --kMer_size $kMer_size --interesting_kMers $interesting_kMers_file);
		my $output = `$cmd`;
		unless(-e $kMers_reference_counts_file)
		{
			die "Could not count reference kMers\n\nCommand:\n $cmd \n\nOutput:\n\n$output";
		}	
	}

	print "Reading $kMers_reference_counts_file ... \n";
	my %kMers_reference_counts;
	open(REFERENCECOUNTS, '<', $kMers_reference_counts_file) or die "Cannot open $kMers_reference_counts_file";
	while(<REFERENCECOUNTS>)
	{
		my $line = $_;
		chomp($line);
		die Dumper("RegExp", $kMers_reference_counts_file, $., $line) unless($line =~ /^(\w+)\s+(\d+)$/);
		my $kMer = $1;
		my $count = $2;
		$kMers_reference_counts{$kMer} = $count;
	}
	close(REFERENCECOUNTS);

	foreach my $kMer (keys %interesting_kMers)
	{
		unless(exists $kMers_reference_counts{$kMer})
		{
			die Dumper("Missing kMer", $kMer);
		}	
	}		

	my %kMers_reference_counts_corrected;
	unless(-e $kMers_reference_counts_corrected)
	{
		my %loci_in_graph;
		open(GRAPH, "<", $expected_normal_graph) or die "Cannot open $expected_normal_graph\n";
		my $firstLine = <GRAPH>;
		chomp($firstLine);
		unless($firstLine eq 'CODE:')
		{
			die "Expect first line to be CODE: - graph format changed, file $expected_normal_graph?";
		}
		while(<GRAPH>)
		{
			my $line = $_;
			chomp($line);
			if($line =~ /\:$/)
			{
				last;
			}
			my @fields = split(/\|\|\|/, $line);
			die unless($#fields == 2);
			my $locusID = $fields[0];
			$loci_in_graph{$locusID} = 1;
		}
		close(GRAPH);
		
		my @positions_in_graph = grep {$_ =~ /^S\d_(\d+)_(\d+)$/} keys %loci_in_graph;
		@positions_in_graph = map {$_ =~ /^S\d_(\d+)_(\d+)$/; $2} @positions_in_graph;
		
		unless(scalar(@positions_in_graph) > 1)
		{
			die "Have less than two loci with integrated coordinates - chr6:1234... - in the graph. Cannot determine graph boundaries";
		}
		
		# find the positions of the first and the last SNP in the graph on reference coordinates
		
		@positions_in_graph = sort {$a <=> $b} @positions_in_graph;
		
		my $first_position_graph_inPGF = $positions_in_graph[0];
		my $last_position_graph_inPGF = $positions_in_graph[$#positions_in_graph];
		
		die unless($first_position_graph_inPGF < $last_position_graph_inPGF);
		
		my %subtract_kMers;
		print "Remove kMers for positions $first_position_graph_inPGF - $last_position_graph_inPGF (in PGF)\n";
		
		# read the pgf (reference) haplotype, extract the stretch spanned by first and last SNP in graph,
		# count kMers, subtract from reference genome kMer count
		
		my $pgf_haplotype = read_PGF();
		(length($pgf_haplotype) > $first_position_graph_inPGF+$kMer_size-1) or die Dumper("PGF coordinate problem!", length($pgf_haplotype), $first_position_graph_inPGF, $kMer_size)."\n\n";
		
		my $kMer = substr($pgf_haplotype, $first_position_graph_inPGF, $kMer_size-1);
		(length($kMer) == ($kMer_size-1)) or die "Problem with extracting position $first_position_graph_inPGF + $kMer_size symbols, first position in graph $first_position_graph_inPGF , last position in graph $last_position_graph_inPGF , total PGF length ".length($pgf_haplotype);
		
		for(my $end_KMer = $first_position_graph_inPGF+$kMer_size-1; $end_KMer <= $last_position_graph_inPGF; $end_KMer++)
		{
			my $char_PGF = substr($pgf_haplotype, $end_KMer, 1);
			(length($char_PGF) == 1) or die;
			$kMer = $kMer.$char_PGF;
			(length($kMer) == $kMer_size) or die;
			$subtract_kMers{$kMer}++;
			substr($kMer, 0, 1) = '';       
		}
		
		# die Dumper(\%subtract_kMers);
		
		open(OUTPUT, '>', $kMers_reference_counts_corrected) or die "Cannot open $kMers_reference_counts_corrected";
		my $have_subtracted = 0;
		foreach my $kMer (sort keys %kMers_reference_counts)
		{
			my $original_count = $kMers_reference_counts{$kMer};
			my $subtract_count = (exists $subtract_kMers{$kMer}) ? $subtract_kMers{$kMer} : 0;
			my $new_count = $original_count - $subtract_count;
			if($subtract_count > 0)
			{
				$have_subtracted = 1;
			}
			die unless($new_count >= 0);
			print OUTPUT join(' ', $kMer, $new_count), "\n";
		}
		close(OUTPUT);   
		
		unless($have_subtracted)
		{
			unlink $kMers_reference_counts_corrected;		
			die "No PGF kMer was subtracted to correct the reference counts. Bad. Fix, and --redo!\n";
		}	
	}
	print "Reading $kMers_reference_counts_corrected ... \n";
	open(REFERENCECOUNTS_CORRECTED, '<', $kMers_reference_counts_corrected) or die "Cannot open $kMers_reference_counts_corrected";
	while(<REFERENCECOUNTS_CORRECTED>)
	{
		my $line = $_;
		chomp($line);
		die Dumper("RegExp (kMers corrected)", $kMers_reference_counts_corrected, $., $line) unless($line =~ /^(\w+)\s+(\d+)$/);
		my $kMer = $1;
		my $count = $2;
		$kMers_reference_counts_corrected{$kMer} = $count;
	}
	close(REFERENCECOUNTS_CORRECTED);

	foreach my $kMer (keys %interesting_kMers)    
	{
		unless(exists $kMers_reference_counts_corrected{$kMer})
		{
			die Dumper("Missing kMer (corrected)", $kMer);
		}	
	}		


	# now process contigs

	print "\nProcessing $contigs_file ... \n";   

	my $this_t_contigs_stats_file = $contigs_stats_file .'_thread_'.$thisThread.'_of_'.$totalThreads;
	my $this_t_contigs_filtered_file = $contigs_filtered_file .'_thread_'.$thisThread.'_of_'.$totalThreads;

	my $status_contigs_stats_file = $this_t_contigs_stats_file.'.status';
	open(STATUS, '>', $status_contigs_stats_file) or die "Cannot open status file $status_contigs_stats_file";
	print STATUS 0;
	close(STATUS);
	
	open(STATS, '>', $this_t_contigs_stats_file) or die "Cannot open stats file: $this_t_contigs_stats_file";
	open(FILTERED, '>', $this_t_contigs_filtered_file) or die "Cannot open $this_t_contigs_filtered_file";
	print STATS join("\t", qw/ID kMers fraction_xMHC fraction_unique_xMHC firstUnique lastUnique uniqueStretch_length uniqueStretch_fractionElsewhere uniqueStretch_fractionUnique/), "\n";

	open(CONTIGS, '<', $contigs_file) or die "Cannot open $contigs_file";
	my $contigI = 0;
	
	while(<CONTIGS>)
	{
		if(($. % 1000) == 0)
		{
			# print "\t$.\n";
		}
		
		my $firstLine = $_;
		my $secondLine = <CONTIGS>;
		
		$contigI++;
		
		if(($contigI % $totalThreads) != $thisThread)
		{	
			next;
		}
		
		chomp($firstLine);
		chomp($secondLine);
		
		die unless(substr($firstLine, 0, 1) eq '>');
		my $ID = substr($firstLine, 1);
		
		my $fraction_xMHC = 0;
		my $fraction_unique_xMHC = 0;
		my $firstUnique = -1;
		my $lastUnique = -1;
		my $uniqueStretch_length = 0;
		my $uniqueStretch_fractionElsewhere = 0;
		my $uniqueStretch_fractionUnique = 0;
		
		my @kMers = &string2Kmers($kMer_size, $secondLine);
		
		if(scalar(@kMers) > 0)
		{
			my $kMers_xMHC = 0;
			my $kMers_unique_xMHC = 0;
			my $kMers_xMHC_elseWhere = 0;
			
			my $first_unique_kMer;
			my $last_unique_kMer;
			
			for(my $i = 0; $i <= $#kMers; $i++)
			{
				my $kMer = $kMers[$i];			
				my $kMer_rC = reverseComplement($kMer);
				
				if($interesting_kMers{$kMer} or $interesting_kMers{$kMer_rC})
				{
					$kMers_xMHC++;
					
					if(($kMers_reference_counts_corrected{$kMer} == 0) and ($kMers_reference_counts_corrected{$kMer_rC} == 0))
					{
						$kMers_unique_xMHC++;
						
						if(not defined $first_unique_kMer)
						{
							$first_unique_kMer = $i;
						}
						
						$last_unique_kMer = $i;
					}
				}
			}
				
			$fraction_xMHC = $kMers_xMHC / scalar(@kMers);
			$fraction_unique_xMHC = $kMers_unique_xMHC / scalar(@kMers);
			if(defined $first_unique_kMer)
			{
				$firstUnique = $first_unique_kMer;
				$lastUnique = $last_unique_kMer;
			}
			$uniqueStretch_length = ($kMers_unique_xMHC > 0) ? ($last_unique_kMer - $first_unique_kMer + 1) : 0;
			$uniqueStretch_fractionUnique = ($uniqueStretch_length > 0) ? ($kMers_unique_xMHC / $uniqueStretch_length) : 0;
			
			if($first_unique_kMer)
			{
				my $firstSeqPos = $first_unique_kMer;
				my $lastSeqPos = $last_unique_kMer + $kMer_size - 1;
						
				my $potentiallyUniqueSequence = substr($secondLine, $firstSeqPos, $lastSeqPos - $firstSeqPos + 1);
				my @kMers_covered = &string2Kmers($kMer_size, $potentiallyUniqueSequence);
				die unless(scalar(@kMers_covered) == $uniqueStretch_length);
				for(my $i = 0; $i <= $#kMers_covered; $i++)
				{
					my $kMer = $kMers_covered[$i];			
					my $kMer_rC = reverseComplement($kMer);
					
					if($kMers_reference_counts_corrected{$kMer} or $kMers_reference_counts_corrected{$kMer_rC})
					{
						$kMers_xMHC_elseWhere++;
					}
				}
						
				$uniqueStretch_fractionElsewhere = ($uniqueStretch_length > 0) ? ($kMers_xMHC_elseWhere / $uniqueStretch_length) : 0;
				
				
				if((not ($ID =~ /^TEST/)) and ($uniqueStretch_length > 50) and ($uniqueStretch_fractionUnique > 0.5) and ($uniqueStretch_fractionElsewhere <= 0.3))
				{
					print FILTERED '>', $ID, "\n";
					print FILTERED $potentiallyUniqueSequence, "\n";
				}
				
			}
		}
		
		print STATS join("\t", $ID, scalar(@kMers), $fraction_xMHC, $fraction_unique_xMHC, $firstUnique, $lastUnique, $uniqueStretch_length, $uniqueStretch_fractionElsewhere, $uniqueStretch_fractionUnique), "\n";
		
	}
	close(CONTIGS);

	close(STATS);
	close(FILTERED);
	
	open(STATUS, '>', $status_contigs_stats_file) or die "Cannot open status file $status_contigs_stats_file";
	print STATUS 1;
	close(STATUS);	
}
elsif($mode eq 'collect')
{
	my $allOK = 1;
	for(my $tI = 0; $tI < $totalThreads; $tI++)
	{
		my $this_t_contigs_stats_file = $contigs_stats_file .'_thread_'.$tI.'_of_'.$totalThreads;
		my $this_t_contigs_filtered_file = $contigs_filtered_file .'_thread_'.$tI.'_of_'.$totalThreads;

		my $status_contigs_stats_file = $this_t_contigs_stats_file.'.status';
		open(STATUS, '<', $status_contigs_stats_file) or die "Cannot open status file $status_contigs_stats_file";
		my $status = <STATUS>;
		close(STATUS);
		
		unless(substr($status, 0, 1) eq '1')
		{
			warn "Thread $tI not completed yet or failed - re-execute? ".&cmd_for_tI($tI);
			$allOK = 0;
		}
	}
	
	unless($allOK)
	{
		exit;
	}
	
	open(STATS, '>', $contigs_stats_file) or die "Cannot open stats file: $contigs_stats_file";
	open(FILTERED, '>', $contigs_filtered_file) or die "Cannot open $contigs_filtered_file";
	
	for(my $tI = 0; $tI < $totalThreads; $tI++)
	{
		my $this_t_contigs_stats_file = $contigs_stats_file .'_thread_'.$tI.'_of_'.$totalThreads;
		my $this_t_contigs_filtered_file = $contigs_filtered_file .'_thread_'.$tI.'_of_'.$totalThreads;

		open(IN_STATS, '<', $this_t_contigs_stats_file) or die "Cannot open $this_t_contigs_stats_file";
		my $header = <IN_STATS>;
		if($tI == 0)
		{
			print STATS $header;
		}
		while(<IN_STATS>)
		{
			print STATS $_;
		}
		close(IN_STATS);
		
		open(IN_CONTIGS, '<', $this_t_contigs_filtered_file) or die "Cannot open $this_t_contigs_filtered_file";
		while(<IN_CONTIGS>)
		{
			print FILTERED $_;
		}		
		close(IN_CONTIGS);
	}
	
	
	close(STATS);
	close(FILTERED);
	
	
	
	
}
else
{
	die "Unknown --m: $mode";
}


sub read_PGF
{
	my %alignment;
	my $file = $path_to_PGF_haplotype;
	open(ALIGNMENT, "<", $file) or die "Cannot open $file";
	my $current_name;
	while(<ALIGNMENT>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		if($line =~ /^>/)
		{
			die if ($current_name);
			$current_name = 'pgf';
		}
		else
		{
			die unless($current_name);
			$alignment{$current_name} .= $line;
		}
	}
	close(ALIGNMENT);
	$alignment{'pgf'} or die;
	return $alignment{'pgf'};
}

sub string2Kmers
{
	my $k = shift;
	my $string = shift;
	
	if(length($string) < $k)
	{	
		return ();
	}
	
	# $i is the last character of the kMer
	my @kMers;
	my $runningkMer;
	for(my $i = ($k - 1); $i <= (length($string) - 1); $i++)
	{
		if($i == ($k - 1))
		{
			my $firstMer = substr($string, 0, $k);
			die unless(length($firstMer) == $k);
			push(@kMers, $firstMer);
			$runningkMer = $firstMer;
		}
		else
		{
			$runningkMer = substr($runningkMer, 1);   
			$runningkMer .= substr($string, $i, 1);
			die unless(length($runningkMer) == $k);
			push(@kMers, $runningkMer);
		}
	}
	
	die Dumper("Wrong kMer count", \@kMers, $string) unless(scalar(@kMers) == (length($string) - $k + 1));
	
	return @kMers;
}

sub cmd_for_tI
{
	my $tI = shift;
	my $cmd = qq(./filterMoleculos.pl --tI $tI --T $totalThreads --m execute &);
	return $cmd;
}

sub reverseComplement
{
	my $input = shift;
	$input =~ tr/ACGTacgt/TGCAtgca/;
	$input = reverse $input;
	return $input;
}	




