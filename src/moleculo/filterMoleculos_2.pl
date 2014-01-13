#!/usr/bin/perl
use strict;
use Data::Dumper;
use File::Copy;
use Getopt::Long;
use List::MoreUtils qw/mesh/;

my %contigs_for_extraction;
my %contigs_kMer_length;

$| = 1;

my $kMer_size = 31;
my $contigs_file = '/Net/banyan/data1/projects/gsk/moleculo/contigs.fasta';
my $contigs_filtered_file = $contigs_file.'.filtered';
my $contigs_stats_file = $contigs_file.'.stats';

my $contigs_extracted_file = '/Net/banyan/data1/projects/gsk/moleculo/contigs_xMHC.fasta';
my $contigs_stats_extracted_file = '/Net/banyan/data1/projects/gsk/moleculo/contigs_xMHC.fasta.stats';

open(CONTIGS_EXTRACTED, '>', $contigs_extracted_file) or die "Cannot open $contigs_extracted_file";
open(STATS_EXTRACTED, '>', $contigs_stats_extracted_file) or die "Cannot open $contigs_stats_extracted_file";

print "Reading file $contigs_stats_file ... \n";

open(STATS, '<', $contigs_stats_file) or die "Cannot open $contigs_stats_file";
my $headerLine = <STATS>;
chomp($headerLine);
my @header_fields = split(/\t/, $headerLine);
print STATS_EXTRACTED $headerLine, "\n";

my $desired_extraction_length = 0;
while(<STATS>)
{
	if(($. % 100000) == 0)
	{
		print "\t", $., "\n";
	}
	
	my $line = $_;
	chomp($line);
	my @line_fields = split(/\t/, $line);
	my %line_hash = (mesh @header_fields, @line_fields);
	
	die unless((exists $line_hash{ID}) and (exists $line_hash{fraction_xMHC}) and (exists $line_hash{uniqueStretch_length}) and (exists $line_hash{uniqueStretch_fractionUnique}) and (exists $line_hash{uniqueStretch_fractionElsewhere}));

	if((not ($line_hash{ID} =~ /^TEST/)) and ($line_hash{fraction_xMHC} > 0.8) and ($line_hash{uniqueStretch_length} > 50) and ($line_hash{uniqueStretch_fractionUnique} > 0.5) and ($line_hash{uniqueStretch_fractionElsewhere} <= 0.3))
	{
		die unless((exists $line_hash{firstUnique}) and (exists $line_hash{lastUnique}));
		die if(($line_hash{firstUnique} eq '-1') and ($line_hash{lastUnique} eq '-1'));
		
		my $firstSeqPos = $line_hash{firstUnique};
		my $lastSeqPos = $line_hash{lastUnique} + $kMer_size - 1;		
		$contigs_for_extraction{$line_hash{ID}} = [$firstSeqPos, $lastSeqPos];
		$contigs_kMer_length{$line_hash{ID}} = $line_hash{uniqueStretch_length};
		
		print STATS_EXTRACTED $line, "\n";
		
		$desired_extraction_length += ($lastSeqPos - $firstSeqPos + 1);
	}
}
close(STATS);
close(STATS_EXTRACTED);

print "\nWill now extract ", scalar(keys %contigs_for_extraction), " contigs; total length: ", $desired_extraction_length, ", avg. length: ", $desired_extraction_length/scalar(keys %contigs_for_extraction), "\n";

print "Reading file $contigs_file ... \n";

my $sum_sequence_length = 0;
my %found_contigs;
open(CONTIGS, '<', $contigs_file) or die "Cannot open $contigs_file";
while(<CONTIGS>)
{
	if(($. % 100000) == 0)
	{
		print "\t", $., "\n";
	}
	
	my $IDline = $_;
	chomp($IDline);
	die unless(substr($IDline, 0, 1) eq '>');	
	my $ID = substr($IDline, 1);
	
	my $nextLine = <CONTIGS>;
	chomp($nextLine);
	
	if($contigs_for_extraction{$ID})
	{
		my $sequence_for_extraction = substr($nextLine, $contigs_for_extraction{$ID}[0], $contigs_for_extraction{$ID}[1] - $contigs_for_extraction{$ID}[0] + 1);			
		$sum_sequence_length += length($sequence_for_extraction);
		print CONTIGS_EXTRACTED $IDline, "\n", $sequence_for_extraction, "\n";		
		$found_contigs{$ID} = 1; 
		
		die unless($contigs_kMer_length{$ID});
		my $sequence_extraction_kMer_length = length($sequence_for_extraction) - $kMer_size + 1;
		die unless($contigs_kMer_length{$ID} == $sequence_extraction_kMer_length);
	}
}
close(CONTIGS);

print "\n";
foreach my $contigID (keys %contigs_for_extraction)
{
	unless($found_contigs{$contigID})
	{
		warn "Contig $contigID coult not be found!\n";
	}
}

print "Extracted ", scalar(keys %found_contigs), " contigs; total length: ", $sum_sequence_length, ", avg. length: ", $sum_sequence_length/scalar(keys %found_contigs), "\n";


