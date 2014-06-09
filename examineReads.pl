#!/usr/bin/perl

my %reads;
open(GOOD_READS, '<', 'test_round2.fastQ_1') or die;
while(<GOOD_READS>)
{
	my $line = $_;
	chomp($line);
	die unless(substr($line, 0, 2) eq '@@');
	die unless(substr($line, length($line) - 2, 1) eq '/');
	substr($line, length($line) - 2, 2) = '';
	my $readID = $line;
	<GOOD_READS>;
	<GOOD_READS>;
	<GOOD_READS>;
	
	$reads{$readID} = 1;
}		

print "Good reads: ", scalar(keys %reads), "\n";

open(GOODREADIDS, '>', 'goodreadIDs.txt') or die;
print GOODREADIDS join("\n", map {substr($_, 0, 1) = ''; $_} keys %reads);
close(GOODREADIDS);

my %other_reads;
open(OTHER_READS, '<', 'round3_positive_25_separateOptim.fastQ_1') or die;
while(<OTHER_READS>)
{
	my $line = $_;
	chomp($line);
	$line = '@'.$line;
	die unless(substr($line, 0, 2) eq '@@');
	die unless(substr($line, length($line) - 2, 1) eq '/');
	substr($line, length($line) - 2, 2) = '';
	my $readID = $line;
	<OTHER_READS>;
	<OTHER_READS>;
	<OTHER_READS>;
	
	$other_reads{$readID} = 1;
}		

print "Good reads in other file: ", scalar(grep {$other_reads{$_}} keys %reads), "\n";
