#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/all/;

my @coverage_targets = (40, 30, 20);
my $iterations = 3;
my %samples = (
	'I1_AA02O9Q_Z2' => 53,
#	'I3_NA12878' => 20,
);

foreach my $sample (keys %samples)
{
	my $current_coverage = $samples{$sample};
	die "Coverage problem" unless(all {$_ < $current_coverage} @coverage_targets);
	
	
}

