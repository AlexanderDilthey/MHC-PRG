#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;


my %samples_to_populations;
open(POPULATIONS, '<', 'forPaper/_1000G_populations.txt') or die;
while(<POPULATIONS>)
{
	my $line = $_;
	chomp($line);
	my @fields = split(/\t/, $line);
	$samples_to_populations{$fields[0]} = $fields[1];
}

my %populations;
my %types_2digit;

my $inputFile = '_Platinum_types_as_validated.txt';
open(F, '<', $inputFile) or die "Cannot open $inputFile";
my $headerLine = <F>;
chomp($headerLine);
my @headerFields = split(/\t/, $headerLine);
my @loci = grep {$_ =~ /HLA/} @headerFields;
while(<F>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	die Dumper("Field mismatch", $#line_fields, $#headerFields, $.) unless($#line_fields == $#headerFields);
	
	my %line = (mesh @headerFields, @line_fields);
	
	my $indivID = $line{'IndividualID'};
	die unless(defined $indivID);
	
	die unless(exists $samples_to_populations{$indivID});
	
	my $population = $samples_to_populations{$indivID};
	
	$populations{$population}++;
	
	foreach my $locus (@loci)
	{
		my $alleles = $line{$locus};
		die unless($alleles);
		my @alleles = split(/\//, $alleles);
		die unless($#alleles == 1);
		
		foreach my $allele (@alleles)
		{
			if($allele =~ /\?/)
			{
				$allele = 'Missing';
			}
			else
			{
				my @components = split(/:/, $allele);
				die unless(scalar(@components) >= 2);
				$allele = $components[0];
			}
			$types_2digit{$locus}{$allele}++;
		}
	}
}
close(F);

my $fn_output_populations = $inputFile . '.populations';
open(OUT, '>', $fn_output_populations) or die "Cannot open $fn_output_populations";
print OUT join("\t", "Population", "ValidationSamples"), "\n";
foreach my $population (sort keys %populations)
{
	print OUT join("\t", $population, $populations{$population}), "\n";
}
close(OUT);

my $fn_output_twodigit = $inputFile . '.groups_2digit';
open(OUT, '>', $fn_output_twodigit) or die "Cannot open $fn_output_twodigit";
print OUT join("\t", "Locus", "Group2Digit", "N"), "\n";
foreach my $locus (sort keys %types_2digit)
{
	foreach my $twoDigGroup (sort keys %{$types_2digit{$locus}})
	{
		print OUT join("\t", $locus, $twoDigGroup, $types_2digit{$locus}{$twoDigGroup}), "\n";
	}
}
close(OUT);
