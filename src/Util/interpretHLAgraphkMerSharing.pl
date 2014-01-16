#!/usr/bin/perl
use Modern::Perl;
use List::MoreUtils qw/mesh all/;
use Data::Dumper;

my $kMer_sharing_file = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/graph_kMersPerLevel_HLAGRAPH_31.txt);
my $output_dir = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/kMerSharing);
mkdir($output_dir) unless (-e $output_dir);

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
	
	open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
	print OUTPUT join("\t", qw/Level Part kMers_unweighted kMers_unweighted_inOtherLevels kMers_unweighted_inReference kMers_weighted kMers_weighted_inOtherLevels kMers_weighted_inReference kMer/), "\n";
	
	for(my $locusLevelI = 0; $locusLevelI <= $#{$data_per_locus{$locus}}; $locusLevelI++)
	{
		my $graphLevel = $data_per_locus{$locus}[$locusLevelI][0];
		my $part = $data_per_locus{$locus}[$locusLevelI][1];
		
		die unless(defined $graphLevel);
		
		my $kMers_unique = 0;
		my $kMers_unique_inOtherLevels = 0;
		my $kMers_unique_inReference = 0;

		my $kMers_weighted = 0;
		my $kMers_weighted_inOtherLevels = 0;
		my $kMers_weighted_inReference = 0;		
		
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
		}
				
		print OUTPUT join("\t",
			$graphLevel,
			$part,
			$kMers_unique,
			$kMers_unique_inOtherLevels,
			$kMers_unique_inReference,
			$kMers_weighted,
			$kMers_weighted_inOtherLevels,
			$kMers_weighted_inReference
		), "\n";		
	}
	
	close(OUTPUT);
}




sub reverseComplement
{
	my $in = shift;
	$in =~ tr/ACGTacgt/TGCAtgca/;
	return scalar(reverse $in);
}
