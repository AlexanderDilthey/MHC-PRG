#!/usr/bin/perl
use Modern::Perl;
my $inputFile = qq(/Net/birch/data/dilthey/MHC-PRG/src/varigraph3_31_compare.txt);

my @categories_chart_1 = qw/ALL ALL CleanSNP CleanSNP INDEL INDEL/;
my @categories_chart_2 = qw/realIndel_L1 realIndel_L5 realIndel_L10-19 realIndel_L100-149 realIndel_L1000-1999 realIndel_L10000-10999 realIndel_L60000-60999  realIndel_L124000-124999/;

my %extraction;
my $inContent = 0;
open(F, '<', $inputFile) or die "Cannot open $inputFile";
while(<F>)
{
	my $line = $_;
	chomp($line);
	if($line =~ /Total summary/)
	{
		$inContent = 1;
	}
	elsif($inContent)
	{
		die "Cannot parse line: $line" unless($line =~ /^\s+([A-Za-z_\d\-]+)\s*$/);
		my $cat = $1;
		my @l;
		for(my $i = 1; $i <= 9; $i++)
		{
			push(@l, scalar(<F>));
		}
		
		my $alleles_OK = 0;
		my $alleles_notOK = 0;
		
		for(my $alleleOK = 0; $alleleOK <= 2; $alleleOK++)
		{
			my $alleles_nonOK = 2 - $alleleOK;
			die unless($l[$alleleOK] =~ /(agreement_\d) (\d+)\//);
			die unless($1 eq 'agreement_'.$alleleOK);
			my $alleles = $2;
			
			$alleles_OK += ($alleles * $alleleOK);
			$alleles_notOK += ($alleles * $alleles_nonOK);
		}
		
		my $rate_OK = $alleles_OK/($alleles_OK + $alleles_notOK);
		print $cat, ": ", $rate_OK, "\n";
		
		$extraction{$cat} = [$rate_OK, ($alleles_OK + $alleles_notOK)];
	}
}
close(F);

# chart 1
open(C1, '>', 'data/forR_simulationResults_chart_1.tab') or die;
print C1 join("\t", qw/Category RateOK N RateOK_HET N_HET/), "\n";
foreach my $catName (@categories_chart_1)
{
	my $catName_HET = $catName.'_HET';
	die unless($extraction{$catName});
	die unless($extraction{$catName_HET});
	
	print C1 join("\t", $catName, $extraction{$catName}[0], $extraction{$catName}[1], $extraction{$catName_HET}[0], $extraction{$catName_HET}[1]), "\n";
}
close(C1);

# chart 2
open(C2, '>', 'data/forR_simulationResults_chart_2.tab') or die;
print C2 join("\t", qw/Category RateOK N RateOK_HET N_HET/), "\n";
foreach my $catName (@categories_chart_2)
{
	my $catName_HET = $catName.'_HET';
	die unless($extraction{$catName});
	die unless($extraction{$catName_HET});
	
	print C2 join("\t", $catName, $extraction{$catName}[0], $extraction{$catName}[1], $extraction{$catName_HET}[0], $extraction{$catName_HET}[1]), "\n";
}
	
close(C2);


