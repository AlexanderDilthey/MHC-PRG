#!/usr/bin/perl
use Modern::Perl;
$| = 1;

my $k = 25;

my $output_file = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/classI_exon_kMers.txt);

my @loci = qw/A B C/;
my @exons = (2, 3);
my %kMers;

foreach my $locus (@loci)
{
	foreach my $exon (@exons)
	{
		my $fn_search = "/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/${locus}_*_exon_${exon}.txt";
		my @f = glob($fn_search);
		die unless($#f == 0);
		my $fn = $f[0];
		open(F, '<', $fn) or die "Cannot open $fn";
		<F>; 
		while(<F>)
		{
			my $l = $_;
			chomp($l);
			$l =~ s/\n//g;
			$l =~ s/\r//g;
			$l =~ s/\*/N/g;
			next unless($l);
			my @f = split(/ /, $l);
			shift(@f);
			my $seq = join('', @f);
			die "Weird sequuence: $seq" unless($seq =~ /^[ACGT_N]+$/);
			
			my $seq_noGaps = $seq;
			$seq_noGaps =~ s/_//g;
			my @kMers = kmers($seq_noGaps, $k);
			foreach my $kMer (@kMers)
			{
				$kMers{$kMer}++;
			}
		}
		close(F);
	}
}

open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
foreach my $kMer (sort {$kMers{$b} <=> $kMers{$a}} keys %kMers)
{
	print OUTPUT join("\t", $kMer, $kMers{$kMer}), "\n";
}
close(OUTPUT);
print "\nWritten $output_file\n";

sub kmers
{
	my $s = shift;
	my $k = shift;
	
	die unless($s =~ /^[ACGTN]+$/);
	if(length($s) < $k)
	{
		return ();
	}
	my @forReturn;
	my $expectedMers = length($s) - $k + 1;
	for(my $i = 0; $i < $expectedMers; $i++)
	{	
		my $kMer = substr($s, $i, $k);
		die unless(length($kMer) == $k);
		if($kMer !~ /N/)
		{
			push(@forReturn, $kMer);
		}
	}
	return @forReturn;
}