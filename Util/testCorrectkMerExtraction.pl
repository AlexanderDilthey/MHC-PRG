#!/usr/bin/perl
use Modern::Perl;
$| = 1;

my $k = 25;

my %sequences = (
'>Ex2_0201' => 
	'GCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG',
'>Ex2_2301' => 
	'GCTCCCACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGAGGAGACAGGGAAAGTGAAGGCCCACTCACAGACTGACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCG',

'>Ex3_0201' => 
	'GTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGG',
'>Ex3_2301' => 
	'GTTCTCACACCCTCCAGATGATGTTTGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGG',
)
;
my $input = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph);

my $reference_href = readFASTA($input);
my $reference_href_complement = complementFASTA($reference_href);


foreach my $seqID (keys %sequences)
{
	print $seqID, "\n";
	my @kMers = kmers($sequences{$seqID}, $k);
	my $kMer_in = 0;
	my $kMer_notIn = 0;
	foreach my $kMer (@kMers)
	{
		print "\t\t", $kMer, " ";
		if(isIn($reference_href, $kMer) or isIn($reference_href_complement, $kMer))
		{
			$kMer_in++;
			print "1";
		}
		else
		{
			$kMer_notIn++;
			print "0";
		}
		print "\n";
	}
	
	my $perc_unique = sprintf("%.2f", ($kMer_notIn / ($kMer_in + $kMer_notIn)) * 100).'%';
	
	print "\t Unique: ", $perc_unique, "\n";
	print "\t\t", $kMer_in, "\t", $kMer_notIn, "\n";	
}

sub isIn
{
	my $inWhat_href = shift;
	my $kMer = shift;
	
	my @results_inWhat;
	
	foreach my $key (keys %$inWhat_href)
	{
		if(index($inWhat_href->{$key}, $kMer) != -1)
		{
			push(@results_inWhat, $key);
		}
	}
	
	return @results_inWhat;
}


sub complementFASTA
{
	my $in_href = shift;
	my %forReturn;
	
	foreach my $key (keys %$in_href)
	{
		my $revComp = reverseComplement($in_href->{$key});
		$forReturn{$key} = $revComp;
	}
	return \%forReturn;
}


sub reverseComplement
{
	my $input = shift;
	$input =~ tr/ACGTacgt/TGCAtgca/;
	$input = reverse $input;
	return $input;
}	

sub readFASTA
{
	my $file = shift;
	
	print "Reading $file\n";
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000) == 0)
		{
			print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		
		if(substr($line, 0, 1) eq '>')
		{
			$currentSequence = substr($line, 1);
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
	
	print "\n\n";
	
	return \%R;
}


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