#!/usr/bin/perl
use Modern::Perl;
$| = 1;

my @kMers = qw/AATTTTCTCCCATTTTGTAGGTTGC AGGTTGCGAAAATTTTCTCCCATTT GTTGCGAAAATTTTCTCCCATTTTG TAGGTTGCGAAAATTTTCTCCCATT TCTCCCATTTTGTAGGTTGCCTGTT TCTTGTAAATTTGTTTGAGTTCATT TGTAGGTTGCCTGTTCACTCTGATG TTCTCCCATTTTGTAGGTTGCCTGT TTGCTGTGCAGAAGCTCTTTAGTTT TTTCTCCCATTTTGTAGGTTGCCTG/;
my $input = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/derived/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph);

my $reference_href = readFASTA($input);
my $reference_href_complement = complementFASTA($reference_href);


foreach my $kMer (@kMers)
{
	print "Searching $kMer\n";
	print "\tIn reference       : ", join(', ', isIn($reference_href, $kMer)), "\n";
	print "\tIn reference compl.: ", join(', ', isIn($reference_href_complement, $kMer)), "\n";
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

