#!/usr/bin/perl
use Modern::Perl;
$| = 1;

my %block_positions_chr6;

my $genomicMapping = qq(/Net/birch/data/dilthey/MHC-PRG/tmp2/GS_nextGen/hla/genomicMapping.txt);
my $inputGenome = qq(/gpfs1/well/gsk_hla/GRCh37.60/fasta/combined/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.fa);
my $outputGenome = qq(Homo_sapiens.GRCh37.60.dna.chromosome.ALL.blockedHLAgraph);

# $inputGenome = qq(/gpfs1/well/gsk_hla/GRCh37.60/fasta/Homo_sapiens.GRCh37.60.dna.chromosome.20.fa);
# printFASTA('temp.txt', $reference_href);

print "Reading $genomicMapping\n";
open(MAPPING, '<', $genomicMapping) or die "Cannot open $genomicMapping";
while(<MAPPING>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	
	next unless($line);
	
	my @references = split(/,/, $line);
	@references = grep {$_ =~ /chr6/} @references;
		
	die if($#references > 0);
	next if($#references == -1);
	
	die "Can't parse $references[0]" unless($references[0] =~ /^chr6:(\d+)$/);
	
	my $position = $1;
	$block_positions_chr6{$position} = 1;
}
close(MAPPING);

my $reference_href = readFASTA($inputGenome);

print "Now block ", scalar(keys %block_positions_chr6), " positions.\n";

my $found_6 = 0;
foreach my $chromosomeID (keys %$reference_href)
{
	if($chromosomeID =~ /^6\s/)
	{
		die if($found_6);
		$found_6++;
		for(my $i = 0; $i < length($reference_href->{$chromosomeID}); $i++)
		{
			if($block_positions_chr6{$i})
			{
				substr($reference_href->{$chromosomeID}, $i, 1) = 'N';
			}
		}
	}
}	
unless($found_6)
{
	die "Could not find chromosome 6.\n";
}	

print "Produce output file $outputGenome ...\n";

printFASTA($outputGenome, $reference_href);

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

sub printFASTA
{
	my $outputFile = shift;
	my $sequences_href = shift;
		
	open(F, '>', $outputFile) or die "Cannot open $outputFile";
	foreach my $sequenceID (sort keys %$sequences_href)
	{
		print F '>', $sequenceID, "\n";
		my $sequence = $sequences_href->{$sequenceID};
		while(length($sequence) > 0)
		{
			my $want_length = 100; 
			if(length($sequence) < $want_length)
			{
				$want_length = length($sequence);
			}
			my $line = substr($sequence, 0, $want_length);
			
			die unless(length($line) == $want_length);
			
			print F $line, "\n";
			
			substr($sequence, 0, $want_length) = '';
		}
	}
	close(F);
}