#!/usr/bin/perl
use Modern::Perl;
use Getopt::Long;   
use List::MoreUtils qw/mesh/;
use Data::Dumper;

$| = 1;

my $inputGenome = qq(/gpfs1/well/gsk_hla/GRCh37.60/fasta/combined/Homo_sapiens.GRCh37.60.dna.chromosome.ALL.fa);
my $refSeq_UCSC = qq(/gpfs1/well/gsk_hla/debugNA12891/RefSeq_B37_July2015_fromUCSC.txt);
my $addGenePositions = qq(/gpfs1/well/gsk_hla/debugNA12891/manualGenes_UCSCFormat.txt);
my $padding = 30;

my $outputGenomeFile;
my $geneList;

GetOptions ('outputGenomeFile:s' => \$outputGenomeFile,
 'geneList:s' => \$geneList, 
);

unless($outputGenomeFile)
{
	die "Please specify --outputGenomeFile";
}

unless($geneList)
{
	die "Please specify --geneList";
}

# read positions

my %gene_positions;
my %double_genes;
foreach my $file ($refSeq_UCSC, $addGenePositions)
{
	open(GENEPOS, '<', $file) or die "Cannot open $file";
	my $headerLine = <GENEPOS>;
	chomp($headerLine);
	$headerLine =~ s/[\n\r]//g;
	my @header_fields = split(/\t/, $headerLine);
	while(<GENEPOS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless($#line_fields == $#header_fields);
		
		my %line = (mesh @header_fields, @line_fields);
		
		my $geneID = $line{'name2'};
		my $chr = $line{'chrom'};
		
		next if($chr =~ /hap/);
		$chr =~ s/chr//;
		
		my $start = $line{'cdsStart'};
		my $stop = $line{'cdsEnd'};
		
		my $NM = $line{'name'};
		die "Weird name: $NM" unless($NM =~ /(NM_)|(NR_)/);
		$NM =~ s/NM_//;
		$NM =~ s/NR_//;
		
		die "Weird coordinates for $geneID : $start $stop line $." unless($start <= $stop);
		
		if (exists $gene_positions{$geneID})
		{
			unless(	($gene_positions{$geneID}[0] eq $chr) and
					($gene_positions{$geneID}[1] eq $start) and
					($gene_positions{$geneID}[2] eq $stop))
			{
				$double_genes{$geneID}++;
			}
			# warn "Double gene ID $geneID";
			next if($gene_positions{$geneID}[2] < $NM);
		}
		$gene_positions{$geneID} = [$chr, $start, $stop, $NM];
	}
	close(GENEPOS);
}

# read blocked genes

my @blockGenes; 
print "Reading $geneList\n";
open(GENELIST, '<', $geneList) or die "Cannot open $geneList";
my $headerLine = <GENELIST>;
chomp($headerLine);
$headerLine =~ s/\r//g;
$headerLine =~ s/\n//g;
my @header_fields = split(/\t/, $headerLine);
while(<GENELIST>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	
	next unless($line);
	
	my @fields = split(/\t/, $line);
	my %line = (mesh @header_fields, @fields);
	
	my $geneName = $line{'Locus'};
	die unless($geneName);
	
	if((defined $line{'B37_Start_1based'}) and ($line{'B37_Start_1based'}  ne ''))
	{
		die unless($line{'Chromosome'});
		die unless($line{'B37_Stop_1based'});
		
		die unless($line{'B37_Stop_1based'} >= $line{'B37_Start_1based'});
		push(@blockGenes, [$geneName, $line{'Chromosome'}, $line{'B37_Start_1based'} - 1, $line{'B37_Stop_1based'} - 1]);
		
	}
	else
	{
		unless(exists $gene_positions{$geneName})
		{
			warn "No position for $geneName - ignore.";
			next;
		}
		
		if(exists $double_genes{$geneName})
		{
			warn "Ambiguous position for $geneName";
		}
		
		my $genePos = $gene_positions{$geneName};
		die "Weird positions for $geneName $genePos->[2] $genePos->[1]" unless($genePos->[2] >= $genePos->[1]);
		
		push(@blockGenes, [$geneName, @$genePos]);
	}
}
close(GENELIST);

print "Now block ", scalar(@blockGenes), " genes.\n";

print "Reading $inputGenome\n";
my $reference_href = readFASTA($inputGenome);
print "\t... done!\n";

my $blocked = 0;
foreach my $genePos (@blockGenes)
{
	my $chr = $genePos->[1];
	die "No chromosomal sequence for $chr" unless(exists $reference_href->{$chr});
	
	my $startPos_N = $genePos->[2];
	my $stopPos_N = $genePos->[3];
	die "Weird positions $startPos_N $stopPos_N" unless($stopPos_N >= $startPos_N);
	
	$startPos_N  += $padding;
	$stopPos_N  -= $padding;
	next unless($stopPos_N >= $startPos_N);
	my $substitute_L = $stopPos_N - $startPos_N + 1;
	
	print "Gene $genePos->[0]: Substitute $substitute_L characters with Ns\n";
	
	my $Nstring = ('N' x $substitute_L);
	die unless(length($Nstring) == $substitute_L);
	
	substr($reference_href->{$chr}, $startPos_N, $substitute_L) = $Nstring;
	
	$blocked += $substitute_L;
}

print "\nTotal blocked characters: $blocked\n";

print "Produce output file $outputGenomeFile ...\n";

printFASTA($outputGenomeFile, $reference_href);

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
			$currentSequence =~ s/\s+.+//;
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