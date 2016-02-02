#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

my $MFA;
my $graphName;

GetOptions ('MFA:s' => \$MFA,
'graphName:s' => \$graphName,
);

unless($MFA)
{
	die "Please specify input MFA in FASTA format with --MFA";
}

unless($graphName)
{
	die "Please specify output graph name --graphName";
}

unless(-e $MFA)
{
	die "Specified MFA $MFA not existing!";
}

my $baseDir_graphs = '../tmp2/GS_nextGen';
unless(-e $baseDir_graphs)
{
	die "Please call me from MHC-PRG src directory - $baseDir_graphs not existing.";
}

my $seq_href = readFASTA($MFA);

my $MFA_basename = basename($MFA);

unless(scalar(keys %$seq_href))
{
	die "No sequences in $MFA?";
}
my $L;
my %S_for_output;
foreach my $seqID (keys %$seq_href)
{
	my $S = $seq_href->{$seqID};
	unless(defined $L)
	{
		$L = length($S);
	}
	die "Length problem with $seqID - not equal to $L, which is assumed consensus length" unless($L == length($S));
	
	die "Weird characteters in $seqID" unless($S =~ /^[ACGTN\-]+$/i);
	$S = uc($S);
	$S =~ s/-/_/g;
	my $id_for_output = $seqID;
	$id_for_output =~ s/\s/_/g;	
	$S_for_output{$id_for_output} = $S;
}

my $outputDir = $baseDir_graphs . '/' . $graphName;
mkdir($outputDir);
die unless(-e $outputDir);

my @positionIDs;
my $f_positions = $outputDir . '/positions.txt';
open(POSITIONS, '>', $f_positions) or die "Cannot open $f_positions";
for(my $pI = 0; $pI < $L; $pI++)
{
	my $pName = 'P_' . $pI . '_' . ($pI+1);
	print POSITIONS join(' ', $pName, ($pI+1)), "\n";
	push(@positionIDs, $pName);
}
close(POSITIONS);

my $f_segments = $outputDir . '/segments.txt';
my $fn_alignment = 'PRG_' . $MFA_basename;

open(SEGMENTS, '>', $f_segments) or die "Can't open $f_segments";
print SEGMENTS $fn_alignment, "\n";
close(SEGMENTS);

my $f_alignment = $outputDir . '/' . $fn_alignment;
open(ALIGNMENTOUT, '>', $f_alignment) or die "Cannot open $f_alignment";
print ALIGNMENTOUT join(' ', 'IndividualID', @positionIDs), "\n";
foreach my $key (keys %S_for_output)
{
	my @f = split(//, $S_for_output{$key});
	die unless($#f == $#positionIDs);
	print ALIGNMENTOUT join(' ', $key, @f), "\n";
}
close(ALIGNMENTOUT);

print "Generated input for graph construction in $outputDir\n";

my $cmd_suggestion = qq(../bin/MHC-PRG domode createConcatenatedVariationGraphs ${baseDir_graphs}/${graphName} --noPGFprotection);

print "\nSuggested command\n";
print $cmd_suggestion, "\n\n";

sub readFASTA
{
	my $file = shift;	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
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
		
	return \%R;
}
