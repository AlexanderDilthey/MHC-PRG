#!/usr/bin/perl

use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use Data::Dumper;
use Getopt::Long;   
use File::Basename;

my $k = 31;

my $file_target_kMers_details = qq(/gpfs1/well/gsk_hla/debugNA12891/DQB_target_kMers.txt);
my $outputfile = qq(/gpfs1/well/gsk_hla/debugNA12891/extended_kMers.txt);
my @cortex_binaries = (
	qq(/gpfs1/well/gsk_hla/CortexGraphs/referenceGenome_genesBlocked.ctx),
	qq(/gpfs1/well/gsk_hla/CortexGraphs/NA12891.ctx),
	qq(/gpfs1/well/gsk_hla/CortexGraphs/NA12878.ctx)
);

my $cortex_binaries_str;
GetOptions ('file_target_kMers_details:s' => \$file_target_kMers_details,
 'cortex_binaries:s' => \$cortex_binaries_str, 
 'outputfile:s' => \$outputfile, 
);

if($cortex_binaries_str)
{
	@cortex_binaries = split(/,/, $cortex_binaries_str);
}

die "Please specify --file_target_kMers_details" unless($file_target_kMers_details);
die "Please specify --outputfile" unless($outputfile);
die "Please specify --cortex_binaries bin1[,bin2] .." unless(scalar(@cortex_binaries));

open(INKMERS, '<', $file_target_kMers_details) or die "Cannot open $file_target_kMers_details";
my $headerLine = <INKMERS>;
chomp($headerLine);
$headerLine =~ s/\n//g;
$headerLine =~ s/\r//g;
my @header_fields = split(/\t/, $headerLine);
my @inkMers_contentLines;
my %kMers;
while(<INKMERS>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	next unless($line);
	
	my @fields = split(/\t/, $line, -1);
	die "Wrong field length $#fields $#header_fields" unless($#fields == $#header_fields);
	my %line = (mesh @header_fields, @fields);
	
	die "No kMer?" unless($line{'kMer'});
	die "Wrong kMer length ".length(length($line{'kMer'}))." ($line{'kMer'})in $file_target_kMers_details"  unless(length($line{'kMer'}) == $k);
	$kMers{$line{'kMer'}}++;
	
	push(@inkMers_contentLines, \%line);
}

my %kMers_by_ctxFile;
my @kMers = keys %kMers;

foreach my $ctxFile (@cortex_binaries)
{
	my $kMers_file = $ctxFile . '._kmers';
	open(KMERS, '>', $kMers_file) or die "Cannot open $kMers_file";
	print KMERS join("\n", @kMers);
	close(KMERS);
	
	my $binaryCount = $kMers_file . '.binaryCount';
	if(-e $binaryCount)
	{
		unlink($binaryCount) or die "Cannot unlink $binaryCount";
	}
	
	unless(-e '../readCortexCoverage.pl')
	{
		die '../readCortexCoverage.pl not existing - call me from right directory!';
	}
	
	my $count_command = qq(perl ../readCortexCoverage.pl --kMer_size $k --cortex_bin $ctxFile --interesting_kMers $kMers_file);
	print "Now executing:\n\t$count_command\n";
	my $ret = system($count_command);
	die "Command $count_command failed" unless($ret == 0);
	
	unless(-e $binaryCount)
	{
		die "Can't find $binaryCount";
	}
	
	unlink($kMers_file) or die "Cannot unlink $kMers_file";
	
	open(KMERCOUNTS, '<', $binaryCount) or die "Cannot open $binaryCount";
	while(<KMERCOUNTS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		
		next if($line =~ /Mean/);
		next if($line =~ /Total/);
		
		die "Weird format -${line}- in $binaryCount" unless($line =~ /^(\w+)\s+(\d+)$/);
		my $kMer = $1;
		my $count = $2;
		
		unless($kMers{$kMer})
		{
			die "Weird kMer -${kMer}- in $binaryCount";
		}
		
		$kMers_by_ctxFile{$ctxFile}{$kMer} = $count;
	}
	close(KMERCOUNTS);
}

my %__existing_header_fields = map {$_ => 1} @header_fields;
open(OUTKMERS, '>', $outputfile) or die "Cannot open $outputfile";

my @output_header_fields = @header_fields;
foreach my $ctxFile (@cortex_binaries)
{
	my $basename = (fileparse($ctxFile))[0];
	$basename =~ s/\..+//;
	push(@output_header_fields, $basename);
}
print OUTKMERS join("\t", @output_header_fields), "\n";

foreach my $line_href (@inkMers_contentLines)
{
	my @line_outputfields = map {$line_href->{$_}} @header_fields;
	
	foreach my $ctxFile (@cortex_binaries)
	{	
		die "No kMer count for $$line_href->{'kMer'} and $ctxFile" unless(defined $kMers_by_ctxFile{$ctxFile}{$line_href->{'kMer'}});
		push(@line_outputfields, $kMers_by_ctxFile{$ctxFile}{$line_href->{'kMer'}});
	}
	
	print OUTKMERS join("\t", @line_outputfields),  "\n";
}
close(OUTKMERS);

print "\nDone. Generated $outputfile \n";