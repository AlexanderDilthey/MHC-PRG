#!/usr/bin/perl
use Modern::Perl;
use Getopt::Long;   
use List::MoreUtils qw/mesh/;
use Data::Dumper;

my $graph = qq(hla);

my $graph_dir = qq(../../tmp2/GS_nextGen/) . $graph . '/';
my $positions_file = $graph_dir . 'positions.txt';
unless(-e $graph_dir)
{
	die "Inferred graph directory $graph_dir not existing"
}
unless(-e $positions_file)
{
	die "Inferred graph directory file $positions_file not existing"
}

my %genes;
open(POSITIONS, '<', $positions_file) or die "Cannot open $positions_file";
while(<POSITIONS>)
{
	my $line = $_;
	die unless($line =~ /^(.+?)\/(.+)/);
	my $gene = $1;
	my $components = $2;
	next if ($gene =~ /before/);
	my @parts = split(/_/, $components);
	my $name_component = $parts[0];
	if($name_component =~ /((exon)|(intron))\/(\d)/)
	{
		my $intronExon = $1;
		my $number = $4;
		$genes{$gene}{$intronExon}{$number}++;
		# die Dumper(scalar(keys %{$genes{$gene}{$intronExon}}), $number, $line, $genes{$gene}) unless(scalar(keys %{$genes{$gene}{$intronExon}}) == $number);
	}
	

}
close(POSITIONS);

# print join("\t", "Gene", "Introns", "Exons"), "\n";
foreach my $gene (sort keys %genes)
{
	my @components = sort keys %{$genes{$gene}};
	print $gene, "\n";
	print join("\n", map {' - ' . $_ . ' ' . join(', ', sort keys %{$genes{$gene}{$_}})} @components), "\n";
	print "\n";
	
	# print join("\t", $gene, scalar(keys %{$genes{$gene}{'intron'}}), scalar(keys %{$genes{$gene}{'exon'}})), "\n";
}
