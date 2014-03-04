#!/usr/bin/perl

use strict;
use Data::Dumper;

# params

my $VCF_file = '../data/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.xMHC';
my $path_to_PGF_haplotype = qq(/Net/x9000/projects/gsk_hla/shared/mhc_ref_8_haplotypes/pgf_ref.fasta);
my $graph_file = '../tmp2/GS_nextGen/varigraph_new/graph.txt';

# pgf

my $pgf_start = 28702185;
my $pgf_length = length(read_PGF());
die unless($pgf_length);
my $pgf_end = $pgf_start + $pgf_length - 1;

# variables

my %known_variant_alleles;
my %have_variant_alleles;
my %snp_pos;

# get SNPs

read_SNP_VCF();

# check graph

my $graph_string;
open(GRAPH, '<', $graph_file) or die "Cannot open $graph_file";
my $mode = -1;
while(<GRAPH>)
{
	my $line = $_;
	if($line =~ /CODE\:/)
	{
		$mode = 1;
	}
	elsif($line =~ /NODES\:/)
	{
		$mode = 2;
	}
	elsif($line =~ /EDGES\:/)
	{
		$mode = 3;
	}
	else
	{
		if($mode == 1)
		{
			chomp($line);
			next unless($line =~ /^rs/);			
			my @f = split(/\|\|\|/, $line);
			die unless($#f == 2);
			$f[0] =~ s/\_.+//g;
			$have_variant_alleles{$f[0]}{$f[1]}++;
		}
	}
}
close(GRAPH);

my $snps_found;
my $snps_searched;   
my $alleles_found;
my $alleles_searched; 

foreach my $rsID (sort {$snp_pos{$a} <=> $snp_pos{$b}} grep {($snp_pos{$_} > 2537944) and ($snp_pos{$_} < 2619711)} keys %known_variant_alleles)
{
	$snps_searched++;

	if(exists $have_variant_alleles{$rsID})
	{
		$snps_found++;
	}
	else
	{
		next;
	}
	
	foreach my $allele (keys %{$known_variant_alleles{$rsID}})
	{
		if(exists $have_variant_alleles{$rsID}{$allele})
		{
			$alleles_found++;
		}		
		else
		{
			print "Missing allele $allele for SNP $rsID [$snp_pos{$rsID}]\n";
		}
		$alleles_searched++;
	}
}

print "\n\nFound $snps_found / $snps_searched SNPs.\n";
print "Found $alleles_found / $alleles_searched alleles for these SNPs.\n";

sub read_SNP_VCF
{
	open(VCF, '<', $VCF_file) or die "Cannot open $VCF_file";
	my $id_counter = 0;
	
	while(<VCF>)
	{
		my $line = $_;
		chomp($line);
		next if(substr($line, 0, 2 ) eq '##');
		if(substr($line, 0, 1) eq '#')
		{
			substr($line, 0, 1) = '';
			my @header_fields = split(/\t/, $line);
			die unless(($header_fields[0] eq 'CHROM') and ($header_fields[1] eq 'POS') and ($header_fields[2] eq 'ID') and ($header_fields[3] eq 'REF') and ($header_fields[4] eq 'ALT'));
		}
		else
		{
			my @fields = split(/\t/, $line);
			my $ref_allele = $fields[3];
			my $alt_alleles = $fields[4];
			my @alt_alleles = split(/,/, $alt_alleles);
			die unless($fields[0] eq '6');
			
			next unless(($fields[1] >= $pgf_start) and ($fields[1] <= $pgf_end));

			next unless(length($ref_allele) == 1);
			my @alt_alleles_SNP = grep {length($_) == 1} @alt_alleles;
			
			if($#alt_alleles_SNP > -1)
			{
				my $id = $fields[2];
				if($id eq '.')
				{
					$id_counter++;
					$id = 'G1000_'.$id_counter;
					next;
				}
				
				$snp_pos{$id} = $fields[1] - $pgf_start;
		
				if($known_variant_alleles{$id})
				{
					warn "$id appears more than once in VCF!";
				}
				
				
				$known_variant_alleles{$id}{$ref_allele}++;
				foreach my $allele (@alt_alleles)
				{
					$known_variant_alleles{$id}{$allele}++;
				}
			}
				
		}
	}
	close(VCF);
}



sub read_PGF
{
	my %alignment;
	my $file = $path_to_PGF_haplotype;
	open(ALIGNMENT, "<", $file) or die "Cannot open $file";
	my $current_name;
	while(<ALIGNMENT>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		if($line =~ /^>/)
		{
			die if ($current_name);
			$current_name = 'pgf';
		}
		else
		{
			die unless($current_name);
			$alignment{$current_name} .= $line;
		}
	}
	close(ALIGNMENT);
	$alignment{'pgf'} or die;
	return $alignment{'pgf'};
}
