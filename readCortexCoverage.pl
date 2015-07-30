#!/usr/bin/perl

BEGIN {
	push(@INC, 'Toolbox');
} 

# Comments denoted by ALL-HLA refer to changes which would become necessary if we wanted
# to create a graph for all HLA loci at the same time.

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;
use Sys::Hostname;

my $cortex_bin; # = '/Net/x9000/projects/gsk_hla/CortexGraphs/AA02O9Q_A1_55.ctx';
# $cortex_bin = '/Net/x9000/projects/gsk_hla/CortexGraphs/AA02O9Q_A1_55.ctx'; # debug
my $interestingkMers; # = '../tmp2/GS_nextGen/training_A.txt.leucodog_p_probabilistic.finalgraph.augmented.kmers.requiredKMers';
my $kMer_size = 55;

GetOptions ( 
	'cortex_bin:s' => \$cortex_bin,
	'interesting_kMers:s' => \$interestingkMers,
	'kMer_size:s' => \$kMer_size
);

my $output_file = $interestingkMers.'.binaryCount';

my %ACGT2INT = ('A' => 0, 'C' => 1, 'G' => 2, 'T' => 3);
my %INT2ACGT;
foreach my $nucleotide (keys %ACGT2INT) {$INT2ACGT{$ACGT2INT{$nucleotide}} = $nucleotide} 

my $sixtyFourBitUnits;
if($kMer_size > 0 and $kMer_size <= 31)
{
	$sixtyFourBitUnits = 1;
}
elsif($kMer_size > 31 and $kMer_size <= 63)
{
	$sixtyFourBitUnits = 2;
}
else
{
	die "Not yet implemented";
}

my $nullBitString = 0;
for(my $bit = 0; $bit < (64*$sixtyFourBitUnits)/2; $bit++)
{
	vec($nullBitString, $bit, 2) = 0;
}
my $maxK = (64*$sixtyFourBitUnits)/2;
my $firstPositionInBinary = $maxK - $kMer_size;
my $bitPrintFormat = '%#0'.(64*$sixtyFourBitUnits).'b';

die "Cortex binary not existing: $cortex_bin" unless (-e $cortex_bin);
die "Interesting kMers file not existing: $interestingkMers" unless (-e $interestingkMers); # debug

# some testing
foreach my $c (qw/A C G T/)
{
	my $kMer = ($c x $kMer_size);
	my $binary_kMer = binaryRepresentation($kMer);
	my $recovered_kMer = normalRepresentation($binary_kMer);
	die unless($recovered_kMer eq $kMer);
	
	#print "Original kMer: ", $kMer, "\n";
	#print "Binary representation: ", unpack("b*",$binary_kMer), "\n";
	#print "Recovered kMer: ", $recovered_kMer, "\n";
	#print "\t\tATTENTION!\n\n" unless($recovered_kMer eq $kMer);	
}

my @random_kMers;
my @nucleotides = qw/A C G T/;
for(my $i = 0; $i < 100; $i++)
{
	my $kMer = '';
	for(my $j = 0; $j < $kMer_size; $j++)
	{
		my $random_choice = int(rand(4));
		$kMer .= $nucleotides[$random_choice];
	}
	push(@random_kMers, $kMer);
}
foreach my $kMer (@random_kMers)
{
	my $binary_kMer = binaryRepresentation($kMer);
	my $recovered_kMer = normalRepresentation($binary_kMer);
	die unless($recovered_kMer eq $kMer);
}

#die unless(reverseComplement('AAAAA') eq 'TTTTT');
#die unless(reverseComplement('ACGT') eq 'ACGT');
#die unless(reverseComplement('CCCC') eq 'GGGG');
#die unless(reverseComplement('GGGG') eq 'CCCC');
#die unless(reverseComplement('TTTT') eq 'AAAA');

my @interesting_kMers;
my %interesting_kMer_results;
my %interesting_binary_kMers;
#if(1 == 0) debug
#{  debug
my $focal_kMer = '';

print "Reading $interestingkMers\n";
open(KMERS, '<', $interestingkMers) or die "Cannot open $interestingkMers";
while(<KMERS>)
{
	chomp;
	$_ =~ s/\n//g;
	$_ =~ s/\r//g;
	
	my $kMer = $_;
	push(@interesting_kMers, $kMer);
	my $binary_kMer_1 = binaryRepresentation($kMer);
	my $binary_kMer_2 = binaryRepresentation(reverseComplement2($kMer));
	# die if((exists $interesting_binary_kMers{$binary_kMer_1}) || (exists $interesting_binary_kMers{$binary_kMer_2}));
	push(@{$interesting_binary_kMers{$binary_kMer_1}}, $kMer);
	push(@{$interesting_binary_kMers{$binary_kMer_2}}, $kMer);
	if($kMer eq $focal_kMer)
	{
		print "Found kMer of interest in interesting kMers file...\n";
	}
}
close(KMERS);
#} # debug

print "\t...done. Looking for ", scalar(@interesting_kMers), " kMers. \n\n";

# read cortex binary

my $int_length = length(pack( 'I', 0 ));
my $ulongint_length = length(pack( 'L', 0));
my $longlongint_length = length(pack( 'q', 0));

my $buffer;

open(CORTEX, '<', $cortex_bin) or die "Cannot open $cortex_bin";
my $first_buffer_read = read(CORTEX, $buffer, 100);
die "Buffer reading problem: $first_buffer_read" unless($first_buffer_read == 100);
my $cortex_first = index($buffer, 'CORTEX');
my $cortex_second = index($buffer, 'CORTEX', $cortex_first+1);
die unless(($cortex_first != -1) and ($cortex_second != -1));

my $first_header_fields = $cortex_first+length('CORTEX');
my $last_header_fields = $cortex_second-1;
 
# check some header information


my $header_fields = substr($buffer, $first_header_fields , $last_header_fields-$first_header_fields+1);

#for(my $i = 0; $i < length($header_fields)/2; $i++)
#{
#	print "F ", $i, ": ", unpack('I', substr($header_fields, $i*$int_length, $int_length)), "\n";
#}


my $version = unpack('I', substr($header_fields, 0*$int_length, $int_length));
die "Version $version is not supported - I want version >= 4 binaries!\n" unless ($version >= 4);
my $kmer_size_in_file = unpack('I', substr($header_fields, 1*$int_length, $int_length));
die "Wrong k in file: $kmer_size_in_file -- I want $kMer_size!\n" unless ($kmer_size_in_file == $kMer_size);
my $num_of_bitfields_file = unpack('I', substr($header_fields, 2*$int_length, $int_length));
die "Wrong bitfield number $num_of_bitfields_file -- I want $sixtyFourBitUnits!\n" unless ($num_of_bitfields_file == $sixtyFourBitUnits);

my $colors = unpack('I', substr($header_fields, 3*$int_length, $int_length));
die "Wrong number of colors in file: $colors -- I want 1!\n" unless ($colors == 1);

my $mean_read_len = unpack('I', substr($header_fields, 4*$int_length, $int_length));
my $total_seq = unpack('q', substr($header_fields, 5*$int_length, $longlongint_length));

my $beginning_of_cortex_table = $last_header_fields + 1 + length('CORTEX');

my $len_of_kMer_field = (64*$sixtyFourBitUnits)/8;
my $len_of_coverage_field = ($version == 4) ? $ulongint_length : $ulongint_length;
my $len_of_edges_field = 1;
my $total_length_per_kMer = $len_of_kMer_field + $len_of_coverage_field + $len_of_edges_field;

seek CORTEX, $beginning_of_cortex_table, 0;
my $read_at_once = $total_length_per_kMer*10000000; 
my $total_kMer_counter = 0;
$| = 1;
my $total_kmer_coverage = 0;
while(read(CORTEX, $buffer, $read_at_once))
{
	
	if((length($buffer) % $total_length_per_kMer) != 0)
	{
		warn "Weird return value from read: not equal to multiplicity of per-kMer length!\n";
	}
	
	my $received_kMers = int(length($buffer) / $total_length_per_kMer);
	
	for(my $kMerI = 0; $kMerI < $received_kMers; $kMerI++)
	{
		$total_kMer_counter++;
		
		my $kMer_string = substr($buffer, $kMerI*$total_length_per_kMer, $len_of_kMer_field);
		if($sixtyFourBitUnits > 1)
		{
			if ($sixtyFourBitUnits == 2)
			{
				$kMer_string = substr($kMer_string, ($len_of_kMer_field/2)).substr($kMer_string, 0, $len_of_kMer_field/2);
				(length($kMer_string) == $len_of_kMer_field) or die(length($kMer_string));
			}
			else
			{
				die;
			}
		}
		
		my $coverage_string = substr($buffer, $kMerI*$total_length_per_kMer+$len_of_kMer_field, $len_of_coverage_field);
		my $coverage = unpack(($version == 4) ? 'I' : 'L', $coverage_string);
		$total_kmer_coverage += $coverage; 
		
		#my $kMer_decoded = normalRepresentation($kMer_string);
		#print $kMer_decoded, " ", $coverage, "\n";
		#print reverseComplement($kMer_decoded), " ", $coverage, "\n";
		#print reverseComplement2($kMer_decoded), " ", $coverage, "\n";
		#print "\n";
				
		if(exists $interesting_binary_kMers{$kMer_string})
		{
			# my $edges_string = substr($buffer, $kMerI*$total_length_per_kMer+$len_of_kMer_field+$len_of_coverage_field, $len_of_edges_field);
			# my $kMer_decoded = normalRepresentation($kMer_string);
					
			foreach my $unencoded_kMer (@{$interesting_binary_kMers{$kMer_string}})
			{
				$interesting_kMer_results{$unencoded_kMer} = $coverage;
			}
			
			#if($interesting_binary_kMers{$kMer_string} eq $focal_kMer)
			#{
			#	print "Found kMer of interest in the cortex binary with coverage $coverage\n";
			#}
		}
		
		if(($total_kMer_counter % 5000000) == 0)
		{
			print "\r $total_kMer_counter    ";
		}
		 		
	}	
}
close(CORTEX);

open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
foreach my $kMer (@interesting_kMers)
{
	# next unless(exists $interesting_kMer_results{$kMer});
	my $coverage = (exists $interesting_kMer_results{$kMer}) ? $interesting_kMer_results{$kMer} : '00';
	print OUTPUT $kMer, ' ', $coverage, "\n";
}
print OUTPUT 'MeanReadLen', ' ', $mean_read_len, "\n";
print OUTPUT 'TotalKMerCoverage', ' ', $total_kmer_coverage, "\n";
print OUTPUT 'TotalSeq', ' ', $total_seq, "\n";

close(OUTPUT);

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	
	#return reverse($kMer);
	
	return $kMer;
}

sub reverseComplement2
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	
	return scalar(reverse($kMer));
	
	# return $kMer;
}



sub binaryRepresentation
{
	my $kMer = shift;
	die "Received kMer of length ".length($kMer)." while expecting $kMer_size" unless(length($kMer) == $kMer_size);

	my $binary_kMer = $nullBitString;
	
	$kMer = reverse($kMer);
	
	for(my $i = 0; $i < $kMer_size; $i++)
	{
		my $nucleotide = substr($kMer, $i, 1);
		die unless (exists $ACGT2INT{$nucleotide});
		my $binary_value = $ACGT2INT{$nucleotide};
		vec($binary_kMer, $i, 2) = $binary_value;
	}
	
	return $binary_kMer;
}

sub normalRepresentation
{
	my $binary_kMer = shift;
	my $kMer = '';
	
	for(my $i = 0; $i < $kMer_size; $i++)
	{
		my $binary_value = vec($binary_kMer, $i, 2);
		die unless (exists $INT2ACGT{$binary_value});
		my $nucleotide =  $INT2ACGT{$binary_value};
		$kMer .= $nucleotide;
	}
	
	die unless(length($kMer) == $kMer_size);
	
	return reverse($kMer);
}