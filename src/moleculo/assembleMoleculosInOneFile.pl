#!/usr/bin/perl
use strict;
use Data::Dumper;

my @input_directories = map {'/Net/banyan/data2/projects/1000genomes/NA12878/moleculo/'.$_} qw/good bad ok/;
my $original_alignment = '/gpfs1/well/gsk_hla/shared/mhc_ref_8_haplotypes/alignment/all_aligned.fasta';
my $alignment = read_alignment($original_alignment);

my $output_file = '/Net/banyan/data1/projects/gsk/moleculo/contigs.fasta';
open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";

foreach my $haploID (keys %$alignment)
{	
	print OUTPUT '>TEST_', $haploID, "\n";
	print OUTPUT $alignment->{$haploID}, "\n";		
}

my %haveIDs;
foreach my $dir (@input_directories)
{
	my @files = glob($dir.'/*.fq');
	
	foreach my $file (@files)
	{
		print "Reading $file ... \n";
		
		open(INPUT, '<', $file) or die "Cannot open $file";
		while(<INPUT>)
		{
			my $firstLine = $_;   
			my $secondLine = <INPUT>;
			my $thirdLine = <INPUT>;
			my $fourthLine = <INPUT>;
			
			$secondLine =~ s/X/N/g;
			
			die Dumper('Length', $file, $., length($secondLine), length($fourthLine)) unless(length($secondLine) == length($fourthLine));
			die Dumper('Plus', $file, $.)  unless(substr($thirdLine, 0, 1) eq '+');
			
			die Dumper('At', $file, $.)  unless(substr($firstLine, 0, 1) eq '@');
			die Dumper('RegExp', $file, $.)  unless($firstLine =~ /^\@(.+?)\s+length=/);
			
			my $ID = $1;
			
			$haveIDs{$ID}++;
			
			if($haveIDs{$ID} > 1)
			{
				warn "Double ID $ID in file $file";
				$ID .= '_'.$haveIDs{$ID};
			}
			chomp($secondLine);
			unless($secondLine =~ /^[acgtnACGTN]+$/)
			{
				for(my $i = 0; $i < length($secondLine); $i++)
				{
					my $C = substr($secondLine, $i, 1);
					my $OK = ($C =~ /^[acgtnACGTN]$/);
					print join("\t", $i, $C, $OK), "\n";
					if(not $OK)
					{
						last;
					}
				}	
				die "Can't match offending!" unless($secondLine =~ /([^acgtnACGTN])/);
				my $offending = $1;
				die Dumper('Characters!', $file, $., $offending)  			
			}
			
			print OUTPUT '>', $ID, "\n";
			print OUTPUT $secondLine, "\n";			
		}
		
		close(INPUT);
	}
}

close(OUTPUT);


print "\n\nWritten file $output_file \n\n";

sub read_alignment
{
	# this is code duplicated from documents\analysis\04 Januar 2012\MHC alignment prepare for variation graph 2.pl , and should stay
	# consistent with that! (apart from the bit which determines alignment positions)
	
	my $file = shift;

	my %alignment;
	my %alignment_positions;

	open(ALIGNMENT, "<", $file) or die "Cannot open $file";
	my $current_name;
	while(<ALIGNMENT>)
	{
		my $line = $_;
		chomp($line);

		if($line =~ /^>/)
		{
			if($line =~ /chr6_(\w+?)_/)
			{
				$current_name = $1;
			}
			else
			{
				$current_name = 'pgf';
			}
			die if(exists $alignment{$current_name});
		}
		else
		{
			$line =~ tr/acgtn/ACGTN/;
			die unless($current_name);
			$line =~ s/\./_/g;
			$line =~ s/\-/_/g;
			if($current_name eq 'pgf')
			{
				if($line =~ /n|N/)
				{
					die $line;
				}
			}
			
			$line =~ s/N/\*/g;			
			
			die unless ($line =~ /^[ACGT\*\_]+$/);			
			$line =~ s/_//g;
			$line =~ s/\*//g;
			$alignment{$current_name} .= $line;
		}
	}
	close(ALIGNMENT);
	
	return \%alignment;
}
