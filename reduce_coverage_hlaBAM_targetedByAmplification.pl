#!/usr/bin/perl

# Always assume paired-end reads!

use strict;
use Data::Dumper;
use Getopt::Long;   

my $samtools_bin = find_present_alternative('/opt/sw/software/Bio/samtools/0.1.19/bin/samtools', '/home/dilthey/samtools-0.1.18/samtools', '/apps/well/samtools/0.1.19/bin/samtools');


my $maximum_number_starting_reads = 250;
#my $input_BAM = '/gpfs1/well/gsk_hla/temp_mapping_2/AM_WTSI19556/merged.bam';
#my $output_BAM = '/gpfs1/well/gsk_hla/Alexander_Mentzer_BAMs_downsampled/AM_WTSI19556.bam';

my $input_BAM;
my $output_BAM;


GetOptions (
 'input_BAM:s' => \$input_BAM,  
 'output_BAM:s' => \$output_BAM,  
);

die "Please specify --input_BAM" unless($input_BAM);
die "Please specify --output_BAM" unless($output_BAM);
die "Please specify an existing --input_BAM" unless(-e $input_BAM);
die unless($input_BAM ne $output_BAM);

print "\nTarget (maximum) rate of reads starting per position: ", $maximum_number_starting_reads, "\n\n";

my $tF = getTempFile();

my $cmd_extract = qq($samtools_bin view -h $input_BAM > $tF);

system($cmd_extract) and die "Samtools extraction command failed: $cmd_extract";

print "Extracted SAM into $tF\n";

my $fn_samout = $tF . '.filtered';

my $examined_reads = 0;
my $printed_reads = 0;

open(SAMOUT, '>', $fn_samout) or die "Cannot open $fn_samout";

my %seen_read_IDs;
my %read_IDs_decisions;

my $current_start_position;
my @current_start_position_reads;


my $process_open_reads = sub {

	if(not defined $current_start_position)
	{
		return;
	}

	$current_start_position = undef;	
	return unless(scalar(@current_start_position_reads));
	
	my %read_IDs;
	my $positionID;
	foreach my $read (@current_start_position_reads)
	{
		my $readID = $read->[0];
		$read_IDs{$readID}++;
		
		my $thisRead_posID = $read->[2] . ':' . $read->[3];
		unless(defined $positionID)
		{
			$positionID = $thisRead_posID;
		}
		
		die unless($thisRead_posID eq $positionID);
	}
	
	my %take_read_IDs;
	my $take_reads_rate = $maximum_number_starting_reads / scalar(keys %read_IDs);
	
	if($positionID eq '*:0')
	{
		$take_reads_rate = 1;
	}
	
	if($take_reads_rate >= 1)
	{
		%take_read_IDs = map {$_ => 1} (keys %read_IDs);
	}
	else
	{
		foreach my $readID (keys %read_IDs)
		{
			my $r = rand(1);
			die unless(($r >= 0) and ($r <= 1));
			$take_read_IDs{$readID} = ($r <= $take_reads_rate);
		}
	}
	
	# print "Position ${positionID}, original reads ", scalar(keys %read_IDs), ", take rate ${take_reads_rate}, taken reads ", scalar(grep {$take_read_IDs{$_}} keys %take_read_IDs), "\n";
	
	foreach my $read (@current_start_position_reads)
	{
		my $readID = $read->[0];
		die unless(defined $take_read_IDs{$readID});
		if($take_read_IDs{$readID})
		{
			print SAMOUT join("\t", @$read), "\n";
			$printed_reads++;
		}
	}
	
	foreach my $readID (keys %take_read_IDs)
	{
		die if(defined $read_IDs_decisions{$readID});
		$read_IDs_decisions{$readID} = $take_read_IDs{$readID};
	}
	
	@current_start_position_reads = ();
};	

open(SAMIN, '<', $tF) or die "Cannot open $tF";
while(<SAMIN>)
{	
	my $line = $_;
	chomp($line);
	next unless($line);
	
	if(substr($line, 0, 1) eq '@')
	{
		print SAMOUT $line, "\n";
	}
	else
	{
		$examined_reads++;
		
		my @line_fields = split(/\t/, $line);
		my $readName = $line_fields[0];
		my $reference = $line_fields[2];
		my $position = $line_fields[3];
		my $posID = $reference . ':' . $position;
		
		if($posID ne $current_start_position)
		{
			$process_open_reads->();
		}
		
		$current_start_position = $posID;
		
		$seen_read_IDs{$readName}++;
		if(exists $read_IDs_decisions{$readName})
		{	
			if($read_IDs_decisions{$readName})
			{
				print SAMOUT $line, "\n";
				$printed_reads++;
			}
		}
		else
		{
			push(@current_start_position_reads, \@line_fields);
		}
	}
}

$process_open_reads->();

close(SAMOUT);
close(SAMIN);

print "\nProcessed reads: ${examined_reads}; taken: $printed_reads - rate ", sprintf("%.2f", $printed_reads / $examined_reads ), "\n";

my $lost_read_IDs = 0;
my $split_read_IDs = 0;
foreach my $readID (keys %seen_read_IDs)
{
	my $c = $seen_read_IDs{$readID};
	if($c <= 1)
	{
		$lost_read_IDs++;
	}
	elsif($c >= 3)
	{
		$split_read_IDs++;
	}
}

print "\nRead ID statistics:\n";
print "\tTotal: ", scalar(keys %seen_read_IDs), "\n";
print "\tApparently non-paired: ", $lost_read_IDs, "\n";
print "\tApparently split: ", $split_read_IDs, "\n";

my $cmd_transform_BAM = qq(cat $fn_samout | ${samtools_bin} view -Sb - > ${output_BAM});
system($cmd_transform_BAM) and "BAM conversion command $cmd_transform_BAM failed";

my $cmd_index_BAM = qq(${samtools_bin} index ${output_BAM});
system($cmd_index_BAM) and "BAM conversion command $cmd_index_BAM failed";

print "\nProduced BAM: $output_BAM \n\n";

unlink($tF) or warn "Could not delete $tF";
unlink($fn_samout) or warn "Could not delete $fn_samout";

sub getTempFile
{
	my $tmp_dir = "../tmp";
	die unless(-e $tmp_dir);
	die unless(-d $tmp_dir);
	
	my $fnTmp;
	do {
		my $k = genkey();
		$fnTmp = $tmp_dir . '/' . $k;		
	} while (-e $fnTmp);
	
	return $fnTmp;
}
sub find_present_alternative
{
        foreach my $a (@_)
        {
                if(-e $a)
                {
                        return $a;
                }
        }
        die "Could not find a present alternative from list:\n".join("\n", @_);
}

sub genkey
{
        my $key;
        my $num = $_[0] ? $_[0] : 25;
        for (my $i = 0; $i <= $num; $i++)
        {
                my $val =  int(rand(9)+1);
                $val = ($val > 9) ? 9 : $val;
                $key .= $val;
        }
        return $key;
}


