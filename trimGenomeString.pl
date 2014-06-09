#!/usr/bin/perl
use strict;
use Getopt::Long;   

my $file;
my $start = 2000000;
my $stop = 2500000;
GetOptions (
	'file:s' => \$file,
	'start:s' => \$start,
	'stop:s' => \$stop,
);

my $accumulated_length = 0;
my $printed_length = 0;
open(F, '<', $file) or die "Cannot open $file";
while(<F>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @fields = split(/,/, $line, -1);
	die unless($#fields == 1);
	if(length($fields[1]) > 0)
	{
		die unless(length($fields[0]) == length($fields[1]));
	}
	my $L = length($fields[0]);
	my $start_this_segment = $accumulated_length;
	my $stop_this_segment = $accumulated_length + $L - 1;
	
	if((($start_this_segment >= $start) and ($start_this_segment <= $stop)) or
	   (($stop_this_segment >= $start) and ($stop_this_segment <= $stop)) or
	   (($start_this_segment <= $start) and ($stop_this_segment >= $stop)))
	{
		my $remove_front_characters = ($start_this_segment < $start) ? ($start - $start_this_segment) : 0;
		my $remove_back_characters = ($stop_this_segment > $stop) ? ($stop_this_segment - $stop) : 0;
		
		die unless($remove_front_characters >= 0);
		die unless($remove_back_characters >= 0);
		
		my $f_1 = $fields[0];
		my $f_2 = $fields[1];
		
		$f_1 = substr($f_1, $remove_front_characters, length($f_1) - $remove_front_characters - $remove_back_characters);
		if($f_2)
		{
			$f_2 = substr($f_2, $remove_front_characters, length($f_2) - $remove_front_characters - $remove_back_characters);		
		}
		
		print $f_1, ',', $f_2, "\n";	
		
		$printed_length += length($f_1);
	}
	
	$accumulated_length += $L;
}
close(F);

print STDERR "Printed $printed_length characters\n";



