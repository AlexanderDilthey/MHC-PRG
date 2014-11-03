#!/usr/bin/perl

use strict;
use Data::Dumper;

my $output_directory = '_displayReadsAlignmentOutput';
unless(-e $output_directory)
{
	mkdir($output_directory) or die;
}	

my %interestingReadsBySample = (	
	# 'I1_AA02O9Q_B6' => {
		# 'Unknown' => [
			# '@@A8197TABXX:4:2208:4575:106554#GCCAATAT/1:FROM:6:31323963:FROM:'
		# ],
	# },
	# 'I1_AA02O9Q_B1' => {
		# 'mapperSaysThere' => [
			# '@@B819B8ABXX:2:2203:4373:99705#CGATGTAT/2:FROM:6:32489862:FROM:',
			# '@@B819B8ABXX:2:1203:10388:38266#CGATGTAT/2:FROM:6:32552015:FROM:',
			# '@@B819B8ABXX:2:1205:11083:89762#CATGTATC/2:FROM:6:32552024:FROM:',
			# '@@A8176BABXX:1:1101:11794:63897#CGATGTAT/2:FROM:6:32489763:FROM:',
			# '@@B819B8ABXX:2:1205:11083:89762#CATGTATC/1:FROM:6:32551957:FROM:',
		# ],
		# 'trueHLAsaysThere' => [
			# '@@B819B8ABXX:2:1204:2500:66314#CGATGTAT/2:FROM:6:32552025:FROM:',
			# '@@A8176BABXX:1:1206:2915:25661#CGATGTAT/1:FROM:6:32552024:FROM:',
			# '@@A8176BABXX:1:1104:18107:177078#CGATGTAT/2:FROM:6:32525833:FROM:',
			# '@@FCC00N2ABXX:8:1203:18982:24309#CGATGTAT/2:FROM:6:32489672:FROM:',
		# ],
	# },
	
	# 'I1_AA02O9Q_B3' => {
		# 'mapperSaysThere' => [
			# '@@B81998ABXX:6:1108:16939:81798#GCCAATAT/1:FROM:6:32489762:FROM:',
			# '@@B81998ABXX:6:1108:16939:81798#GCCAATAT/2:FROM:6:32489698:FROM:',
			# '@@B819B9ABXX:8:1105:13324:120581#CGCCNNNN/1:FROM:6:32551913:FROM:',
			# '@@B819B9ABXX:8:1105:13324:120581#CGCCNNNN/1:FROM:6:32551913:FROM:',
			# '@@B81998ABXX:6:1207:2569:100263#GCCNATAT/2:FROM:6:32551912:FROM:',
			# '@@B81998ABXX:6:1108:16939:81798#GCCAATAT/2:FROM:6:32489698:FROM:',
			# '@@B819B9ABXX:8:1105:13324:120581#CGCCNNNN/1:FROM:6:32551913:FROM:',
			# '@@B819B9ABXX:8:2208:3453:9743#GCCAATAT/1:FROM:6:32551903:FROM:',
			# '@@B81998ABXX:6:1108:16939:81798#GCCAATAT/2:FROM:6:32489698:FROM:',
			
		# ],
		# 'trueHLAsaysThere' => [	
			# '@@B81EP5ABXX:7:1107:20989:92286#GCCAATAT/1:FROM:6:32489739:FROM:',
			# '@@A819GPABXX:5:2107:18518:148281#GCCAATAT/1:FROM:6:32551968:FROM:',
			# '@@B819B9ABXX:8:1108:4618:200340#GCCANNNN/1:FROM:6:32489707:FROM:',
		# ],
	# },
	'I1_AA02O9Q_A4' => {
		'mapperSaysThere' => [
			'@@A819BJABXX:3:2206:5747:176757#ATCACGAT:normalAlignment=8.70858e-05[6:29856481-29856570];1.13446e-18[6:29856377-29856466];1;15/2:FROM:6:29856377:FROM:',
			'@@A819BJABXX:2:1205:8561:17512#ATCACGAT:normalAlignment=2.43113e-07[6:29856511-29856600];4.29519e-20[6:29856423-29856512];1;-1/2:FROM:6:29856423:FROM:',
			'@@A819GPABXX:4:2108:3436:12394#ATCACGAT:normalAlignment=1.17013e-29[6:29856451-29856540];2.70154e-08[6:29856537-29856626];1;-3/1:FROM:6:29856451:FROM:',
			
		],
		'trueHLAsaysThere' => [	
		],
	},	
);

# pre-processing

my %reads_2_group;
foreach my $sampleID (keys %interestingReadsBySample)
{
	foreach my $groupName (keys %{$interestingReadsBySample{$sampleID}})
	{
		my @readIDs = @{$interestingReadsBySample{$sampleID}{$groupName}};
		foreach my $readID (@readIDs)
		{
			$reads_2_group{$readID} = $groupName;
		}
	}
}

# C++ part
{
	my $cpp_output_path = $output_directory . '/' . 'cppmap.txt';
	open(CPP, '>', $cpp_output_path) or die "Cannot open $cpp_output_path";
	print CPP qq(std::map<std::string, std::string> readID_2_group;), "\n";			
	foreach my $readID (keys %reads_2_group)
	{
		print CPP qq(readID_2_group["$readID"] = "$reads_2_group{$readID}";), "\n";
	}
	close(CPP);
}  

# main part
my $output_file_sequenceOnly = $output_directory . '/' . 'sequenceOnly.txt';
open(SEQUENCEONLY, '>', $output_file_sequenceOnly) or die "Cannot open $output_file_sequenceOnly";

print "Printing to $output_file_sequenceOnly\n";

foreach my $sampleID (keys %interestingReadsBySample)
{
	my $input_filename = qq(../tmp/hla/${sampleID}/reads.p.n.aligned);
	
	my %output_lines_perGroup;
	
	my %specifiedRead_2_firstReadID;
	my %firstReadIDs_2_specifiedReadIDs;
	
	my %firstReadID_2_sequences;
	
	my $input_fh;
	open($input_fh, '<', $input_filename) or die "Cannot open $input_filename";
	
	my $firstLine = <$input_fh>;
	die unless($firstLine =~ /^IS/);
	
	while(<$input_fh>)
	{
		my $alignedPair_firstLine = $_;
		chomp($alignedPair_firstLine);
		
		my $alignment_startLine = $.;
		
		die "Cannot parse (first) line $. of file $input_filename\n\n$alignedPair_firstLine\n\n" unless($alignedPair_firstLine =~ /Aligned pair (\d+)/);
		
		my @p1_lines = getNlines($input_fh, 9);
		my @p2_lines = getNlines($input_fh, 9);
		
		die "Cannot parse first line of first read of alignment starting line $alignment_startLine of $input_filename" unless($p1_lines[0] =~ /\s+Read (.+)/); 
		my $read1_ID = $1;
		
		die "Cannot parse first line of second read of alignment starting line $alignment_startLine of $input_filename" unless($p2_lines[0] =~ /\s+Read (.+)/); 
		my $read2_ID = $1;				
		
		if($reads_2_group{$read1_ID} or $reads_2_group{$read2_ID})
		{
			if($reads_2_group{$read1_ID} and $reads_2_group{$read2_ID})  
			{
				die Dumper("Double group for reads $read1_ID and $read2_ID") unless($reads_2_group{$read1_ID} eq $reads_2_group{$read2_ID});
			}
			
			if(exists $reads_2_group{$read1_ID})
			{
				$specifiedRead_2_firstReadID{$read1_ID} = $read1_ID;
				push(@{$firstReadIDs_2_specifiedReadIDs{$read1_ID}}, $read1_ID);
			}
			
			if(exists $reads_2_group{$read2_ID})
			{
				$specifiedRead_2_firstReadID{$read2_ID} = $read1_ID;
				push(@{$firstReadIDs_2_specifiedReadIDs{$read1_ID}}, $read2_ID);
			}
			
			my $group = (exists $reads_2_group{$read1_ID}) ? $reads_2_group{$read1_ID} : $reads_2_group{$read2_ID};
			die unless($group);
			
			push(@{$output_lines_perGroup{$group}}, ($alignment_startLine, @p1_lines, @p2_lines));
			
			die "Cannot parse sequence line of first read of alignment starting line $alignment_startLine of $input_filename" unless($p1_lines[7] =~ /\s+(\w+)/); 
			my $sequence_R1 = $1;

			die "Cannot parse sequence line of second read of alignment starting line $alignment_startLine of $input_filename" unless($p2_lines[7] =~ /\s+(\w+)/); 
			my $sequence_R2 = $1;	
			
			$firstReadID_2_sequences{$read1_ID} = [$sequence_R1, $sequence_R2];
		}
	}
	
	close(INPUT);	
	
	foreach my $foundGroup (keys %output_lines_perGroup)
	{
		my $output_filename = $output_directory . '/' . 'fullAlignments_ ' . $sampleID . '_' . $foundGroup . '.txt';
		open(OUTPUT, '>', $output_filename) or die "Cannot open $output_filename";
		print OUTPUT join("\n", @{$output_lines_perGroup{$foundGroup}}), "\n";
		close(OUTPUT);
	}
	
	my %printed_read1_IDs;
	print SEQUENCEONLY $sampleID, "\n";
	foreach my $groupName (keys %{$interestingReadsBySample{$sampleID}})
	{
		print SEQUENCEONLY "\t", $groupName, "\n";
		
		foreach my $specifiedReadID (@{$interestingReadsBySample{$sampleID}{$groupName}})
		{
			if(exists $specifiedRead_2_firstReadID{$specifiedReadID})
			{
				my $read1_ID = $specifiedRead_2_firstReadID{$specifiedReadID};
				next if($printed_read1_IDs{$read1_ID});
				
				my @specifiedReadIDs = @{$firstReadIDs_2_specifiedReadIDs{$read1_ID}};
				
				foreach my $rID (@specifiedReadIDs)
				{
					print SEQUENCEONLY "\t\t", $rID, "\n";
				}
				
				foreach my $sequence (@{$firstReadID_2_sequences{$read1_ID}})
				{
					print SEQUENCEONLY "\t\t\t", $sequence, "\n";
				}
				
				$printed_read1_IDs{$read1_ID} = 1;
			}
			else
			{
				print SEQUENCEONLY "\t\t", $specifiedReadID, "\n";
				print SEQUENCEONLY "\t\t\t", "Not found", "\n";
			}
		}
	}
	
}

close(SEQUENCEONLY);

sub getNlines
{
	my $input_fh = shift;
	my $n = shift;
	my @lines;
	for(my $i = 1; $i <= $n; $i++)
	{
		my $line = <$input_fh>;
		die unless(defined $line);
		chomp($line);
		push(@lines, $line);
	}
	die unless(scalar(@lines) == $n);
	return @lines;
}	
