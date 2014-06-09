#!/usr/bin/perl -w

use strict;
use 5.010;
use List::MoreUtils qw/all mesh any /;
use Data::Dumper;
use Getopt::Long;   
use Sys::Hostname;
use File::Copy;
use File::Basename;
use Storable;
  
die unless(more_max(1,2,3) == 3);
die unless(more_max(3,2,1) == 3);
die unless(more_max(1,3,2) == 3);

my $kMer_size = 55;  

# input parameters

my $graph;   
my $sample;
my $redo = 0;
my $collect = 0;
my $labelOnly = 1;
my $vcfPos;
my $mapping_display_only_instructions = 0;
my $localOxford = 0;
if(hostname() =~ /(sequoia)|(elm)|(birch)|(banyan)|(cluster3)/)
{
	$localOxford = 1;
}

GetOptions ('graph:s' => \$graph,
 'sample:s' => \$sample, 
 'redo:s' => \$redo,  
 'collect:s' => \$collect, 
 'labelOnly:s' => \$labelOnly,
 'kmer:s' => \$kMer_size,
 'vcfPos:s' => \$vcfPos,
 'mapping_display_only_instructions:s' => \$mapping_display_only_instructions,
 );         

 
# external programs path

my $stampy_bin = qq(/home/dilthey/Stampy/stampy-1.0.20/stampy.py);
my $PICARD_SAM2FASTQ = '/Net/fs1/home/dilthey/picard/picard-tools-1.83/SamToFastq.jar';
my $BWA_bin = qq(/home/dilthey/BWA/bwa-0.6.2/bwa);
my $samtools_bin = qq(/home/dilthey/samtools-0.1.18/samtools);
my $platypus_executable = qq(/home/dilthey/Platypus/Platypus_0.2.0/Platypus.py);

# external data paths

my $path_to_PGF_haplotype = qq(../data/mhc_ref_8_haplotypes/pgf_ref.fasta);
my $reference_genome_FASTA_path = qq(../data/GRCh37.60/fasta);
my $path_to_HumanGenome_graph = qq(../data/GRCh37.60/graph/alex_grc37_auto_X_Y.k${kMer_size}.ctx);

# temporary directory for read re-mapping

my $remapping_basedir = qq(../tmp/readReMapping/);

# data lookup paths

my $sample_path = $localOxford ? qq(/gpfs1/well/gsk_hla/CortexGraphs/) : qq(../data/samples/CortexGraphs);
my $base_path_BAMs = $localOxford ? qq(/gpfs1/well/gsk_hla/bam_output/) : qq(../data/samples/BAMs/);

### do not modify anything below this line

# my $genomic_SNP_information = qq(../data/GS v3 SNP information.txt);
# die "Cannot access $genomic_SNP_information" unless (-e $genomic_SNP_information);

# my $base_path_reads = qq(/Net/x9000/projects/gsk_hla/);


my $xMHC_reference = $path_to_PGF_haplotype;
die "Cannot access $path_to_PGF_haplotype" unless (-e $path_to_PGF_haplotype);

# my $tmpPGF = read_PGF();
# die "PGF spanning from 28702185 to ".(28702185+length($tmpPGF)-1)."\n";
# die read_PGF();

my $memory_height = 26;
my $memory_width = 100;


 
my $pgf_start = 28702185;
my $remapping_flank_length = 100;


 die "Cannot access $path_to_HumanGenome_graph" unless (-e $path_to_HumanGenome_graph);
die unless($kMer_size);

my $expected_normal_graph = $graph.'/graph.txt';
die "Normal graph $expected_normal_graph not there!" unless(-e $expected_normal_graph);

my $expected_kmer_graph_file = $expected_normal_graph.".kmers_".$kMer_size;
die "Augmented kMer graph $expected_kmer_graph_file not there!" unless(-e $expected_kmer_graph_file);

die "Normal graph more recent than kMer graph" unless((stat ($expected_normal_graph))[9] < (stat ($expected_kmer_graph_file))[9]);

my $sample_graph_file = $sample_path.'/'.$sample.'_'.$kMer_size.'.ctx';
#my $sample_graph_glob = $sample_path.'/'.$sample.'_'.$kMer_size.'_*.ctx';
#my @potential_sample_graph_files = glob($sample_graph_glob);
#unless($#potential_sample_graph_files == 0)
#{#
#	die "Failed to find out which binary to use.\nGlob: $sample_graph_glob\n".Dumper(@potential_sample_graph_files);
#}
#my $sample_graph_file = $potential_sample_graph_files[0];

if($sample ne 'N')
{
	die "Sample graph $sample_graph_file not there!" unless(-e $sample_graph_file)
}

# Load list of required kMers

my $expected_required_kMers_file = $expected_kmer_graph_file.".requiredKMers";
unless((-e $expected_required_kMers_file) and (!$redo) and ((stat ($expected_kmer_graph_file))[9] < (stat ($expected_required_kMers_file))[9]))
{
	print "Determining required kMers\n";
	my $cmd = qq(../bin/MHC-PRG domode determineRequiredKMers $expected_kmer_graph_file $expected_required_kMers_file);
	my $output = `$cmd`;    
	unless(-e $expected_required_kMers_file)
	{
		die "Could not determine required kMers.\n\nCommand:\n$cmd\n\nOutput:\n\n$output";
	}
	
	open(KMERS, "<", $expected_required_kMers_file) or die "Cannot open $expected_required_kMers_file";
	while(<KMERS>)
	{
		my $line = $_;
		chomp($line);
		last if ($. > 5);
		die "Wrong length of kMers in graph $expected_kmer_graph_file: expect $kMer_size, got ".length($line) unless(length($line) == $kMer_size);
	}
	close(KMERS);
}

# Count occurences of required kMers in the reference genome (cortex-based, expect output file)

my %kMer_reference_count;
my $kMer_count_reference = $expected_required_kMers_file.'.binaryCount';
unless((-e $kMer_count_reference) and (!$redo))
{
	print "Count reference kMers\n";
	#my $cmd = qq($cortex_binary --multicolour_bin $path_to_HumanGenome_graph  --kmer_size $kMer_size --mem_height ${memory_height} --mem_width ${memory_width} --dilthey2 $expected_required_kMers_file);
	my $cmd = qq(./readCortexCoverage.pl --cortex_bin $path_to_HumanGenome_graph  --kMer_size $kMer_size --interesting_kMers $expected_required_kMers_file);
	my $output = `$cmd`;
	unless(-e $kMer_count_reference)
	{
		die "Could not count reference kMers\n\nCommand:\n $cmd \n\nOutput:\n\n$output";
	}
}

open(KMERS, "<", $kMer_count_reference) or die "Cannot open $kMer_count_reference";
while(<KMERS>)
{#
	my $line = $_;
	chomp($line);
	next unless($line);
	next if($line =~ /^Mean/);
	next if($line =~ /^Total/);
	
	my @f = split(/ /, $line);
	unless($#f == 1)
	{
		die "Strange line in $kMer_count_reference -- expect two fields, separated by a space.\n\nL:$line";
	}
	die "length ".length($f[0])." vs expected $kMer_size in file $kMer_count_reference" unless(length($f[0])==$kMer_size);
	$kMer_reference_count{$f[0]} = $f[1];
}
close(KMERS);

# Count occurences of required kMers in the sample genome (cortex-based, expect output file)

my $graphForFileName = $graph;
$graphForFileName =~ s/^.+tmp2//;
$graphForFileName =~ s/\W/_/g;
my $kMer_count_sample_required = qq(../tmp/kMerCount_).join('_', $graphForFileName, $sample, $kMer_size,  'required');
my $kMer_count_sample = $kMer_count_sample_required.'.binaryCount';
if($sample ne 'N')
{
	unless((-e $kMer_count_sample) and (!$redo))
	{
		print "Count sample kMers\n";
		copy($expected_required_kMers_file, $kMer_count_sample_required) or die "Cannot copy $expected_required_kMers_file to $kMer_count_sample_required";	
		#my $cmd = qq($cortex_binary --multicolour_bin $sample_graph_file  --kmer_size $kMer_size --mem_height ${memory_height} --mem_width ${memory_width} --dilthey2 $kMer_count_sample_required);
		my $cmd = qq(./readCortexCoverage.pl --cortex_bin $sample_graph_file  --kMer_size $kMer_size --interesting_kMers $kMer_count_sample_required);
		
		my $output = `$cmd`;
		unless(-e $kMer_count_sample)
		{
			die "Could not count sample kMers\n\nCommand:\n$cmd\n\nOutput:\n\n$output";
		}
	}
}
else
{
	my @kMers;
	open(KMERSREQ, "<", $expected_required_kMers_file) or die "Cannot open $expected_required_kMers_file";
	while(<KMERSREQ>)
	{
		my $line = $_;
		chomp($line);
		push(@kMers, $line);
	}
	close(KMERSREQ);
		
	open(KMERS, '>', $kMer_count_sample) or die "Cannot open $kMer_count_sample";
	foreach my $kMer(@kMers)
	{
		print KMERS $kMer, " ", 1, "\n";
	}
	print KMERS 'MeanReadLen', ' ', 100, "\n";
	print KMERS 'TotalKMerCoverage', ' ', 3000000000*30 * (55/100), "\n";
	print KMERS 'TotalSeq', ' ', 3000000000*30, "\n";	
	
	close(KMERS);
}

# Correct the reference genome-based kMer count. How many of these kMers are  in the area that is covered by the haplotype
# graph? This number needs to be subtracted.

my $output_corrected_kMers_reference = $expected_kmer_graph_file.".kMersReferenceCorrected";

if((not -e $output_corrected_kMers_reference) or ($redo))
{
	# determine SNPs in graph file
	my %loci_in_graph;
	open(GRAPH, "<", $expected_normal_graph) or die "Cannot open $expected_normal_graph\n";
	my $firstLine = <GRAPH>;
	chomp($firstLine);
	unless($firstLine eq 'CODE:')
	{
		die "Expect first line to be CODE: - graph format changed, file $expected_normal_graph?";
	}
	while(<GRAPH>)
	{
		my $line = $_;
		chomp($line);
		if($line =~ /\:$/)
		{
			last;
		}
		my @fields = split(/\|\|\|/, $line);
		die unless($#fields == 2);
		my $locusID = $fields[0];
		$loci_in_graph{$locusID} = 1;
	}
	close(GRAPH);
	
	my @positions_in_graph = grep {$_ =~ /^S\d_(\d+)_(\d+)$/} keys %loci_in_graph;
	@positions_in_graph = map {$_ =~ /^S\d_(\d+)_(\d+)$/; $2} @positions_in_graph;
	
	unless(scalar(@positions_in_graph) > 1)
	{
		die "Have less than two loci with integrated coordinates - chr6:1234... - in the graph. Cannot determine graph boundaries";
	}
	
	# find the positions of the first and the last SNP in the graph on reference coordinates
	
	@positions_in_graph = sort {$a <=> $b} @positions_in_graph;
	
	my $first_position_graph_inPGF = $positions_in_graph[0];
	my $last_position_graph_inPGF = $positions_in_graph[$#positions_in_graph];
	
	die unless($first_position_graph_inPGF < $last_position_graph_inPGF);
	
	my %subtract_kMers;
	print "Remove kMers for positions $first_position_graph_inPGF - $last_position_graph_inPGF (in PGF)\n";
	
	# read the pgf (reference) haplotype, extract the stretch spanned by first and last SNP in graph,
	# count kMers, subtract from reference genome kMer count
	
	my $pgf_haplotype = read_PGF();
	(length($pgf_haplotype) > $first_position_graph_inPGF+$kMer_size-1) or die Dumper("PGF coordinate problem!", length($pgf_haplotype), $first_position_graph_inPGF, $kMer_size)."\n\n";
	
	my $kMer = substr($pgf_haplotype, $first_position_graph_inPGF, $kMer_size-1);
	(length($kMer) == ($kMer_size-1)) or die "Problem with extracting position $first_position_graph_inPGF + $kMer_size symbols, first position in graph $first_position_graph_inPGF , last position in graph $last_position_graph_inPGF , total PGF length ".length($pgf_haplotype);
	
	for(my $end_KMer = $first_position_graph_inPGF+$kMer_size-1; $end_KMer <= $last_position_graph_inPGF; $end_KMer++)
	{
		my $char_PGF = substr($pgf_haplotype, $end_KMer, 1);
		(length($char_PGF) == 1) or die;
		$kMer = $kMer.$char_PGF;
		(length($kMer) == $kMer_size) or die;
		$subtract_kMers{$kMer}++;
		substr($kMer, 0, 1) = '';       
	}
	
	# die Dumper(\%subtract_kMers);
	
	open(OUTPUT, '>', $output_corrected_kMers_reference) or die "Cannot open $output_corrected_kMers_reference";
	my $have_subtracted = 0;
	foreach my $kMer (sort keys %kMer_reference_count)
	{
		my $original_count = $kMer_reference_count{$kMer};
		my $subtract_count = (exists $subtract_kMers{$kMer}) ? $subtract_kMers{$kMer} : 0;
		my $new_count = $original_count - $subtract_count;
		if($subtract_count > 0)
		{
			$have_subtracted = 1;
		}
		die unless($new_count >= 0);
		print OUTPUT join(' ', $kMer, $new_count), "\n";
	}
	close(OUTPUT);
	
	unless($have_subtracted)
	{
		unlink $output_corrected_kMers_reference;		
		die "No PGF kMer was subtracted to correct the reference counts. Bad. Fix, and --redo!\n";
	}
}

my $label_part = ($labelOnly) ? ('--labelonly') : '';
my $final_command = qq(../bin/MHC-PRG domode nextGenInference $expected_kmer_graph_file $output_corrected_kMers_reference $kMer_count_sample $label_part --genotypingMode 8);
print "Execute:\n", $final_command, "\n\n";

if($collect eq '3')
{
	my $needGZ = 0;
	my $expected_output_haplotypes = $kMer_count_sample.'.viterbiHaplotypes';
	
	# get graph information
	
	my @loci_in_graph;
	my $graph_segments_file = $graph.'/segments.txt';
	die "File not there: $graph_segments_file" unless (-e $graph_segments_file);
	open(SEGMENTS, '<', $graph_segments_file) or die;
	while(<SEGMENTS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
		next unless($line);
		
		my $segment_file = $graph.'/'.$line;

		open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
		my $firstLine = <SEGMENT>;
		chomp($firstLine);
		$firstLine =~ s/\n//g;
		$firstLine =~ s/\r//g;			
		my @line_fields = split(/ /, $firstLine);
		shift(@line_fields); # kick out individual ID
		push(@loci_in_graph, @line_fields);
		close(SEGMENT);
	}
	close(SEGMENTS);
	
	my %locus_position;
	my $positions_file = $graph.'/positions.txt';
	die "File not there: $positions_file" unless (-e $positions_file);
	open(POSITIONS, '<', $positions_file) or die;	
	while(<POSITIONS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;		
		my @fields = split(/ /, $line);
		my $field_id = $fields[0];
		my $field_position = $fields[1];
		$locus_position{$field_id} = $field_position;
	}
	close(POSITIONS);
	
	for(my $lI = 1; $lI <= $#loci_in_graph; $lI++)
	{
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI]", $lI);
		}
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI-1]", $lI-1);
		}		
		die unless($locus_position{$loci_in_graph[$lI]} > $locus_position{$loci_in_graph[$lI-1]});
	}
	
	# read xMHC reference
	
	open(REF, "<", $xMHC_reference) or die "Cannot open $xMHC_reference";
	my $pgf_ref_seq;
	while(<REF>)
	{
		my $line = $_; 
		chomp($line);
		$line =~ s/[\x0A\x0D]+//g;	
		next if ($line =~ /^\>/);
		$pgf_ref_seq .= $line;
	}
	close(REF);
	my $pgf_stop = $pgf_start+length($pgf_ref_seq);
#die $pgf_stop;
	
	print "Reading graph information completed\n";
	
	# get the first two haplotypes
	
	open(HAPLOTYPES, '<', $expected_output_haplotypes) or die "Cannot open $expected_output_haplotypes";
	my $haplotype_1, my $haplotype_2;
	my @haplotype_1_alignment_fields;
	my @haplotype_2_alignment_fields;
	while(<HAPLOTYPES>)
	{
		my $line = $_;
		chomp($line);
		last if ($. > 2);
		my @fields = split(/ /, $line);
		my $haplotype = join('', @fields[2 .. ($#fields-2)]);
		if($. == 1)
		{
			$haplotype_1 = $haplotype;
			@haplotype_1_alignment_fields = @fields[2 .. ($#fields-2)];
		}
		elsif($. == 2)
		{
			$haplotype_2 = $haplotype;	
			@haplotype_2_alignment_fields = @fields[2 .. ($#fields-2)];
		}
	}
	close(HAPLOTYPES);
	
	$| = 1;
	
	die unless($haplotype_1 and $haplotype_2);
	die unless($#haplotype_1_alignment_fields == $#loci_in_graph);
	die unless($#haplotype_2_alignment_fields == $#loci_in_graph);
		
	print "Reading two inferred haplotypes completed\n";
	
	my @haplotype_1_alignment_index;
	my @haplotype_2_alignment_index;
	for(my $i = 0; $i < length($haplotype_1); $i++)
	{
		my $character = substr($haplotype_1, $i, 1);
		if($character ne '_')
		{
			push(@haplotype_1_alignment_index, $i);
		}
	}
	for(my $i = 0; $i < length($haplotype_2); $i++)
	{
		my $character = substr($haplotype_2, $i, 1);
		if($character ne '_')
		{
			push(@haplotype_2_alignment_index, $i);
		}
	}	
	
	# print join(', ', @haplotype_1_alignment_index[0 .. 30]), "\n";
	# print join(', ', @haplotype_2_alignment_index[0 .. 30]), "\n";	
	# exit;

	$haplotype_1 =~ s/\*/N/g;
	$haplotype_1 =~ s/_//g;
	unless($haplotype_1 =~ /^[ACGTNacgtn]+$/)
	{
		my $nonOK = $haplotype_1;
		$nonOK =~ s/[ACGTNacgtn]//g;
		die Dumper('h1', length($nonOK), $nonOK);
	}
	
	$haplotype_2 =~ s/\*/N/g;
	$haplotype_2 =~ s/_//g;
	unless($haplotype_2 =~ /^[ACGTNacgtn]+$/)
	{
		my $nonOK = $haplotype_2;
		$nonOK =~ s/[ACGTNacgtn]//g;
		die Dumper('h2', length($nonOK), $nonOK);
	}
	
	
	die Dumper(length($haplotype_1), scalar(@haplotype_1_alignment_index)) unless(length($haplotype_1) == scalar(@haplotype_1_alignment_index));
	die Dumper(length($haplotype_2), scalar(@haplotype_2_alignment_index)) unless(length($haplotype_2) == scalar(@haplotype_2_alignment_index));
	

	# mapping directory
	my $remapping_dir = $remapping_basedir.'/kMerCount_'.join('_', $graphForFileName, $sample, $kMer_size,  'remapping');
	my $xMHC_VCF_1_readDecision = $remapping_dir.'/mapped_xMHC_1.vcf';
	my $xMHC_VCF_2_readDecision = $remapping_dir.'/mapped_xMHC_2.vcf';
	unless((-e $xMHC_VCF_1_readDecision) and (-e $xMHC_VCF_2_readDecision))
	{
		my $remapping_dir_rawGenome_1 = $remapping_dir."/genomeRef_1";
		my $remapping_dir_rawGenome_2 = $remapping_dir."/genomeRef_2";
		unless(-e $remapping_dir)
		{
			mkdir($remapping_dir) or die "Cannot create $remapping_dir"; 	
		}
			
		
		unless(-e $remapping_dir_rawGenome_1)
		{
			mkdir($remapping_dir_rawGenome_1) or die "Cannot create $remapping_dir_rawGenome_1";
			mkdir($remapping_dir_rawGenome_2) or die "Cannot create $remapping_dir_rawGenome_2";
			
			# copy reference genome
			my @reference_FASTAs = glob($reference_genome_FASTA_path.'/*.fa');
			die "Cannot find reference FASTA files - did you specify the right directory? \[ $reference_genome_FASTA_path \]" unless(@reference_FASTAs);
			foreach my $fF (@reference_FASTAs)
			{
				my $basename = fileparse($fF);
				my $newfilename_1 = $remapping_dir_rawGenome_1.'/'.$basename;
				my $newfilename_2 = $remapping_dir_rawGenome_2.'/'.$basename;			
				print "Copying $fF ... \n";
				copy($fF, $newfilename_1) or die "Cannot copy file $fF to $newfilename_1";
				copy($fF, $newfilename_2) or die "Cannot copy file $fF to $newfilename_2";			
			}
		}
		

		my $remapping_dir_reads = $remapping_dir."/reads/";
		unless(-e $remapping_dir_reads )
		{
			mkdir($remapping_dir_reads) or die;
		}
		my $BAM_file;

		if($localOxford and ($sample =~ /Z[12]/))
		{
			$remapping_dir_reads = '/Net/birch/data/oa/NA12878_ILLUMINA_PLATINUM/';
			$needGZ = 1;
			unless($mapping_display_only_instructions)
			{
				die "We want to deal with the reads of NA12878, but we can do so only in external mapping!\n";
			}
		}	
		else
		{
			$BAM_file = find_individual_BAM($sample);
			my $extraction_dir_OK_flag = $remapping_dir_reads.'/'.'extraction_successful';

			unless(-e $extraction_dir_OK_flag)
			{	
				unless(-e $BAM_file)
				{
					die "BAM file $BAM_file not there, but required!";
				}
				
				my $cmd = qq(java -Xmx2g -jar ${PICARD_SAM2FASTQ} INPUT=${BAM_file} OUTPUT_PER_RG=TRUE OUTPUT_DIR=${remapping_dir_reads});
				print `$cmd`;
				
				open(OK, '>', $extraction_dir_OK_flag) or die "Cannot open $extraction_dir_OK_flag";
				print OK 1;
				close(OK);
				
				print "Reads extracted\n";		
			}
		}
		
		# my $extraction_dir_OK_flag = $remapping_dir_reads.'/'.'unzipping_successful';

		# unless(-e $extraction_dir_OK_flag)
		# {
			# Get and extract reads
			
			# my $directory = &find_individual_readsPath($sample);
			# die "Directory not existing: $directory\n" unless (-e $directory);
			# unless(-e $remapping_dir_reads)
			# {
				# mkdir($remapping_dir_reads) or die "Cannot mkdir $remapping_dir_reads";
			# }
			
			# my @zipped_original_files = glob($directory.'/*.gz');
			# my @zipped_temp_files;
			# foreach my $file (@zipped_original_files)
			# {
				# my ($filename,$filepath) = fileparse($file);
				# my $zip_new_temp_name = $remapping_dir_reads.'/'.$filename;
				# my $copy_try = 1;
				# while($copy_try != 0)
				# {
					# my $success = copy($file, $zip_new_temp_name);
					# if($success)
					# {
						# $copy_try = 0;
						# next;
					# }
					# $copy_try++;
					# if($copy_try > 5)
					# {
						# die "Fifth attempt to copy $file to $zip_new_temp_name has failed!\n";
					# }
				# }
				# my $unzip_command = qq(gunzip $zip_new_temp_name);
				# print "Executing $unzip_command\n";
				# system $unzip_command;
				# push(@zipped_temp_files, $zip_new_temp_name);
			# }
			
			# open(OK, '>', $extraction_dir_OK_flag) or die "Cannot open $extraction_dir_OK_flag";
			# print OK 1;
			# close(OK);
			
			# print "Reads extracted\n";
		# }

		my @chromosomal_files_1 = (glob($remapping_dir_rawGenome_1.'/*.fa'));
		my @chromosomal_files_2 = (glob($remapping_dir_rawGenome_2.'/*.fa'));	
		my $chromosome_6_file_1;
		my $chromosome_6_file_2;	
		foreach my $f (@chromosomal_files_1)
		{
			open(CHROMO, '<', $f) or die "Cannot open $f";
			my $firstLine = <CHROMO>;
			close(CHROMO);
			if(($firstLine =~ /^>chr6/) or ($firstLine =~ />6/))
			{
				$chromosome_6_file_1 = $f;
			}
		}
		foreach my $f (@chromosomal_files_2)
		{
			open(CHROMO, '<', $f) or die "Cannot open $f";
			my $firstLine = <CHROMO>;
			close(CHROMO);
			if(($firstLine =~ /^>chr6/) or ($firstLine =~ />6/))
			{
				$chromosome_6_file_2 = $f;
			}
		}	
		
	
		die "Cannot find file for chromosome 6 [1]"  unless ($chromosome_6_file_1);
		die "Cannot find file for chromosome 6 [2]" unless ($chromosome_6_file_2);	
		
		foreach my $files ([$chromosome_6_file_1, $remapping_dir_rawGenome_1, $haplotype_1], [$chromosome_6_file_2, $remapping_dir_rawGenome_2, $haplotype_2])
		{
			my $chromosome_6_file = $files->[0];
			my $remapping_dir_rawGenome = $files->[1];
			my $haplotype = $files->[2];
					
			my $chromo_6_header;
			my $chromo_6_content;
					
			my $chromo6;
			open($chromo6, '<', $chromosome_6_file) or die "Cannot open $chromosome_6_file";
						
			$chromo_6_header = <$chromo6>;
			chomp($chromo_6_header);	
					
			$chromo_6_header =~ s/\n//g;
			$chromo_6_header =~ s/\r//g;
					
			if($chromo_6_header !~ /xMHC_excised/)
			{		
				while(<$chromo6>)
				{
					my $line = $_;
					chomp($line);	
					die if($line =~ />/);
					$line =~ s/\n//g;
					$line =~ s/\r//g;
					$chromo_6_content .= $line;
				}
				
				# now get positions to blind
				
				my %loci_in_graph;
				open(GRAPH, "<", $expected_normal_graph) or die "Cannot open $expected_normal_graph\n";
				my $firstLine = <GRAPH>;
				chomp($firstLine);
				unless($firstLine eq 'CODE:')
				{
					die "Expect first line to be CODE: - graph format changed, file $expected_normal_graph?";
				}
				while(<GRAPH>)
				{
					my $line = $_;
					chomp($line);
					if($line =~ /\:$/)
					{
						last;
					}
					my @fields = split(/\|\|\|/, $line);
					die unless($#fields == 2);
					my $locusID = $fields[0];
					$loci_in_graph{$locusID} = 1;
				}
				close(GRAPH);
				
				my @positions_in_graph = grep {$_ =~ /^S\d_(\d+)_(\d+)$/} keys %loci_in_graph;
				@positions_in_graph = map {$_ =~ /^S\d_(\d+)_(\d+)$/; $2} @positions_in_graph;
				
				unless(scalar(@positions_in_graph) > 1)
				{
					die "Have less than two loci with integrated coordinates - chr6:1234... - in the graph. Cannot determine graph boundaries";
				}
				
				# find the positions of the first and the last SNP in the graph on reference coordinates
				
				@positions_in_graph = sort {$a <=> $b} @positions_in_graph;
				
				my $first_position_graph_inPGF = $positions_in_graph[0];
				my $last_position_graph_inPGF = $positions_in_graph[$#positions_in_graph];
				
				die unless($first_position_graph_inPGF < $last_position_graph_inPGF);
							
				my $pgf_sequence = read_PGF();
				my $pgf_sequence_coveredByGraph = substr($pgf_sequence, $first_position_graph_inPGF, $last_position_graph_inPGF - $first_position_graph_inPGF + 1);
				my $expected_sequence_in_genome = substr($chromo_6_content, $pgf_start + $first_position_graph_inPGF - 1, $last_position_graph_inPGF - $first_position_graph_inPGF + 1);
				
				# die "PGF does not match up with reference, this is not good!" unless ($pgf_sequence_coveredByGraph eq $expected_sequence_in_genome);
				
				my $lengthN = 'N' x length($pgf_sequence_coveredByGraph);
				die unless(length($lengthN) == length($expected_sequence_in_genome));

				my $genomic_graph_start = $pgf_start + $first_position_graph_inPGF;
				my $genomic_graph_lastPosition = $pgf_start + $first_position_graph_inPGF - 1 + length($lengthN) - 1;
				my $graph_before_flank = substr($chromo_6_content, $genomic_graph_start - 1 - $remapping_flank_length, $remapping_flank_length);
				my $graph_after_flank = substr($chromo_6_content, $genomic_graph_lastPosition + 1, $remapping_flank_length);		
				my $combinedBeforeAfter = $graph_before_flank . $expected_sequence_in_genome . $graph_after_flank;
				my $genomeSurrounding = substr($chromo_6_content, $pgf_start + $first_position_graph_inPGF - 1 - 1000, $last_position_graph_inPGF - $first_position_graph_inPGF + 1 + 2000);
				
				die if(index($genomeSurrounding, $combinedBeforeAfter) == -1);

				substr($chromo_6_content, $pgf_start + $first_position_graph_inPGF - 1, $last_position_graph_inPGF - $first_position_graph_inPGF + 1) = $lengthN;
				
				close($chromo6);
				
				open($chromo6, '>', $chromosome_6_file) or die "Cannot open $chromosome_6_file";
				print { $chromo6 } $chromo_6_header.' xMHC_excised'."\n";
				
				while($chromo_6_content ne "")
				{
					my $consumeLength = 60;
					if($consumeLength > length($chromo_6_content))
					{
						$consumeLength = length($chromo_6_content);
					}
					print { $chromo6 } substr($chromo_6_content, 0, $consumeLength), "\n";
					substr($chromo_6_content, 0, $consumeLength) = '';
				}
				
				my $file_haplo1 = $remapping_dir_rawGenome.'/xMHC_haplo1.fa';
					
				open(CHR1, '>', $file_haplo1) or die "Cannot open $file_haplo1";
							
				print CHR1 ">xMHCpseudoChromHaplotype chr6:${genomic_graph_start}:${genomic_graph_lastPosition}", "\n";
				
				my $written_length = 0;
				
				my $haplotype_to_print = $graph_before_flank.$haplotype.$graph_after_flank;
				
				while($haplotype_to_print ne "")
				{
					my $consumeLength = 60;
					if($consumeLength > length($haplotype_to_print))
					{
						$consumeLength = length($haplotype_to_print);
					}
					my $toWrite = substr($haplotype_to_print, 0, $consumeLength);
					$written_length += length($toWrite);
					print CHR1 $toWrite, "\n";
					substr($haplotype_to_print, 0, $consumeLength) = '';
				}
					
				close(CHR1);
				
			}
			close($chromo6);	
		}
		

		# pseudo genome references
		
		my $xMHC_IGV_pseudoGenome_1 = $remapping_dir.'/IGV_pseudoGenome_1.fa';
		my $xMHC_IGV_pseudoGenome_2 = $remapping_dir.'/IGV_pseudoGenome_2.fa';

		copy($remapping_dir_rawGenome_1.'/xMHC_haplo1.fa', $remapping_dir.'/IGV_pseudoGenome_1.fa');
		copy($remapping_dir_rawGenome_2.'/xMHC_haplo1.fa', $remapping_dir.'/IGV_pseudoGenome_2.fa');
		
		my $cmd = qq(${samtools_bin} faidx ${xMHC_IGV_pseudoGenome_1});
		print $cmd, "\n", `$cmd`, "\n";	
		$cmd = qq(${samtools_bin} faidx ${xMHC_IGV_pseudoGenome_2});
		print $cmd, "\n", `$cmd`, "\n";	
		
		print "Chromosome 6 and pseudo-chromosomes prepared!\n";
		
		my $combined_genome_file_f1 = $remapping_dir.'/combinedGenome_1.fasta';
		my $combined_genome_file_f2 = $remapping_dir.'/combinedGenome_2.fasta';;

		my $igv1_genome;
		open(IGV1, '<', $remapping_dir.'/IGV_pseudoGenome_1.fa') or die;
		while(<IGV1>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			next if(substr($line, 0, 1) eq '>');
			$igv1_genome .= $line;
		}
		close(IGV1);
		$igv1_genome = substr($igv1_genome, $remapping_flank_length);
		$igv1_genome = substr($igv1_genome, 0, length($igv1_genome) - $remapping_flank_length);

		my $igv2_genome;
		open(IGV2, '<', $remapping_dir.'/IGV_pseudoGenome_2.fa') or die;
		while(<IGV2>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			next if(substr($line, 0, 1) eq '>');
			$igv2_genome .= $line;
		}
		close(IGV2);
		$igv2_genome = substr($igv2_genome, $remapping_flank_length);
		$igv2_genome = substr($igv2_genome, 0, length($igv2_genome) - $remapping_flank_length);
		
		my $expected_igv1_genome = join('', @haplotype_1_alignment_fields);
		$expected_igv1_genome =~ s/\*/N/g;
		$expected_igv1_genome =~ s/_//g;

		my $expected_igv2_genome = join('', @haplotype_2_alignment_fields);
		$expected_igv2_genome =~ s/\*/N/g;   
		$expected_igv2_genome =~ s/_//g;

		die unless($expected_igv1_genome eq $igv1_genome);
		die unless($expected_igv2_genome eq $igv2_genome);

		# print Dumper("Position 4634404", substr($expected_igv2_genome, 4634404, 1), substr($expected_igv2_genome, 4634404-15, 31)), "\n";

		if($mapping_display_only_instructions)
		{
			my $fa_idx_2 = $combined_genome_file_f2.'.fai';   
			unless(-e $fa_idx_2)
			{
				my @chromosomal_files_1 = (glob($remapping_dir_rawGenome_1.'/*.fa'));
				my @chromosomal_files_2 = (glob($remapping_dir_rawGenome_2.'/*.fa'));	
					
				open(COMBINED1, '>', $combined_genome_file_f1) or die "Cannot open $combined_genome_file_f1";
				open(COMBINED2, '>', $combined_genome_file_f2) or die "Cannot open $combined_genome_file_f2";	
				foreach my $f (@chromosomal_files_1)
				{
					# next unless($f =~ /xMHC/); # remove later
					
					open(F, '<', $f) or die "Cannot open $f";
					while(<F>)
					{
						print COMBINED1 $_;
					}
					close(F);
				} 
			
				foreach my $f (@chromosomal_files_2)
				{
					# next unless($f =~ /xMHC/); # remove later
				
					open(F, '<', $f) or die "Cannot open $f";
					while(<F>)
					{
						print COMBINED2 $_;
					}
					close(F);
				}	
				close(COMBINED1);
				close(COMBINED2);	
				
				my $cmd = qq(${samtools_bin} faidx ${combined_genome_file_f1});
				print $cmd, "\n", `$cmd`, "\n";	
				$cmd = qq(${samtools_bin} faidx ${combined_genome_file_f2});
				print $cmd, "\n", `$cmd`, "\n";	
			}
			die "No fasta index file $fa_idx_2 ??" unless(-e $fa_idx_2);
		}
		else
		{
			unless($mapping_display_only_instructions or ((-e $combined_genome_file_f1.'.sa') and (-e $combined_genome_file_f2.'.sa')))
			{	
			
				my @chromosomal_files_1 = (glob($remapping_dir_rawGenome_1.'/*.fa'));
				my @chromosomal_files_2 = (glob($remapping_dir_rawGenome_2.'/*.fa'));	
					
				open(COMBINED1, '>', $combined_genome_file_f1) or die "Cannot open $combined_genome_file_f1";
				open(COMBINED2, '>', $combined_genome_file_f2) or die "Cannot open $combined_genome_file_f2";	
				foreach my $f (@chromosomal_files_1)
				{
					# next unless($f =~ /xMHC/); # remove later
					
					open(F, '<', $f) or die "Cannot open $f";
					while(<F>)
					{
						print COMBINED1 $_;
					}
					close(F);
				} 
			
				foreach my $f (@chromosomal_files_2)
				{
					# next unless($f =~ /xMHC/); # remove later
				
					open(F, '<', $f) or die "Cannot open $f";
					while(<F>)
					{
						print COMBINED2 $_;
					}
					close(F);
				}	
				close(COMBINED1);
				close(COMBINED2);	
				
				my $cmd = qq(${samtools_bin} faidx ${combined_genome_file_f1});
				print $cmd, "\n", `$cmd`, "\n";	
				$cmd = qq(${samtools_bin} faidx ${combined_genome_file_f2});
				print $cmd, "\n", `$cmd`, "\n";	
				
				$cmd = qq(${BWA_bin} index -a bwtsw $combined_genome_file_f1);
				print $cmd, "\n\n";
				print `$cmd`, "\n\n";

				$cmd = qq(${BWA_bin} index -a bwtsw $combined_genome_file_f2);
				print $cmd, "\n\n";
				print `$cmd`, "\n\n";
				
				print "Indices for BWA built\n";
			}
		}
		
		my $temp_alignment_1 = $remapping_dir.'/temp_alignment_1';
		my $temp_alignment_2 = $remapping_dir.'/temp_alignment_2';
		unless(-e $temp_alignment_1)
		{
			mkdir($temp_alignment_1) or die;
		}
		unless(-e $temp_alignment_2)
		{
			mkdir($temp_alignment_2) or die;
		}
		
		my $merged_BAM_1 = $remapping_dir.'/merged_1.bam';
		my $merged_BAM_2 = $remapping_dir.'/merged_2.bam';
		
		if($mapping_display_only_instructions)
		{
			my $sorted_bam_1 = "${merged_BAM_1}.sorted.bam";
			my $sorted_bam_2 = "${merged_BAM_2}.sorted.bam";
			if((-e $sorted_bam_1) and (-e $sorted_bam_2))
			{
				print "Found $sorted_bam_1 and $sorted_bam_2 , continue!\n";
			}
			else
			{
				unless($localOxford)
				{ 
					die "You activated the switch \$mapping_display_only_instructions, but the corresponding code makes only sense in our local Oxford environment - track down this error message in nextGenInferenceVariGraph.pl, and modify the code following that line according to your environment.";
				}	
				print "\n\nYou now need to make sure $sorted_bam_1 and $sorted_bam_2 are there, by mapping\n\n";
				print "READS: $remapping_dir_reads \n\n";
				print "REFERENCES: $remapping_dir_rawGenome_1 and $remapping_dir_rawGenome_2\nn";
				print "PAIRED-END: 0\n\n";
				print "Command suggestion:\ncd /gpfs1/well/gsk_hla/readLengthSimulator/src; ./alignGenome.pl --action prepare --reference_genome_FASTA_path $remapping_dir_rawGenome_1 --paired_end 0 --input_directory $remapping_dir_reads --name_for_project REMAPPING_${sample}_1 --gz $needGZ; ./alignGenome.pl --action prepare --reference_genome_FASTA_path $remapping_dir_rawGenome_2 --paired_end 0 --input_directory $remapping_dir_reads --name_for_project REMAPPING_${sample}_2 --gz $needGZ\n\n\n";
				print "\n\nFinally, copy BAMs to\n${sorted_bam_1}\n${sorted_bam_2}\n";
				exit;
			}
		} 
		else
		{
			unless((-e $merged_BAM_1) and (-e $merged_BAM_2))
			{
			
				die if($needGZ);
				
				my @fqfiles = glob($remapping_dir_reads.'/*.fastq');
					
				my %pe_pairs = map {$_ => 1} map {($_ =~ /.+\/(.+)_(1|2).+?$/) or die "Cannot parse filename $_"; $1} @fqfiles;
				my %pe_counts;
				my @pe_lists = ([], []);
				my @se_lists = ();

				foreach my $pe_pair (keys %pe_pairs)
				{
					foreach my $d (1, 2)
					{			
						my $f = $remapping_dir_reads.'/'.$pe_pair.'_'.$d.'.fastq';
						next unless (-e $f);
						push(@{$pe_counts{$pe_pair}}, $f);  
					}
				}
					
				my $found_pe = 0;
				foreach my $present_two_key (grep {scalar(@{$pe_counts{$_}}) == 2} keys %pe_counts)
				{
					$found_pe++;
					push(@{$pe_lists[0]}, $pe_counts{$present_two_key}[0]);
					push(@{$pe_lists[1]}, $pe_counts{$present_two_key}[1]);
				}
				
				foreach my $present_one_key (grep {scalar(@{$pe_counts{$_}}) == 1} keys %pe_counts)
				{
					warn "There is a file without mated reads: $present_one_key -- not good, but ignore and resume.\n";
					push(@se_lists, $pe_counts{$present_one_key}[0]);
				}
				
				my @BAMs_toMerge_1;
				my @BAMs_toMerge_2;
				
				die "No files to merge?" unless (scalar(@se_lists) or scalar(@{$pe_lists[0]}));
				
				for(my $i = 0; $i < scalar(@se_lists); $i++)
				{
					print "Align reads chunk $i (no pairs) using BWA [1]\n";
					
					my $reads_file = $se_lists[$i];
					my $sai_file_1 = $temp_alignment_1.'/'.fileparse($reads_file).'.sai';
					my $sai_file_2 = $temp_alignment_2.'/'.fileparse($reads_file).'.sai';
					my $bam_file_1 = $sai_file_1.'.bam';
					my $bam_file_2 = $sai_file_2.'.bam';
					
					my $cmd = qq(${BWA_bin} aln -q10 -t42 ${combined_genome_file_f1} ${reads_file} > ${sai_file_1});
					print $cmd, "\n", `$cmd`, "\n";
					$cmd = qq(${BWA_bin} samse ${combined_genome_file_f1} ${sai_file_1} ${reads_file} | ${samtools_bin} view -Sb - > ${bam_file_1});
					print $cmd, "\n", `$cmd`, "\n";
					
					$cmd = qq(${BWA_bin} aln -q10 -t42 ${combined_genome_file_f2} ${reads_file} > ${sai_file_2});
					print $cmd, "\n", `$cmd`, "\n";
					$cmd = qq(${BWA_bin} samse ${combined_genome_file_f2} ${sai_file_2} ${reads_file} | ${samtools_bin} view -Sb - > ${bam_file_2});
					print $cmd, "\n", `$cmd`, "\n";
					
					# exit;
					
					push(@BAMs_toMerge_1, $bam_file_1);
					push(@BAMs_toMerge_2, $bam_file_2);
					
					# last; # remove later
				}
				
				for(my $i = 0; $i < scalar(@{$pe_lists[0]}); $i++)
				{
					# last; # remove later
					
					print "Align reads chunk $i (paired reads) using BWA [1]\n";
					
					foreach my $reads_file ($pe_lists[0][$i], $pe_lists[1][$i])
					{
						my $sai_file_1 = $temp_alignment_1.'/'.fileparse($reads_file);
						my $sai_file_2 = $temp_alignment_2.'/'.fileparse($reads_file);
						my $bam_file_1 = $sai_file_1.'.bam';
						my $bam_file_2 = $sai_file_2.'.bam';
						
						my $cmd = qq(${BWA_bin} aln -q10 -t42 ${combined_genome_file_f1} ${reads_file} > ${sai_file_1});
						print $cmd, "\n", `$cmd`, "\n";
						$cmd = qq(${BWA_bin} samse ${combined_genome_file_f1} ${sai_file_1} ${reads_file} | ${samtools_bin} view -Sb - > ${bam_file_1});
						print $cmd, "\n", `$cmd`, "\n";
						
						$cmd = qq(${BWA_bin} aln -q10 -t42 ${combined_genome_file_f2} ${reads_file} > ${sai_file_2});
						print $cmd, "\n", `$cmd`, "\n";
						$cmd = qq(${BWA_bin} samse ${combined_genome_file_f2} ${sai_file_2} ${reads_file} | ${samtools_bin} view -Sb - > ${bam_file_2});
						print $cmd, "\n", `$cmd`, "\n";
						
						push(@BAMs_toMerge_1, $bam_file_1);
						push(@BAMs_toMerge_2, $bam_file_2);			
						
						# exit;
						
					}  
				}

				print "Merging BAMs... [1]\n";
				if(scalar(@BAMs_toMerge_1) > 1)
				{
					my $cmd = $samtools_bin.' merge '.$merged_BAM_1.' '.join(' ', @BAMs_toMerge_1);
					print $cmd, "\n", `$cmd`, "\n";
					die unless(-e $merged_BAM_1);
				}
				else
				{
					copy($BAMs_toMerge_1[0], $merged_BAM_1);
				}
				
				print "Merging BAMs... [2]\n";	
				if(scalar(@BAMs_toMerge_2) > 1)
				{		
					my $cmd = $samtools_bin.' merge '.$merged_BAM_2.' '.join(' ', @BAMs_toMerge_2);
					print $cmd, "\n", `$cmd`, "\n";
					die unless(-e $merged_BAM_2);
				}
				else
				{
					copy($BAMs_toMerge_2[0], $merged_BAM_2);		
				}
				
				print "Sorting BAMs... \n";	
				my $cmd = $samtools_bin." sort ${merged_BAM_1} ${merged_BAM_1}.sorted";
				print $cmd, "\n", `$cmd`, "\n";
				$cmd = $samtools_bin." sort ${merged_BAM_2} ${merged_BAM_2}.sorted";
				print $cmd, "\n", `$cmd`, "\n";		

				print "Indexing BAMs... \n";	
				$cmd = $samtools_bin." index ${merged_BAM_1}.sorted.bam";
				print $cmd, "\n", `$cmd`, "\n";
				$cmd = $samtools_bin." index ${merged_BAM_2}.sorted.bam";
				print $cmd, "\n", `$cmd`, "\n";				
			}
		}
		
		my $xMHC_SAM_1 = $remapping_dir.'/mapped_xMHC_1.sam';
		my $xMHC_SAM_2 = $remapping_dir.'/mapped_xMHC_2.sam';
		my $xMHC_SAM_1_readDecision = $remapping_dir.'/mapped_xMHC_1.sam.exclusiveReads';
		my $xMHC_SAM_2_readDecision = $remapping_dir.'/mapped_xMHC_2.sam.exclusiveReads';
		my $xMHC_BAM_1_readDecision = $remapping_dir.'/mapped_xMHC_1.bam.exclusiveReads';
		my $xMHC_BAM_2_readDecision = $remapping_dir.'/mapped_xMHC_2.bam.exclusiveReads';
		
		unless((-e $xMHC_VCF_1_readDecision) and (-e $xMHC_VCF_2_readDecision))
		{
			print "Extracting re-mapped xMHC\n";
			
			my $cmd = $samtools_bin." view -h ${merged_BAM_1}.sorted.bam xMHCpseudoChromHaplotype > ${xMHC_SAM_1}";
			print $cmd, "\n", `$cmd`, "\n";		
			$cmd = $samtools_bin." view -h ${merged_BAM_2}.sorted.bam xMHCpseudoChromHaplotype > ${xMHC_SAM_2}";
			print $cmd, "\n", `$cmd`, "\n";				
			
			print "Deciding where reads come from...\n";
			makeReadOriginDecision($xMHC_SAM_1, $xMHC_SAM_2);
			
			die "Missing file $xMHC_SAM_1_readDecision" unless(-e $xMHC_SAM_1_readDecision);
			die "Missing file $xMHC_SAM_2_readDecision" unless(-e $xMHC_SAM_2_readDecision);
					
			print "Converting files back to BAM...\n";

			$cmd = $samtools_bin." view -S -b -o ${xMHC_BAM_1_readDecision} ${xMHC_SAM_1_readDecision}";
			print $cmd, "\n", `$cmd`, "\n";		

			$cmd = $samtools_bin." view -S -b -o ${xMHC_BAM_2_readDecision} ${xMHC_SAM_2_readDecision}";
			print $cmd, "\n", `$cmd`, "\n";
					
			die unless((-e $xMHC_BAM_1_readDecision) and (-e $xMHC_BAM_2_readDecision));
			
			$cmd = $samtools_bin." sort ${xMHC_BAM_1_readDecision} ${xMHC_BAM_1_readDecision}.sorted";
			print $cmd, "\n", `$cmd`, "\n";		
			$cmd = $samtools_bin." index ${xMHC_BAM_1_readDecision}.sorted.bam";
			print $cmd, "\n", `$cmd`, "\n";		
			
			$cmd = $samtools_bin." sort ${xMHC_BAM_2_readDecision} ${xMHC_BAM_2_readDecision}.sorted";
			print $cmd, "\n", `$cmd`, "\n";		
			$cmd = $samtools_bin." index ${xMHC_BAM_2_readDecision}.sorted.bam";
			print $cmd, "\n", `$cmd`, "\n";		
			
			print "Executing Platypus...\n";

			$cmd = qq(python ${platypus_executable} callVariants --bamFiles=${xMHC_BAM_1_readDecision}.sorted.bam --output=${xMHC_VCF_1_readDecision} --refFile=${combined_genome_file_f1} --logFileName=${xMHC_VCF_1_readDecision}.log --nCPU=1 --mergeClusteredVariants=1 --region=xMHCpseudoChromHaplotype);			
			print $cmd, "\n", `$cmd`, "\n";
			
			$cmd = qq(python ${platypus_executable} callVariants --bamFiles=${xMHC_BAM_2_readDecision}.sorted.bam --output=${xMHC_VCF_2_readDecision} --refFile=${combined_genome_file_f2} --logFileName=${xMHC_VCF_2_readDecision}.log --nCPU=1 --mergeClusteredVariants=1 --region=xMHCpseudoChromHaplotype);			
			print $cmd, "\n", `$cmd`, "\n";
		}
	}
	die unless((-e $xMHC_VCF_1_readDecision) && (-e $xMHC_VCF_2_readDecision));
	my $modifyHaplotypeAccordingToVCF = sub {
	
		# this function takes a VCF and tries to modify a haplotype array accordingly
		
		my $haplotype_fields_aref = shift;
		my $haplotype_index_aref = shift;
		my $haplotype_alleles_aref = shift;		
		my $alleles_details_per_positions = shift;
		my $VCF_filename = shift;
		my $positions_alleles_href = shift;
		my $haplotype = shift;
		my $alignment_index_toFill_aref = shift;
		my $haplotype_id = shift;
		
		for(my $i = 0; $i <= $#{$haplotype_fields_aref}; $i++)
		{
			$alignment_index_toFill_aref->[$i] = undef;
			$haplotype_alleles_aref->[$i] = {};
		}
		
		my %_touched_position;
		for(my $i = 0; $i <= $#{$haplotype_index_aref}; $i++)
		{
			die if($_touched_position{$haplotype_index_aref->[$i]});
			$alignment_index_toFill_aref->[$haplotype_index_aref->[$i]] = $i;
			# die Dumper($alignment_index_toFill_aref->[$haplotype_index_aref->[$i]], $haplotype_index_aref->[$i], $i);
			$_touched_position{$haplotype_index_aref->[$i]}++;
			#if(($i >= 2330646) and ($i <= 2330678))
			#{
			#	my $pgf_pos = $pgf_start + $haplotype_index_aref->[$i];
			#	print join(' ', "Add", "i: $i", "
			#}
		}
			
		# print Dumper("Position 4876672", $haplotype_fields_aref->[4876672], join('', @$haplotype_fields_aref[4876672-15 .. 4876672+15])), "\n";
		
		my $new_haplotype_fields_aref = [@$haplotype_fields_aref];
		
		my $indivID;

		open(VCF, '<', $VCF_filename) or die "Cannot open $VCF_filename";
		my @header_fields;
		while(<VCF>)
		{
			my $line = $_;
			chomp($line);
			next unless ($line);
			next if ($line =~ /^\#\#/);
			if($line =~ /^\#/)
			{
				substr($line, 0, 1) = '';
				@header_fields = split(/\t/, $line);
				$indivID = $header_fields[-1];
				next;
			}
			
			my $lineNum = $.;
			
			my @fields = split(/\t/, $line);
			my %line = (mesh @header_fields, @fields);
		
			die unless (exists $line{POS});
			
			next unless ($line{FILTER} eq 'PASS');
			my $ref_start_haplo = $line{POS} - $remapping_flank_length - 1;
			next unless($ref_start_haplo >= 0);
			next if($ref_start_haplo >= length($haplotype));
			
			my $reference = $line{REF};
			my $alternatives = $line{ALT};
			
			my $ref_extracted_ref = substr($haplotype, $ref_start_haplo, length($reference));
			unless($ref_extracted_ref eq $reference)
			{
				die Dumper("Reference matching problem!", "VCF $VCF_filename line ".$., $ref_start_haplo, "Reference according to VCF: $reference", "Reference from reference haplotype: $ref_extracted_ref", "Context: ".substr($haplotype, $ref_start_haplo-3, 6+length($reference)));
			}
			
			my @alternatives = split(/,/, $alternatives);
			unshift(@alternatives, $reference);
			
			die unless ($line{$indivID});
			
			my $genotype_data_firstPart = $line{FORMAT};
			my $genotype_data_secondPart =  $line{$indivID};
			my @genotype_specifier = split(/\:/, $genotype_data_firstPart);
			die Dumper($line{FORMAT}, \@genotype_specifier) unless($genotype_specifier[0] eq 'GT');
			my @genotype_actualValues = split(/\:/, $genotype_data_secondPart);
			my $variant_genotype = $genotype_actualValues[0];
			my @alleles = split(/\//, $variant_genotype);
			
			my @variantAlleles = grep {$_ ne $reference}  map {(defined $alternatives[$_]) ? $alternatives[$_] : die} @alleles;
			
			if($line{POS} == 3687165)
			{
				my $position_in_alignment = $haplotype_index_aref->[$ref_start_haplo];			
				print Dumper($line{POS}, $ref_start_haplo, $position_in_alignment, $reference, $alternatives, \@alleles, \@variantAlleles, '$positions_alleles_href->{$position_in_alignment}: '.$positions_alleles_href->{$position_in_alignment}, '$new_haplotype_fields_aref->[$position_in_alignment]: '.$new_haplotype_fields_aref->[$position_in_alignment]), "\n";
			}
			
			my $line_id = 'VCF'.$haplotype_id.'_'.$line{POS}.'_refAllele_'.$reference;
			
			foreach my $insertAllele (@variantAlleles)
			{
				# $position_in_VCF specifies the position in the underlying haplotype estimate - 0-based index (?)
				# ... the flanks thing does not work, therefore ignore
				# my $position_in_VCF = $ref_start_haplo;
				die unless(defined $haplotype_index_aref->[$ref_start_haplo]);
				my $position_in_alignment = $haplotype_index_aref->[$ref_start_haplo];
				die unless($position_in_alignment <= $#{$haplotype_fields_aref});
				# now $position_in_alignment tells us which position in the original alignment this position refers to
				
				$positions_alleles_href->{$position_in_alignment} = $insertAllele;
	
				# make sure that the reference we get from the VCF agrees with the haplotype we called variants from
				for(my $i = 0; $i < length($reference); $i++)
				{
					my $c1 = substr($reference, $i, 1);
					die unless(defined $haplotype_index_aref->[$ref_start_haplo + $i]);
					my $pos_of_this_character_in_alignment = $haplotype_index_aref->[$ref_start_haplo + $i];
					my $c2 = $haplotype_fields_aref->[$pos_of_this_character_in_alignment];
					
					die Dumper($position_in_alignment.' of '.$#{$haplotype_fields_aref}.' (i: '.$i.')', 'c1: '.$c1, 'c2: '.$c2, 'ref: '.$reference, 'var[0]: '.$variantAlleles[0]) unless($c1 eq $c2);
				}
				
				# stratify according to the length of the variant allele											
			
				if(length($insertAllele) == length($reference))
				{
					# length identical, assume SNP
					# die "Variant and reference have same length, but > 1 -- this implementation cannot deal with this. Split into multiple 'SNPs'." unless(length($insertAllele) == 1);
					for(my $j = 0; $j < length($insertAllele); $j++)
					{					
						my $this_position_in_alignment = $haplotype_index_aref->[$ref_start_haplo+$j];
						die unless(defined  $this_position_in_alignment);
						die unless (defined $new_haplotype_fields_aref->[$this_position_in_alignment]);									
						$new_haplotype_fields_aref->[$this_position_in_alignment] = substr($insertAllele, $j, 1);	
						$haplotype_alleles_aref->[$this_position_in_alignment]{substr($insertAllele, $j, 1)}++;
						$alleles_details_per_positions->[$this_position_in_alignment]{substr($insertAllele, $j, 1)}{'SNP_'.$line_id.'_ALIGNMENTPOS_'.$this_position_in_alignment}++;
						# push(@{$alignment_index_toFill_aref->[$position_in_alignment+$j]}, $ref_start_haplo+$j);
					}
				}
				elsif(length($insertAllele) < length($reference))
				{
					# deletion
					
					my $pos_of_deletion_start = $haplotype_index_aref->[$ref_start_haplo+length($insertAllele)]; # first character to be deleted
					my $deletion_length = length($reference) - length($insertAllele); # length of deletion
					
					# make sure that the first shared characters agree - this is not VCF spec, but we impose it anyway
					for(my $i = 0; $i < length($insertAllele); $i++)
					{
						my $c1 = substr($insertAllele, $i, 1);
						my $c2 = substr($reference, $i, 1);
						die unless($c1 eq $c2);
					}
					
					# print join(' ', 'DELetion', $line{POS}, $position_in_alignment, $pos_of_deletion_start, length($reference) - length($insertAllele) ), "\n";
					
					# delete deleted characters from haplotype
					for(my $i = 0; $i < $deletion_length; $i++)
					{					
						my $this_position_in_alignment = $haplotype_index_aref->[$ref_start_haplo+length($insertAllele)+$i];
						die unless(defined $this_position_in_alignment);
						die unless (defined $new_haplotype_fields_aref->[$this_position_in_alignment]);	
					
						my $beforeValue = $new_haplotype_fields_aref->[$this_position_in_alignment];
						$new_haplotype_fields_aref->[$this_position_in_alignment] = '_';
						$haplotype_alleles_aref->[$this_position_in_alignment]{'_'}++;		
						$alleles_details_per_positions->[$this_position_in_alignment]{'_'}{'DEL_'.$line_id.'_ALIGNMENTPOS_'.$this_position_in_alignment}++;						

						if($line{POS} == 4634505)
						{
							# print "\t ", $i, " alignmentPos ", $this_position_in_alignment,  " => ", $new_haplotype_fields_aref->[$pos_of_deletion_start + $i], "(replaces $beforeValue)\n";
						}
						# push(@{$alignment_index_toFill_aref->[$pos_of_deletion_start + $i]}, $ref_start_haplo);						
					}
				}
				elsif(length($insertAllele) > length($reference))
				{
					# insertion
					
					# make sure that the first characters between REF and ALT agree - again, not VCF spec, but 
					# we impose it
					for(my $i = 0; $i < length($reference); $i++)
					{
						my $c1 = substr($insertAllele, $i, 1);
						my $c2 = substr($reference, $i, 1);
						warn Dumper("INDEL - first chars of REF and ALT not agreeing", $VCF_filename, $lineNum, $c1, $c2) unless($c1 eq $c2);
						next;
					}
					
					# We want to find out what region of the alignment is covered
					# by this insertion. We will then extract the corresponding sequence and
					# re-align the insertion to the sequence. Gaps are subsituted with their reference
					# characters (from PGF) if possible.
					
					my $insertion_first_position_in_alignment = $position_in_alignment;
					my $insertion_ref_last_position_in_VCF = $ref_start_haplo + length($reference); # the is one after the last position!
					
					my $insertion_last_position_in_alignment;
			
					if(($insertion_ref_last_position_in_VCF - 1) == $#{$haplotype_index_aref})
					{
						# if the last position of the ref allele in the VCF is at the end of the
						# haplotype, the last position in the alignment for the allele is the total last
						# position
						$insertion_last_position_in_alignment = $#{$haplotype_fields_aref};
					}
					else
					{
						# the last position of the ref allele in the VCF is NOT at the end of the haplotype
						# i.e. we have a defined alignment position for the following ref allele, and this minus 1
						# has to be the last position for the reference allele, possibly expanded by subsequent gaps
						die unless(defined $haplotype_index_aref->[$insertion_ref_last_position_in_VCF]);
						$insertion_last_position_in_alignment = $haplotype_index_aref->[$insertion_ref_last_position_in_VCF] - 1;
					}
					die unless($insertion_last_position_in_alignment >= $insertion_first_position_in_alignment);
					
					
					# print "\nInsertion $line{POS}: starts at $insertion_first_position_in_alignment in alignment and goes to $insertion_last_position_in_alignment\n";					
					
					# this is a hack to keep alignment complexity down
					if(($insertion_last_position_in_alignment - $insertion_first_position_in_alignment) > 2*length($insertAllele))
					{
						$insertion_last_position_in_alignment = 2*length($insertAllele) + $insertion_first_position_in_alignment;
						# print "\t decrease complexity -- set insertion end to $insertion_last_position_in_alignment\n";
					}
				
					# we now know beginning and end in the alignment of the bit of the underlying haplotype which served
					# as a basis for discovering the variant we are dealing with. Some positions in this haplotype are
					# potentially gaps. We will these with PGF positions if possible.
					my %got_PGF_positions;
					my $augmented_reference_sequence;
					my @augmented_reference_sequence_alignmentPositions;
					for(my $i = $insertion_first_position_in_alignment; $i <= $insertion_last_position_in_alignment; $i++)
					{
						my $alignment_char = $haplotype_fields_aref->[$i];
						# print join(' ', "\t\tcheckPGFpos", $i), "\n";
						
						if($i == $insertion_first_position_in_alignment)
						{
							# we need to make sure we don't get the first character en double
							die if($alignment_char eq '_'); # we check this again, a bit below!
							die unless(defined $loci_in_graph[$i]);
							my @locusID_fields = split(/_/, $loci_in_graph[$i]);
							die unless($#locusID_fields == 2);
							my $pgf_position = $locusID_fields[2]; # this is not a genomic position but an index into the PGF string
							$got_PGF_positions{$pgf_position} = 1;
						}

						if($alignment_char eq '_')
						{
							die if ($i == $insertion_first_position_in_alignment); # first position may not be gap
							die unless(defined $loci_in_graph[$i]);
							my @locusID_fields = split(/_/, $loci_in_graph[$i]);
							die unless($#locusID_fields == 2);
							my $pgf_position = $locusID_fields[2]; # this is not a genomic position but an index into the PGF string
							
							if($line{POS} eq '4843752')
							{
								#print join(' ', "\t\tcheckPGF", $i, $pgf_position, $loci_in_graph[$i]), "\n";
							}
								
							unless($got_PGF_positions{$pgf_position})
							{
								# the alignment can have gaps etc, we fill the first occurrence of a PGF position
								# with the respective character
								$got_PGF_positions{$pgf_position} = 1;
								$augmented_reference_sequence .= substr($pgf_ref_seq, $pgf_position, 1);
								push(@augmented_reference_sequence_alignmentPositions, $i);
								
								if($line{POS} eq '4843752')
								{
									#print join(' ', "\t\tFromPGF", $pgf_position,  substr($pgf_ref_seq, $pgf_position, 1)), "\n";
								}
					
							}
						}
						else
						{
							$augmented_reference_sequence .= $alignment_char;
							push(@augmented_reference_sequence_alignmentPositions, $i);							
						}
						
						$new_haplotype_fields_aref->[$i] = undef; # set the target region in the haplotype to undef - will be re-filled in a bit
						# push(@{$alignment_index_toFill_aref->[$i]}, $ref_start_haplo);						
						
					}
					die unless(scalar(@augmented_reference_sequence_alignmentPositions) == length($augmented_reference_sequence));
					
					# transpose alignment into new haplotype
					if($line{POS} eq '4843752')
					{
						#print "$line{POS}\n\taugmented ref: $augmented_reference_sequence\n";
					}								
					
					# print "Needleman-Wunsch:\n\taugmented reference: $augmented_reference_sequence\n\tquery: $insertAllele\n";
					
					# align the variant to the specific reference
					my @aligned_sequences = needleman_wunsch($augmented_reference_sequence, $insertAllele, 2, -1, -1);
					
					# print "\taligned ref seq: $aligned_sequences[0]\n\taligned query: $aligned_sequences[1]\n";
					

					# transpose alignment into new haplotype
					if($line{POS} eq '4843752')
					{
						# print "$line{POS} write into new haplotype!\n";
					}					
					
					my $lastAlignmentPositionInString = -1;
					for(my $i = 0; $i < length($aligned_sequences[0]); $i++)
					{
						my $aligned_ref_char = substr($aligned_sequences[0], $i, 1);
						my $aligned_variant_char = substr($aligned_sequences[1], $i, 1);
						die unless((defined $aligned_ref_char) and (defined $aligned_variant_char));
						
						if($aligned_ref_char ne '_')
						{  
							$lastAlignmentPositionInString++;
						}
						die unless($lastAlignmentPositionInString >= 0);
						
						my $alignmentPos = $augmented_reference_sequence_alignmentPositions[$lastAlignmentPositionInString];
						$new_haplotype_fields_aref->[$alignmentPos] .= $aligned_variant_char;

						# print "\t", join(' ', $i, $alignmentPos, $aligned_variant_char), "\n" if($line{POS} eq '4852866');
						# push(@{$alignment_index_toFill_aref->[$alignmentPos]}, $ref_start_haplo);						
						
					}
					
					# if there are any undef fields left, they have to be gaps
					for(my $i = $insertion_first_position_in_alignment; $i <= $insertion_last_position_in_alignment; $i++)
					{
						if(not defined $new_haplotype_fields_aref->[$i])
						{
							$new_haplotype_fields_aref->[$i] = '_';
							$haplotype_alleles_aref->[$i]{'_'}++;						
							$alleles_details_per_positions->[$i]{'_'}{'IN_DEL_'.$line_id.'_ALIGNMENTPOS_'.$i}++;						

							# print "\t", join(' ', '!!!!', $i, 'GAP'), "\n" if($line{POS} eq '4852866');													
							# push(@{$alignment_index_toFill_aref->[$i]}, $ref_start_haplo);													
						}
						else
						{
							$haplotype_alleles_aref->[$i]{$new_haplotype_fields_aref->[$i]}++;						
							$alleles_details_per_positions->[$i]{$new_haplotype_fields_aref->[$i]}{'IN_'.$line_id.'_ALIGNMENTPOS_'.$i}++;												
						}
					}
					
					# print "Final reocnstructed haplotype:\n", join('', @{$new_haplotype_fields_aref}[$insertion_first_position_in_alignment .. $insertion_last_position_in_alignment]), "\n";
					# print "\t in separate fields: ", join(', ', @{$new_haplotype_fields_aref}[$insertion_first_position_in_alignment .. $insertion_last_position_in_alignment]), "\n";
				}
			}
		}
		close(VCF);	
		
		for(my $i = 0; $i <= $#{$haplotype_fields_aref}; $i++)
		{
			my @existing_alleles = keys %{$haplotype_alleles_aref->[$i]};
			my $allele_sum = 0;
			for(@existing_alleles){$allele_sum += $haplotype_alleles_aref->[$i]{$_};}
			warn "Weird allele sum: $allele_sum" unless(($allele_sum == 0) or ($allele_sum == 1) or ($allele_sum == 2));
			if($allele_sum == 0)
			{
				$haplotype_alleles_aref->[$i]{$haplotype_fields_aref->[$i]} += 2;
				$alleles_details_per_positions->[$i]{$haplotype_fields_aref->[$i]}{'REF_2_VCF'.$haplotype_id.'_ALIGNMENTPOS_'.$i}++;						
				
			}
			elsif($allele_sum == 1)
			{
				$haplotype_alleles_aref->[$i]{$haplotype_fields_aref->[$i]} += 1;
				$alleles_details_per_positions->[$i]{$haplotype_fields_aref->[$i]}{'REF_1_VCF'.$haplotype_id.'_ALIGNMENTPOS_'.$i}++;										
			}			
			# $alignment_index_toFill_aref->[$i] = undef;
		}
		
		return $new_haplotype_fields_aref;
	};
	

	#die unless($haplotype_2_alignment_fields[5137237] eq '_');
	#die unless($haplotype_2_alignment_fields[5137238] eq '_');
	
	# modify haplotypes according to the VCFs
	my $aligment_positions_to_haplotypes = [[], []];
	my $positions_alleles_href = {};
	my $alleles_per_positions_1 = [];
	my $alleles_per_positions_2 = [];	
	my $alleles_details_per_positions = [];


	my $ameneded_haplotype_alignment_fields_1 = $modifyHaplotypeAccordingToVCF->(\@haplotype_1_alignment_fields, \@haplotype_1_alignment_index, $alleles_per_positions_1, $alleles_details_per_positions, $xMHC_VCF_1_readDecision, $positions_alleles_href, $haplotype_1, $aligment_positions_to_haplotypes->[0], 1);
	

	for(my $i = 0; $i < $#{$ameneded_haplotype_alignment_fields_1}; $i++)
	{
		$positions_alleles_href->{$i} = $ameneded_haplotype_alignment_fields_1->[$i];
	}

	
	my $ameneded_haplotype_alignment_fields_2 = $modifyHaplotypeAccordingToVCF->(\@haplotype_2_alignment_fields, \@haplotype_2_alignment_index, $alleles_per_positions_2, $alleles_details_per_positions, $xMHC_VCF_2_readDecision, $positions_alleles_href, $haplotype_2, $aligment_positions_to_haplotypes->[1], 2);
		
	#die if($ameneded_haplotype_alignment_fields_2->[5137237] eq '_');
	#die if($ameneded_haplotype_alignment_fields_2->[5137238] eq '_');


	die unless($#{$alleles_per_positions_1} == $#{$alleles_per_positions_2});
	die unless($#{$alleles_per_positions_1} == $#{$ameneded_haplotype_alignment_fields_1});
	die unless($#{$alleles_per_positions_1} == $#{$ameneded_haplotype_alignment_fields_2});

	for(my $i = 0; $i <= $#{$alleles_per_positions_1}; $i++)
	{
		my %combined_alleles;
		foreach my $key (keys %{$alleles_per_positions_1->[$i]}){$combined_alleles{$key} += $alleles_per_positions_1->[$i]{$key}}
		foreach my $key (keys %{$alleles_per_positions_2->[$i]}){$combined_alleles{$key} += $alleles_per_positions_2->[$i]{$key}}
		
		my @sorted_alleles = sort {$combined_alleles{$b} <=> $combined_alleles{$a}} keys %combined_alleles;
		if(scalar(@sorted_alleles) > 1)
		{
			die unless($combined_alleles{$sorted_alleles[0]} >= $combined_alleles{$sorted_alleles[1]});
		}
		die unless($#sorted_alleles >= 0);
		if($#sorted_alleles > 1)
		{
			print Dumper("Warning - more than two alleles at position!", \@sorted_alleles, \%combined_alleles, $alleles_details_per_positions->[$i]), "\n";
		}
		
		if($#sorted_alleles > 0)
		{
			if((exists $alleles_per_positions_1->[$i]{$sorted_alleles[0]}) and ($alleles_per_positions_1->[$i]{$sorted_alleles[0]} > 0))
			{
				# $sorted_alleles[0] can come from haplo 1 - can it also come from haplo 2?
				if((exists $alleles_per_positions_2->[$i]{$sorted_alleles[0]}) and ($alleles_per_positions_2->[$i]{$sorted_alleles[0]} > 0))
				{	
					# $sorted_alleles[0] can come from both haplotypes - what about $sorted_alleles[1]?
					if((exists $alleles_per_positions_1->[$i]{$sorted_alleles[1]}) and ($alleles_per_positions_1->[$i]{$sorted_alleles[1]} > 0))
					{
						# $sorted_alleles[1] can come from haplo 1 as well. Can it come from haplo 2 as well?
						if((exists $alleles_per_positions_2->[$i]{$sorted_alleles[1]}) and ($alleles_per_positions_2->[$i]{$sorted_alleles[1]} > 0))
						{
							# both alleles can come from both haplotypes. selection arbitrary.
							if(($haplotype_1_alignment_fields[$i] eq $sorted_alleles[0]) or ($haplotype_2_alignment_fields[$i] eq $sorted_alleles[1]))
							{
								die unless($alleles_per_positions_1->[$i]{$sorted_alleles[0]} > 0);			
								die unless($alleles_per_positions_2->[$i]{$sorted_alleles[1]} > 0);												
								$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[0];		
								$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[1];										
							}
							else
							{
								die unless($alleles_per_positions_1->[$i]{$sorted_alleles[1]} > 0);			
								die unless($alleles_per_positions_2->[$i]{$sorted_alleles[0]} > 0);												
								$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[1];		
								$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[0];											
							}
						}
						else
						{
							# $sorted_alleles[1] can NOT come from haplo 2, thus it comes from haplo 1!
							die unless($alleles_per_positions_1->[$i]{$sorted_alleles[1]} > 0);			
							die unless($alleles_per_positions_2->[$i]{$sorted_alleles[0]} > 0);												
							$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[1];		
							$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[0];										
						}
					}
					else
					{
						# $sorted_alleles[1] can come from haplo 2 exclusively
						die unless($alleles_per_positions_2->[$i]{$sorted_alleles[1]} > 0);					
						$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[0];		
						$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[1];					
					}
				}
				else
				{
					# $sorted_alleles[0]  can come ONLY from haplo 1.
					warn "An allele we are going to insert into a haplotype is not supported by the VCF output.\nAlleles $sorted_alleles[0], $sorted_alleles[1].\nSupported: ".join(', ', keys %{$alleles_per_positions_1->[$i]}).' // '.join(', ', keys %{$alleles_per_positions_2->[$i]})."\n\n"	unless($alleles_per_positions_2->[$i]{$sorted_alleles[1]} > 0);					
					$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[0];		
					$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[1];					
				}
			}
			else
			{
				# $sorted_alleles[0] can NOT come from haplo 1 thus it goes to haplo 2. Haplo 2 must be able to support the other allele!
			
				die unless($alleles_per_positions_2->[$i]{$sorted_alleles[0]} > 0);
				$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[1];		
				$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[0];			
			}
		}
		else
		{
			$ameneded_haplotype_alignment_fields_1->[$i] = $sorted_alleles[0];				
			$ameneded_haplotype_alignment_fields_2->[$i] = $sorted_alleles[0];		
		}
	}
	

		
	# save haplotypes to file
	my $output_amended_haplotypes = $kMer_count_sample.'.amendedHaplotypes';
	open(AMENDEDHAPLOTYPES, '>', $output_amended_haplotypes) or die "Cannot open $output_amended_haplotypes";
	print AMENDEDHAPLOTYPES join(' ', $sample, 1, @$ameneded_haplotype_alignment_fields_1), "\n";
	print AMENDEDHAPLOTYPES join(' ', $sample, 2, @$ameneded_haplotype_alignment_fields_2), "\n";
	close(AMENDEDHAPLOTYPES);
	
	# save haplotype positions to file
	my $haplo1_to_pgf_aref = alignHaplotypeCharactersToPGF(\@haplotype_1_alignment_fields, \@loci_in_graph);
	die unless(scalar(@$haplo1_to_pgf_aref) == length($haplotype_1));
	my $haplo1_within_alignment_positions = haplotypeCharacterAlignmentPositions(\@haplotype_1_alignment_fields, \@loci_in_graph);
	die unless(scalar(@$haplo1_within_alignment_positions) == length($haplotype_1));	
	my $output_amended_haplotype_1_positions = $kMer_count_sample.'.rawAmendedHaplotype_1.positionsToReference';
	open(AMENDEDHAPLOTYPESPOS, '>', $output_amended_haplotype_1_positions) or die "Cannot open $output_amended_haplotype_1_positions";
	print AMENDEDHAPLOTYPESPOS join(' ', 'pos_in_h_1', 'pos_genomic', 'pos_pgf', 'pos_alignment'), "\n";
	for(my $i = 0; $i <= $#{$haplo1_to_pgf_aref}; $i++)
	{
		my $position_in_pgf = $haplo1_to_pgf_aref->[$i];
		die unless(defined $position_in_pgf);
		my $position_in_alignment = $haplo1_within_alignment_positions->[$i];
		die unless(defined $position_in_alignment);		
		print AMENDEDHAPLOTYPESPOS join(' ', $i, $position_in_pgf, $position_in_pgf,-$pgf_start, $position_in_alignment), "\n";	
	} 
	close(AMENDEDHAPLOTYPESPOS);	
	
	my $haplo2_to_pgf_aref = alignHaplotypeCharactersToPGF(\@haplotype_2_alignment_fields, \@loci_in_graph);
	die unless(scalar(@$haplo2_to_pgf_aref) == length($haplotype_2));
	my $haplo2_within_alignment_positions = haplotypeCharacterAlignmentPositions(\@haplotype_2_alignment_fields, \@loci_in_graph);
	die unless(scalar(@$haplo2_within_alignment_positions) == length($haplotype_2));	
	my $output_amended_haplotype_2_positions = $kMer_count_sample.'.rawAmendedHaplotype_2.positionsToReference';
	open(AMENDEDHAPLOTYPESPOS, '>', $output_amended_haplotype_2_positions) or die "Cannot open $output_amended_haplotype_2_positions";
	print AMENDEDHAPLOTYPESPOS join(' ', 'pos_in_h_1', 'pos_genomic', 'pos_pgf', 'pos_alignment'), "\n";
	for(my $i = 0; $i <= $#{$haplo2_to_pgf_aref}; $i++)
	{
		my $position_in_pgf = $haplo2_to_pgf_aref->[$i];
		die unless(defined $position_in_pgf);
		my $position_in_alignment = $haplo2_within_alignment_positions->[$i];
		die unless(defined $position_in_alignment);		
		print AMENDEDHAPLOTYPESPOS join(' ', $i, $position_in_pgf, $position_in_pgf,-$pgf_start, $position_in_alignment), "\n";	
	} 
	close(AMENDEDHAPLOTYPESPOS);	
	
	
	# generate VCF
	my $output_amendedVCF = $kMer_count_sample.'.amendedHaplotypes.VCF';
	haplotypesToVCF($ameneded_haplotype_alignment_fields_1, $ameneded_haplotype_alignment_fields_2, $sample, \@loci_in_graph, $pgf_ref_seq, $output_amendedVCF, $aligment_positions_to_haplotypes);

	exit;
}
elsif($collect eq '2')
{
	die "Legacy code. Use --collect 2viterbi instead.";
	
	my $expected_output_filename = $kMer_count_sample.'.viterbiHaplotypes';
	
	# read loci from graph specification
	
	my @loci_in_graph;
	my $graph_segments_file = $graph.'/segments.txt';
	die "File not there: $graph_segments_file" unless (-e $graph_segments_file);
	open(SEGMENTS, '<', $graph_segments_file) or die;
	while(<SEGMENTS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
		next unless($line);
		
		my $segment_file = $graph.'/'.$line;

		open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
		my $firstLine = <SEGMENT>;
		chomp($firstLine);
		$firstLine =~ s/\n//g;
		$firstLine =~ s/\r//g;			
		my @line_fields = split(/ /, $firstLine);
		shift(@line_fields); # kick out individual ID
		push(@loci_in_graph, @line_fields);
		close(SEGMENT);
	}
	close(SEGMENTS);
	
	my %locus_position;
	my $positions_file = $graph.'/positions.txt';
	die "File not there: $positions_file" unless (-e $positions_file);
	open(POSITIONS, '<', $positions_file) or die;	
	while(<POSITIONS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;		
		my @fields = split(/ /, $line);
		my $field_id = $fields[0];
		my $field_position = $fields[1];
		$locus_position{$field_id} = $field_position;
	}
	close(POSITIONS);
	
	for(my $lI = 1; $lI <= $#loci_in_graph; $lI++)
	{
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI]", $lI);
		}
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI-1]", $lI-1);
		}		
		die unless($locus_position{$loci_in_graph[$lI]} > $locus_position{$loci_in_graph[$lI-1]});
	}
	
	# read xMHC reference
	
	open(REF, "<", $xMHC_reference) or die "Cannot open $xMHC_reference";
	my $pgf_ref_seq;
	while(<REF>)
	{
		my $line = $_; 
		chomp($line);
		$line =~ s/[\x0A\x0D]+//g;	
		next if ($line =~ /^\>/);
		$pgf_ref_seq .= $line;
	}
	close(REF);
	my $pgf_stop = $pgf_start+length($pgf_ref_seq);
	
	my %identifier_to_position;
	my %identifier_to_refalleles;
	if(-e $vcfPos)
	{
		open(VCFSNPS, '<', $vcfPos) or die "Cannot open $vcfPos";
		my $header_line = <VCFSNPS>;
		chomp($header_line);
		$header_line =~ s/\n//g;
		$header_line =~ s/\r//g;	
		
		my @header_fields = split(/\t/, $header_line);
		while(<VCFSNPS>)
		{
			my $line = $_;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;	
		
			my @line_fields = split(/\t/, $line);
			my %line = (mesh @header_fields, @line_fields);
		
			die unless(defined $line{Position});
			
			die unless($line{RefAllele});
			
			my $reference_charPos_PGF = $line{Position} - $pgf_start;
			
			next if($reference_charPos_PGF < 0);
			next if($reference_charPos_PGF > (length($pgf_ref_seq) - 1));
			
			$identifier_to_position{$line{ID}} = $line{Position};
			$identifier_to_refalleles{$line{ID}} = $line{RefAllele};
						
			my $reference_char = substr($pgf_ref_seq, $reference_charPos_PGF, length($line{RefAllele}));
			
			die "$reference_char ne $line{RefAllele}, $reference_charPos_PGF, line $." unless($reference_char eq $line{RefAllele});
		}			
		close(VCFSNPS);
	}
	
	if((-e $expected_output_filename))
	{
		print "Analyze $expected_output_filename\n";
		
		my %inference;
		
		my %inference_PP;
		my %inference_seen_pairs;
		
		my %inference_firstAllele;
		my %inference_secondAllele; 
		my %firstAllele_called;
		my %secondAllele_called;		
		
		my %seen_PGF_positions;
		my %pgf_genotypes_labels;
		
		open(INFERENCE, '<', $expected_output_filename) or die "Cannot open $expected_output_filename";				
		while(<INFERENCE>)
		{
			my $line = $_;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;
			
			my @fields = split(/ /, $line);
			
			my $indiv_spec = shift(@fields);
			shift(@fields); # chromosome
			
			die "Error indivID not parseable: $indiv_spec" unless($indiv_spec =~ /^(.+?)-x-(.+?)$/);
			my $individualID = $1;
			my $iteration = $2;
			
			die Dumper($#fields, $#loci_in_graph) unless($#fields == ($#loci_in_graph+2));
			@fields = @fields[0 .. ($#fields-2)];
			
			my %pgf_genotypes;
			
			for(my $pI = 0; $pI <= $#fields; $pI++)
			{
				my $locusID = $loci_in_graph[$pI];
				my @locusID_fields = split(/_/, $locusID);
				die unless($#locusID_fields == 2);

				my $pgf_position = $pgf_start + $locusID_fields[2];
				
				if($locusID_fields[0] =~ /^rs/)
				{
					$pgf_genotypes_labels{$pgf_position}{$locusID_fields[0]}++;
				}
				
				my $char_from_inference = $fields[$pI];
				
				$pgf_genotypes{$pgf_position} .= $char_from_inference;
				
				$seen_PGF_positions{$pgf_position}++;
			}
			
			foreach my $pgfPos (keys %pgf_genotypes)
			{
				push(@{$inference{$individualID}{$pgfPos}[$iteration-1]}, $pgf_genotypes{$pgfPos});
			}
		}
		close(INFERENCE);
			
		foreach my $individualID (keys %inference)
		{
			foreach my $pgfPos (keys %seen_PGF_positions)
			{
				foreach my $inference_combination (@{$inference{$individualID}{$pgfPos}})
				{
					warn "Empty field in inference combinations!" unless ($inference_combination);
					my @alleles = sort @{$inference_combination};
					
					my $alleles = join('_', @alleles);
					$inference_seen_pairs{$pgfPos}{$alleles} = 1;

		
					$inference_PP{$individualID}{$pgfPos}{$alleles}++;
					$inference_PP{$individualID}{$pgfPos}{'ALL'}++;			
		
					
					my %tmp_alleles = map {$_ => 1} @alleles;
					

					foreach my $allele (keys %tmp_alleles)
					{
						$inference_firstAllele{$individualID}{$pgfPos}{$allele}++;
					}
					$inference_firstAllele{$individualID}{$pgfPos}{'ALL'}++;		

				}
			}
		}
		
		foreach my $individualID (keys %inference)
		{
			foreach my $pgfPos (keys %seen_PGF_positions)
			{					
				my $call = maxAllele($inference_firstAllele{$individualID}{$pgfPos});
				$firstAllele_called{$individualID}{$pgfPos} = $call;
				$firstAllele_called{$individualID}{$pgfPos}[2] = $firstAllele_called{$individualID}{$pgfPos}[1];
			}
		}
			
		foreach my $individualID (keys %inference)
		{
			foreach my $pgfPos (keys %seen_PGF_positions)
			{		
				foreach my $inference_combination (@{$inference{$individualID}{$pgfPos}})
				{
					warn "Empty field in inference combinations!" unless ($inference_combination);
					my @alleles = sort @{$inference_combination};
					my $firstAllele = $firstAllele_called{$individualID}{$pgfPos}[0];
					if($firstAllele ~~ @alleles)
					{
						my $gotFirst = 0;
						my $otherAllele;
						for(@alleles)
						{
							if(! $gotFirst)
							{
								if($_ eq $firstAllele)
								{
									$gotFirst = 1;
									next;
								}
							}
							$otherAllele = $_;
						}
						die unless ($gotFirst);
						
						$inference_secondAllele{$individualID}{$pgfPos}{$otherAllele}++;
						$inference_secondAllele{$individualID}{$pgfPos}{'ALL'}++;	
					}
				}
			}
		}
			
		foreach my $individualID (keys %inference)
		{
			foreach my $pgfPos (keys %seen_PGF_positions)
			{		
				my $call = maxAllele($inference_secondAllele{$individualID}{$pgfPos});
				$secondAllele_called{$individualID}{$pgfPos} = $call;
				$secondAllele_called{$individualID}{$pgfPos}[2] = $secondAllele_called{$individualID}{$pgfPos}[1] * $firstAllele_called{$individualID}{$pgfPos}[2];
			}	
		}		
		
		
		my %max_pair_value;
		my %max_pair_alleles;
		my $output_file = $expected_output_filename.".PP";
		open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
		my @allele_headers;
		my %pairs_in_order;
			foreach my $pgfPos (keys %seen_PGF_positions)
		{	
			my @pairs_in_order = sort keys %{$inference_seen_pairs{$pgfPos}};
			$pairs_in_order{$pgfPos} = \@pairs_in_order;
			push(@allele_headers, @pairs_in_order);
		}
		print OUTPUT join(' ', 'IndividualID', @allele_headers), "\n";
		foreach my $individualID (sort keys %inference_PP)
		{
			my @allele_fields;
			foreach my $pgfPos (keys %seen_PGF_positions)
			{
				push(@allele_fields,
					map {sprintf("%.5f", $inference_PP{$individualID}{$pgfPos}{$_}/$inference_PP{$individualID}{$pgfPos}{'ALL'} )} @{$pairs_in_order{$pgfPos}}
				);
			}
			print OUTPUT join(' ', $individualID, @allele_fields), "\n";
			print OUTPUT join(' ', $individualID, ), "\n";
		}
		close(OUTPUT);
		
		$output_file = $expected_output_filename.".bestguess";
		open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
		
		@allele_headers = ();
		foreach my $pgfPos (keys %seen_PGF_positions)
		{	
			my @headers = ($pgfPos, qw/Q Q2 P/);			
			push(@allele_headers, @headers);
		}
		
		my %bestguess_exist_alleles;
		
		print OUTPUT join(' ', 'IndividualID', 'Chromosome', @allele_headers), "\n";
		foreach my $individualID (keys %inference_firstAllele)
		{	
			my $likelihood = 'NA';
			my @f_1;
			my @f_2;
			foreach my $pgfPos (keys %seen_PGF_positions)
			{	
				$bestguess_exist_alleles{$pgfPos}{$firstAllele_called{$individualID}{$pgfPos}[0]}++;
				$bestguess_exist_alleles{$pgfPos}{$secondAllele_called{$individualID}{$pgfPos}[0]}++;
				
				push(@f_1,
					$firstAllele_called{$individualID}{$pgfPos}[0], $firstAllele_called{$individualID}{$pgfPos}[1], $firstAllele_called{$individualID}{$pgfPos}[2], $likelihood
				
				);
				push(@f_2,
					$secondAllele_called{$individualID}{$pgfPos}[0], $secondAllele_called{$individualID}{$pgfPos}[1], $secondAllele_called{$individualID}{$pgfPos}[2], $likelihood
				);				
			}
			print OUTPUT join(' ', $individualID, 1, @f_1), "\n";
			print OUTPUT join(' ', $individualID, 2, @f_2), "\n";
		}
		close(OUTPUT);
		
		if($vcfPos)
		{
			my $output_file = $expected_output_filename.".VCF";
			open(VCF, '>', $output_file) or die "Cannot open $output_file";
my $VCF_header = qq(##fileformat=VCFv4.0
##fileDate=27/04/12
##phasing=none, though some calls involve phasing clustered variants
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=COV,Number=2,Type=Integer,Description="Number of reads on ref and alt alleles">
##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">
##FORMAT=<ID=SITE_CONF,Number=1,Type=Float,Description="Probabilitic site classification confidence. Difference in log likelihood of most likely and next most likely model (models are variant, repeat and error)">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=SNP_FROM_COMPLEX,Description="SNP called from a cluster of phased SNPs or complex SNP/indel , split out for easier comparison with other SNP call sets">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INV_INDEL,Description="Inversion+indel - this script overcalls these, so worth checking">
##ALT=<ID=DEL_INV,Description="Deletion + Inversion">
##ALT=<ID=INS_INV,Description="Insertion + Inversion">
##ALT=<ID=PH_SNPS,Description="Phased SNPs">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##FILTER=<ID=MAPQ,Description="5prime flank maps to reference with mapping quality below 40">);
			my @header_fields = qw/CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/;
			my @individualIDs = sort (keys %inference_firstAllele);
			print VCF $VCF_header, "\n";
			print VCF '#'.join("\t", @header_fields, @individualIDs), "\n";
			
			my @output_variants = (keys %seen_PGF_positions);
			foreach my $pgfPos (sort {$a <=> $b} @output_variants)
			{
				my $variantID = '.';
				my $rsID;
				
				if(scalar( keys %{$pgf_genotypes_labels{$pgfPos}}) == 1)
				{
					$rsID = ((keys %{$pgf_genotypes_labels{$pgfPos}})[0]);
					die Dumper('$pgfPos != $identifier_to_position{$rsID}', $rsID, $pgfPos, $identifier_to_position{$rsID}) unless($pgfPos == $identifier_to_position{$rsID});
					$variantID = $rsID;
				}
				
				my @output_fields;
				push(@output_fields, 6);
				push(@output_fields, $pgfPos);
				push(@output_fields, $variantID);
				
				my $reference_charPos_PGF = $pgfPos - $pgf_start;
				die "???" if($reference_charPos_PGF < 0);
				die "???" if($reference_charPos_PGF > (length($pgf_ref_seq) - 1));
				my $reference_char = substr($pgf_ref_seq, $reference_charPos_PGF, 1);

				if($rsID)
				{
					die Dumper($., $pgfPos, $rsID, $identifier_to_refalleles{$rsID}, $reference_char) unless($identifier_to_refalleles{$rsID} eq $reference_char);
				}
				
				my $refAllele = $reference_char;
				push(@output_fields, $refAllele);
				
				die unless($bestguess_exist_alleles{$pgfPos});
				delete $bestguess_exist_alleles{$pgfPos}{$refAllele};
				my @nonref_alleles = sort keys %{$bestguess_exist_alleles{$pgfPos}};
				
				my %allele_2_num = ($refAllele => 0);
				my $allele_counter = 1;
				for(@nonref_alleles)
				{
					$allele_2_num{$nonref_alleles[$allele_counter-1]} = $allele_counter;
					$allele_counter++
				}
				if(scalar(@nonref_alleles) == 0)
				{
					push(@output_fields, '.');
				}
				else
				{
					push(@output_fields, join(',', @nonref_alleles));
				}
				push(@output_fields, '.');
				# this is a hack!
				die if($#individualIDs > 0);
				if(not ($secondAllele_called{$individualIDs[0]}{$pgfPos}[2] > 0.7))
				{
					push(@output_fields, 'FAIL');
				}
				else
				{
					push(@output_fields, 'PASS');
				}
				push(@output_fields, 'SVTYPE=SNP;SVLEN=0');
				push(@output_fields, 'GT');
				
				foreach my $indivID (@individualIDs)
				{
					(exists $allele_2_num{$firstAllele_called{$indivID}{$pgfPos}[0]}) or die;
					(exists $allele_2_num{$secondAllele_called{$indivID}{$pgfPos}[0]}) or die;					
					my $gt = join('/', $allele_2_num{$firstAllele_called{$indivID}{$pgfPos}[0]}, $allele_2_num{$secondAllele_called{$indivID}{$pgfPos}[0]});
					push(@output_fields, $gt);
				}
				
				print VCF join("\t", @output_fields), "\n";
			}
			
			close(VCF);
			
			print "\n\nVCF file $output_file produced.\n";
		}
		
		print "\n\nOutput collected in\n$output_file\n\n";		
		
		print "\n\n";
		
		# print `cat $output_file`;
	}
	else
	{
		print "Cannot find inference file $expected_output_filename\n";
	}
}
elsif($collect eq '2viterbi')
{
	my $expected_output_filename = $kMer_count_sample.'.viterbiHaplotypes';
	
	# read loci from graph specification
	
	my @loci_in_graph;
	my $graph_segments_file = $graph.'/segments.txt';
	die "File not there: $graph_segments_file" unless (-e $graph_segments_file);
	open(SEGMENTS, '<', $graph_segments_file) or die;
	while(<SEGMENTS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
		next unless($line);
		
		my $segment_file = $graph.'/'.$line;

		open(SEGMENT, '<', $segment_file) or die "Cannot open $segment_file";
		my $firstLine = <SEGMENT>;
		chomp($firstLine);
		$firstLine =~ s/\n//g;
		$firstLine =~ s/\r//g;			
		my @line_fields = split(/ /, $firstLine);
		shift(@line_fields); # kick out individual ID
		push(@loci_in_graph, @line_fields);
		close(SEGMENT);
	}
	close(SEGMENTS);
	
	my %locus_position;
	my $positions_file = $graph.'/positions.txt';
	die "File not there: $positions_file" unless (-e $positions_file);
	open(POSITIONS, '<', $positions_file) or die;	
	while(<POSITIONS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;		
		my @fields = split(/ /, $line);
		my $field_id = $fields[0];
		my $field_position = $fields[1];
		$locus_position{$field_id} = $field_position;
	}
	close(POSITIONS);
	
	for(my $lI = 1; $lI <= $#loci_in_graph; $lI++)
	{
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI]", $lI);
		}
		if(not defined $locus_position{$loci_in_graph[$lI]})
		{
			die Dumper("No position information for locus $loci_in_graph[$lI-1]", $lI-1);
		}		
		die unless($locus_position{$loci_in_graph[$lI]} > $locus_position{$loci_in_graph[$lI-1]});
	}
	
	# read xMHC reference
	
	open(REF, "<", $xMHC_reference) or die "Cannot open $xMHC_reference";
	my $pgf_ref_seq;
	while(<REF>)
	{
		my $line = $_; 
		chomp($line);
		$line =~ s/[\x0A\x0D]+//g;	
		next if ($line =~ /^\>/);
		$pgf_ref_seq .= $line;
	}
	close(REF);
	my $pgf_stop = $pgf_start+length($pgf_ref_seq);
	

	if((-e $expected_output_filename))
	{
		print "Analyze $expected_output_filename\n";
		
		my %inference;
		
		my %inference_PP;
		my %inference_seen_pairs;
		
		my %inference_firstAllele;
		my %inference_secondAllele; 
		my %firstAllele_called;
		my %secondAllele_called;		
		
		my %seen_PGF_positions;
		my %pgf_genotypes_labels;
		
		open(HAPLOTYPES, '<', $expected_output_filename) or die "Cannot open $expected_output_filename";
		my $haplotype_1, my $haplotype_2;
		my @haplotype_1_alignment_fields;
		my @haplotype_2_alignment_fields;
		while(<HAPLOTYPES>)
		{
			my $line = $_;
			chomp($line);
			last if ($. > 2);
			my @fields = split(/ /, $line);
			my $haplotype = join('', @fields[2 .. ($#fields-2)]);
			if($. == 1)
			{
				$haplotype_1 = $haplotype;
				@haplotype_1_alignment_fields = @fields[2 .. ($#fields-2)];
			}
			elsif($. == 2)
			{
				$haplotype_2 = $haplotype;	
				@haplotype_2_alignment_fields = @fields[2 .. ($#fields-2)];
			}
		}
		close(HAPLOTYPES);
		
	
		my $output_file = $expected_output_filename.".viterbiVCF";		
		haplotypesToVCF(\@haplotype_1_alignment_fields, \@haplotype_2_alignment_fields, $sample, \@loci_in_graph, $pgf_ref_seq, $output_file);
	}
	else
	{
		print "Cannot find inference file $expected_output_filename\n";
	}
}
elsif($collect eq '1')
{
	my $expected_output_filename = $kMer_count_sample.'.inference';
	
	open(REF, "<", $xMHC_reference) or die "Cannot open $xMHC_reference";
	my $pgf_ref_seq;
	while(<REF>)
	{
		my $line = $_; 
		chomp($line);
		$line =~ s/[\x0A\x0D]+//g;	
		next if ($line =~ /^\>/);
		$pgf_ref_seq .= $line;
	}
	close(REF);
	my $pgf_stop = $pgf_start+length($pgf_ref_seq);
	
	my %identifier_to_position;
	my %identifier_to_refalleles;
	if(-e $vcfPos)
	{
		open(VCFSNPS, '<', $vcfPos) or die "Cannot open $vcfPos";
		my $header_line = <VCFSNPS>;
		chomp($header_line);
		my @header_fields = split(/\t/, $header_line);
		while(<VCFSNPS>)
		{
			my $line = $_;
			chomp($line);
			my @line_fields = split(/\t/, $line);
			my %line = (mesh @header_fields, @line_fields);
		
			die unless(defined $line{Position});
			
			die unless($line{RefAllele});
			
			my $reference_charPos_PGF = $line{Position} - $pgf_start;
			
			next if($reference_charPos_PGF < 0);
			next if($reference_charPos_PGF > (length($pgf_ref_seq) - 1));
			$identifier_to_position{$line{ID}} = $line{Position};
			$identifier_to_refalleles{$line{ID}} = $line{RefAllele};
						
			my $reference_char = substr($pgf_ref_seq, $reference_charPos_PGF, length($line{RefAllele}));
			
			die "$reference_char ne $line{RefAllele}, $reference_charPos_PGF, line $." unless($reference_char eq $line{RefAllele});
		}			
		close(VCFSNPS);
	}
	
	if((-e $expected_output_filename) && (! $redo))
	{
		print "Analyze $expected_output_filename\n";
		
		my %inference;
		
		my %inference_PP;
		my %inference_seen_pairs;
		
		my %inference_firstAllele;
		my %inference_secondAllele; 
		my %firstAllele_called;
		my %secondAllele_called;		
		
		open(INFERENCE, '<', $expected_output_filename) or die "Cannot open $expected_output_filename";
		my $firstLine_inference = <INFERENCE>;
		chomp($firstLine_inference);
		$firstLine_inference =~ s/\n//g;
		$firstLine_inference =~ s/\r//g;		
		my @header_fields = split(/ /, $firstLine_inference);
		die unless($header_fields[0] eq 'IndividualID');
		my @inference_fields =  @header_fields[2 .. $#header_fields];
				
		while(<INFERENCE>)
		{
			my $line = $_;
			chomp($line);
			$line =~ s/\n//g;
			$line =~ s/\r//g;
			
			my @fields = split(/ /, $line);
			my %linehash = (mesh @header_fields, @fields);
			
			my $individualID_mix = $linehash{'IndividualID'};
						
			die "Error indivID not parseable: $individualID_mix" unless($individualID_mix =~ /^(.+?)-x-(.+?)$/);
			my $individualID = $1;
			my $iteration = $2;
			
			foreach my $fN (@inference_fields)
			{
				die "Found question marks in output for the desired fieldname $fN -- not good!\nValue: $linehash{$fN}" if ($linehash{$fN} =~ /\?/);
				
				if($fN =~ /HLA/)
				{
					if(length($linehash{$fN}) > 5)
					{
						$linehash{$fN} = substr($linehash{$fN}, 0, 5);
					}	
				}
				
				push(@{$inference{$individualID}{$fN}[$iteration-1]}, $linehash{$fN});
			}			
		}
		close(INFERENCE);
			
		foreach my $individualID (keys %inference)
		{
			foreach my $fN (@inference_fields)
			{
				foreach my $inference_combination (@{$inference{$individualID}{$fN}})
				{
					warn "Empty field in inference combinations!" unless ($inference_combination);
					my @alleles = sort @{$inference_combination};
					
					my $alleles = join('_', @alleles);
					$inference_seen_pairs{$fN}{$alleles} = 1;

		
					$inference_PP{$individualID}{$fN}{$alleles}++;
					$inference_PP{$individualID}{$fN}{'ALL'}++;			
		
					
					my %tmp_alleles = map {$_ => 1} @alleles;
					

					foreach my $allele (keys %tmp_alleles)
					{
						$inference_firstAllele{$individualID}{$fN}{$allele}++;
					}
					$inference_firstAllele{$individualID}{$fN}{'ALL'}++;		

				}
			}
		}
		
		foreach my $individualID (keys %inference)
		{
			foreach my $fN (@inference_fields)
			{					
				my $call = maxAllele($inference_firstAllele{$individualID}{$fN});
				$firstAllele_called{$individualID}{$fN} = $call;
				$firstAllele_called{$individualID}{$fN}[2] = $firstAllele_called{$individualID}{$fN}[1];
			}
		}
			
		foreach my $individualID (keys %inference)
		{
			foreach my $fN (@inference_fields)
			{		
				foreach my $inference_combination (@{$inference{$individualID}{$fN}})
				{
					warn "Empty field in inference combinations!" unless ($inference_combination);
					my @alleles = sort @{$inference_combination};
					my $firstAllele = $firstAllele_called{$individualID}{$fN}[0];
					if($firstAllele ~~ @alleles)
					{
						my $gotFirst = 0;
						my $otherAllele;
						for(@alleles)
						{
							if(! $gotFirst)
							{
								if($_ eq $firstAllele)
								{
									$gotFirst = 1;
									next;
								}
							}
							$otherAllele = $_;
						}
						die unless ($gotFirst);
						
						$inference_secondAllele{$individualID}{$fN}{$otherAllele}++;
						$inference_secondAllele{$individualID}{$fN}{'ALL'}++;	

					}
				}
			}
		}
			
		foreach my $individualID (keys %inference)
		{
			foreach my $fN (@inference_fields)
			{		
				my $call = maxAllele($inference_secondAllele{$individualID}{$fN});
				$secondAllele_called{$individualID}{$fN} = $call;
				$secondAllele_called{$individualID}{$fN}[2] = $secondAllele_called{$individualID}{$fN}[1] * $firstAllele_called{$individualID}{$fN}[2];
			}	
		}		
		
		
		my %max_pair_value;
		my %max_pair_alleles;
		my $output_file = $expected_output_filename.".PP";
		open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
		my @allele_headers;
		my %pairs_in_order;
		foreach my $fN (@inference_fields)
		{	
			my @pairs_in_order = sort keys %{$inference_seen_pairs{$fN}};
			$pairs_in_order{$fN} = \@pairs_in_order;
			push(@allele_headers, @pairs_in_order);
		}
		print OUTPUT join(' ', 'IndividualID', @allele_headers), "\n";
		foreach my $individualID (sort keys %inference_PP)
		{
			my @allele_fields;
			foreach my $fN (@inference_fields)
			{
				push(@allele_fields,
					map {sprintf("%.5f", $inference_PP{$individualID}{$fN}{$_}/$inference_PP{$individualID}{$fN}{'ALL'} )} @{$pairs_in_order{$fN}}
				);
			}
			print OUTPUT join(' ', $individualID, @allele_fields), "\n";
			print OUTPUT join(' ', $individualID, ), "\n";
		}
		close(OUTPUT);
		
		$output_file = $expected_output_filename.".bestguess";
		open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
		
		@allele_headers = ();
		foreach my $fN (@inference_fields)
		{	
			my @headers = ($fN, qw/Q Q2 P/);			
			push(@allele_headers, @headers);
		}
		
		my %bestguess_exist_alleles;
		
		print OUTPUT join(' ', 'IndividualID', 'Chromosome', @allele_headers), "\n";
		foreach my $individualID (keys %inference_firstAllele)
		{	
			my $likelihood = 'NA';
			my @f_1;
			my @f_2;
			foreach my $fN (@inference_fields)
			{	
				$bestguess_exist_alleles{$fN}{$firstAllele_called{$individualID}{$fN}[0]}++;
				$bestguess_exist_alleles{$fN}{$secondAllele_called{$individualID}{$fN}[0]}++;
				
				push(@f_1,
					$firstAllele_called{$individualID}{$fN}[0], $firstAllele_called{$individualID}{$fN}[1], $firstAllele_called{$individualID}{$fN}[2], $likelihood
				
				);
				push(@f_2,
					$secondAllele_called{$individualID}{$fN}[0], $secondAllele_called{$individualID}{$fN}[1], $secondAllele_called{$individualID}{$fN}[2], $likelihood
				);				
			}
			print OUTPUT join(' ', $individualID, 1, @f_1), "\n";
			print OUTPUT join(' ', $individualID, 2, @f_2), "\n";
		}
		close(OUTPUT);
		
		if($vcfPos)
		{
			my $output_file = $expected_output_filename.".VCF";
			open(VCF, '>', $output_file) or die "Cannot open $output_file";
my $VCF_header = qq(##fileformat=VCFv4.0
##fileDate=27/04/12
##phasing=none, though some calls involve phasing clustered variants
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=COV,Number=2,Type=Integer,Description="Number of reads on ref and alt alleles">
##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">
##FORMAT=<ID=SITE_CONF,Number=1,Type=Float,Description="Probabilitic site classification confidence. Difference in log likelihood of most likely and next most likely model (models are variant, repeat and error)">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=SNP_FROM_COMPLEX,Description="SNP called from a cluster of phased SNPs or complex SNP/indel , split out for easier comparison with other SNP call sets">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INV_INDEL,Description="Inversion+indel - this script overcalls these, so worth checking">
##ALT=<ID=DEL_INV,Description="Deletion + Inversion">
##ALT=<ID=INS_INV,Description="Insertion + Inversion">
##ALT=<ID=PH_SNPS,Description="Phased SNPs">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##FILTER=<ID=MAPQ,Description="5prime flank maps to reference with mapping quality below 40">);
			my @header_fields = qw/CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/;
			my @individualIDs = sort (keys %inference_firstAllele);
			print VCF $VCF_header, "\n";
			print VCF '#'.join("\t", @header_fields, @individualIDs), "\n";
			
			my @output_variants = grep {exists $identifier_to_position{$_}} @inference_fields;
			foreach my $variantID (@output_variants)
			{
				my @output_fields;
				push(@output_fields, 6);
				push(@output_fields, $identifier_to_position{$variantID});
				push(@output_fields, $variantID);
				
				my $refAllele = $identifier_to_refalleles{$variantID};
				push(@output_fields, $refAllele);
				
				die unless($bestguess_exist_alleles{$variantID});
				delete $bestguess_exist_alleles{$variantID}{$refAllele};
				my @nonref_alleles = sort keys %{$bestguess_exist_alleles{$variantID}};
				
				my %allele_2_num = ($refAllele => 0);
				my $allele_counter = 1;
				for(@nonref_alleles)
				{
					$allele_2_num{$nonref_alleles[$allele_counter-1]} = $allele_counter;
					$allele_counter++
				}
				if(scalar(@nonref_alleles) == 0)
				{
					push(@output_fields, '.');
				}
				else
				{
					push(@output_fields, join(',', @nonref_alleles));
				}
				push(@output_fields, '.');
				push(@output_fields, 'PASS');
				push(@output_fields, 'SVTYPE=SNP;SVLEN=0');
				push(@output_fields, 'GT');
				
				foreach my $indivID (@individualIDs)
				{
					(exists $allele_2_num{$firstAllele_called{$indivID}{$variantID}[0]}) or die;
					(exists $allele_2_num{$secondAllele_called{$indivID}{$variantID}[0]}) or die;					
					my $gt = join('/', $allele_2_num{$firstAllele_called{$indivID}{$variantID}[0]}, $allele_2_num{$secondAllele_called{$indivID}{$variantID}[0]});
					push(@output_fields, $gt);
				}
				
				print VCF join("\t", @output_fields), "\n";
			}
			
			close(VCF);
			
			print "\n\nVCF file $output_file produced.\n";
		}
		
		print "\n\nOutput collected in\n$output_file\n\n";		
		
		print "\n\n";
		
		# print `cat $output_file`;
	}
	else
	{
		print "Cannot find inference file $expected_output_filename\n";
	}
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

sub maxAllele
{
	my $indiv_href = shift;
	my $max_allele;
	my $max_allele_Q = -1; 
	foreach my $allele (keys %$indiv_href)
	{
		next if ($allele eq 'ALL');
		if($indiv_href->{$allele} > $max_allele_Q)
		{
			$max_allele = $allele;
			$max_allele_Q = $indiv_href->{$allele};
		}
	}
	
	die "Field $indiv_href->{'ALL'} missing\n\n".Dumper($indiv_href) unless($indiv_href->{'ALL'});
	
	return [$max_allele, $max_allele_Q/$indiv_href->{'ALL'}];
}


sub find_individual_readsPath
{
	# my $individualID = shift;
	# my @possibilities = map {$base_path_reads.'/drive_'.$_.'/'.$individualID} qw/f l v w x/;
	# my @found = grep {(-e $_) && (-d $_)} @possibilities;
	# unless ($#found == 0)
	# {
		# die "Problem finding $individualID. Details:\n".Dumper(\@found);
	# }
	# return $found[0];
}

sub find_individual_BAM
{
	my $individualID = shift;
	my $BAM = $base_path_BAMs.'/'.$individualID.'.bam';
	unless(-e $BAM)
	{
		# die "Cannot find BAM for individual $individualID. Expected path:\n$BAM";
	}
	return $BAM;
}


sub makeReadOriginDecision
{
	my $input_1 = shift;
	my $input_2 = shift;
	
	my $output_1 = $input_1.'.exclusiveReads';
	my $output_2 = $input_2.'.exclusiveReads';
	
	my $mapped_flag = 2**2;
	die unless($mapped_flag == 0x4);
	
	my @header_lines_1;
	my @header_lines_2;
	
	my %reads;
	
	open(INPUT1, '<', $input_1) or die "Cannot open $input_1";
	open(INPUT2, '<', $input_2) or die "Cannot open $input_2";
	open(OUTPUT1, '>', $output_1) or die "Cannot open $output_1";
	open(OUTPUT2, '>', $output_2) or die "Cannot open $output_2";
	while(<INPUT1>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 1) eq '@')
		{
			push(@header_lines_1, $line);
		}
		else
		{
			my @fields = split(/\t/, $line);
			
			my $name = $fields[0];
			my $flags = $fields[1];
			my $quality = $fields[4];
			my $sequence = $fields[9];
			
			next if($flags & $mapped_flag);
			next if($quality == 255);
			
			my $edit_distance;
			unless($line =~ /NM:i:(\d+)/)
			{
				die "Weird line $. in $input_1 :\n$line";
			}
			$edit_distance = $1;
								
			$reads{$sequence}{$name}[0]{q} = $quality;
			$reads{$sequence}{$name}[0]{edit_distance} = $edit_distance;						
			$reads{$sequence}{$name}[0]{line} = $line;			
		}
	}
	
 
	while(<INPUT2>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 1) eq '@')
		{
			push(@header_lines_2, $line);
		}
		else
		{
			my @fields = split(/\t/, $line);
			
			my $name = $fields[0];
			my $flags = $fields[1];
			my $quality = $fields[4];
			my $sequence = $fields[9];
			
			next if($flags & $mapped_flag);
			next if($quality == 255);
			
			my $edit_distance;
			unless($line =~ /NM:i:(\d+)/)
			{
				die "Weird line $. in $input_2 :\n$line";
			}
			$edit_distance = $1;  
							
			$reads{$sequence}{$name}[1]{q} = $quality;
			$reads{$sequence}{$name}[1]{edit_distance} = $edit_distance;			
			$reads{$sequence}{$name}[1]{line} = $line;			
		}
	}	
	
	print OUTPUT1 join("\n", @header_lines_1), "\n";
	print OUTPUT2 join("\n", @header_lines_2), "\n";
	
	foreach my $sequence (keys %reads)
	{
		my @names = keys %{$reads{$sequence}};
		if($#names > 0)
		{
			print "Sequence $sequence, have ", scalar(@names), " names\n";
		}
		foreach my $name (@names)
		{
			my $q1 = $reads{$sequence}{$name}[0]{q};
			my $q2 = $reads{$sequence}{$name}[1]{q};
			my $e1 = $reads{$sequence}{$name}[0]{edit_distance};
			my $e2 = $reads{$sequence}{$name}[1]{edit_distance};
			
			my $num_of_qualities = (( (defined $q1) ? 1 : 0) + ( (defined $q2) ? 1 : 0 ) );
			my $num_of_editDistances = (( (defined $e1) ? 1 : 0) + ( (defined $e2) ? 1 : 0 ) );
			
			# print "\tName $name with $num_of_qualities qualities.\n" if ($#names > 0);
			
			$e1 = 1000 unless(defined $e1);
			$e2 = 1000 unless(defined $e2);
			die if(($e1 == -1) and ($e2 == -1));
			
			$q1 = -1 unless(defined $q1);
			$q2 = -1 unless(defined $q2);
			die if(($q1 == -1) and ($q2 == -1));
			
			# print $q1, " vs ", $q2, "\n";
			
			# if($e1 < $e2)
			# {
					# print OUTPUT1 $reads{$sequence}{$name}[0]{line}, "\n";			
			# }
			# elsif($e1 > $e2)
			# {
					# print OUTPUT2 $reads{$sequence}{$name}[1]{line}, "\n";			
			# }
			# else
			# {
				# die Dumper('$e1 != $e2', $e1, $e2) unless($e1 == $e2);
				
				if($q1 > $q2)
				{
					print OUTPUT1 $reads{$sequence}{$name}[0]{line}, "\n";
				}
				elsif($q1 == $q2)
				{
					if(rand() < 0.5)
					{
						die Dumper($q1, $q2)  unless(defined $reads{$sequence}{$name}[0]{line});
						print OUTPUT1 $reads{$sequence}{$name}[0]{line}, "\n";				
					}
					else
					{
						die Dumper($q1, $q2) unless(defined $reads{$sequence}{$name}[1]{line});
						print OUTPUT2 $reads{$sequence}{$name}[1]{line}, "\n";				
					}
				}
				else
				{
					print OUTPUT2 $reads{$sequence}{$name}[1]{line}, "\n";
				}
			# }
		}
	}
	
	
	close(INPUT1);
	close(INPUT2);
	close(OUTPUT1);
	close(OUTPUT2);
}
 

sub needleman_wunsch
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $score_match = shift;
	my $score_mismatch = shift;
	my $score_gap = shift;
	
	my $length_x = length($seq1);
	my $length_y = length($seq2);
	
	my $score_firstGap = $score_gap - 1000;
	
	my $t = [];
	for(my $x = 0; $x <= $length_x; $x++)
	{
		$t->[$x][0] = $x * $score_firstGap;
	}
	
	for(my $y = 0; $y <= $length_y; $y++)
	{
		$t->[0]->[$y] = $y * $score_firstGap;
	}
	
	for(my $y = 1; $y <= $length_y; $y++)
	{
		for(my $x = 1; $x <= $length_x; $x++)
		{
			die unless(defined substr($seq1, $x-1, 1));
			die unless(defined substr($seq2, $y-1, 1));
			my @alt = ($t->[$x-1][$y] + $score_gap, $t->[$x][$y-1] + $score_gap, $t->[$x-1][$y-1] + (((substr($seq1, $x-1, 1) eq substr($seq2, $y-1, 1)) or (substr($seq1, $x-1, 1) eq '*') or (substr($seq1, $x-1, 1) eq 'N')) ? $score_match : $score_mismatch));
				
			$t->[$x][$y] = more_max(@alt);
		}
	}
	
	my $backtrack_x = $length_x, my $backtrack_y = $length_y;
	
	my $seq1_aligned, my $seq2_aligned;
	while(($backtrack_x != 0) or ($backtrack_y != 0))
	{
		die unless($backtrack_x >= 0);
		die unless($backtrack_y >= 0);
		
		
		die if( (($backtrack_x == 0) and ($backtrack_y == 1)) or (($backtrack_x == 1) and ($backtrack_y == 0)));
		
		my $v = $t->[$backtrack_x][$backtrack_y];
		my @alts = ($t->[$backtrack_x-1][$backtrack_y] + $score_gap, $t->[$backtrack_x][$backtrack_y-1] + $score_gap, $t->[$backtrack_x-1][$backtrack_y-1] + (((substr($seq1, $backtrack_x-1, 1) eq substr($seq2, $backtrack_y-1, 1)) or (substr($seq1, $backtrack_x-1, 1) eq '*') or (substr($seq1, $backtrack_x-1, 1) eq 'N')) ? $score_match : $score_mismatch));
		my $whereFrom = more_max_from(@alts);
		#print Dumper(@alts, $whereFrom), "\n\n";
		
		if($backtrack_x == 0)
		{
			$whereFrom = 1;
		}
		if($backtrack_y == 0)
		{
			$whereFrom = 0;
		}			
		
		if($whereFrom == 0)
		{
			$seq2_aligned .= '_';
			$seq1_aligned .= substr($seq1, $backtrack_x-1, 1);
			$backtrack_x--;
		}
		elsif($whereFrom == 1)
		{
			$seq1_aligned .= '_';
			$seq2_aligned .= substr($seq2, $backtrack_y-1, 1);
			$backtrack_y--;
		}
		elsif($whereFrom == 2)
		{
			$seq1_aligned .= substr($seq1, $backtrack_x-1, 1);
			$seq2_aligned .= substr($seq2, $backtrack_y-1, 1);
			$backtrack_x--;
			$backtrack_y--;
		}
		else
		{
			die;
		}
	}
	
	$seq1_aligned = reverse($seq1_aligned);
	$seq2_aligned = reverse($seq2_aligned);
	
	return ($seq1_aligned, $seq2_aligned);
}

sub more_max
{
	return (reverse sort {$a <=> $b}(@_))[0];
}

sub more_max_from
{
	my @alts = @_;
	my $m = more_max(@alts);
	my %from;
	for(my $i = 0; $i <= $#alts; $i++)
	{
		if($alts[$i] eq $m)
		{
			return $i;
		}
	}
	die;
}	

sub alignHaplotypeCharactersToPGF
{
	my $haplotype_aref = shift;
	my $loci_in_graph_aref = shift;
	
	die unless(scalar(@$haplotype_aref) == scalar(@$loci_in_graph_aref));
	
	my @forReturn;
	for(my $pI = 0; $pI <= $#{$loci_in_graph_aref}; $pI++)
	{
		my $locusID = $loci_in_graph_aref->[$pI];
		my @locusID_fields = split(/_/, $locusID);
		die unless($#locusID_fields == 2);
		my $haplo_char = $haplotype_aref->[$pI];
		if($haplo_char ne '_')
		{
			my $pgf_position = $pgf_start + $locusID_fields[2];
			push(@forReturn, $pgf_position);
		}	
	}
	
	return \@forReturn;
}

sub haplotypeCharacterAlignmentPositions
{
	my $haplotype_aref = shift;
	my $loci_in_graph_aref = shift;
	
	die unless(scalar(@$haplotype_aref) == scalar(@$loci_in_graph_aref));
	
	my @forReturn;
	for(my $pI = 0; $pI <= $#{$loci_in_graph_aref}; $pI++)
	{
		my $locusID = $loci_in_graph_aref->[$pI];
		my @locusID_fields = split(/_/, $locusID);
		die unless($#locusID_fields == 2);
		my $haplo_char = $haplotype_aref->[$pI];
		if($haplo_char ne '_')
		{
			push(@forReturn, $pI);
		}	
	}
	
	return \@forReturn;
}


sub haplotypesToVCF
{
	my $haplo1_fields_aref = shift;
	my $haplo2_fields_aref = shift;
	my $individualID = shift;
	my $loci_in_graph_aref = shift;
	my $pgf_ref_seq = shift;
	my $output_file = shift;
	my $individual_haplotypes_positions = shift;
	
	my @loci_in_graph = @$loci_in_graph_aref;
	 
	my %identifier_to_position;   
	my %identifier_to_refalleles;
	
	die "Please suppy parameter --vcfPos" unless ((defined $vcfPos));
	die "Specified --vcfPos ($vcfPos) not existing" unless(-e $vcfPos);

	open(VCFSNPS, '<', $vcfPos) or die "Cannot open $vcfPos";
	my $header_line = <VCFSNPS>;
	chomp($header_line);
	$header_line =~ s/\n//g;
	$header_line =~ s/\r//g;	
	
	my @header_fields = split(/\t/, $header_line);
	while(<VCFSNPS>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\n//g;
		$line =~ s/\r//g;	
	
		my @line_fields = split(/\t/, $line);
		my %line = (mesh @header_fields, @line_fields);
	
		die unless(defined $line{Position});
		
		die unless($line{RefAllele});
		
		my $reference_charPos_PGF = $line{Position} - $pgf_start;
		
		next if($reference_charPos_PGF < 0);
		next if($reference_charPos_PGF > (length($pgf_ref_seq) - 1));
		
		$identifier_to_position{$line{ID}} = $line{Position};
		$identifier_to_refalleles{$line{ID}} = $line{RefAllele};
					
		my $reference_char = substr($pgf_ref_seq, $reference_charPos_PGF, length($line{RefAllele}));
		
		die "$reference_char ne $line{RefAllele}, $reference_charPos_PGF, line $." unless($reference_char eq $line{RefAllele});
	}			
	close(VCFSNPS);
	
	my %inference;
	
	my %inference_PP;
	my %inference_seen_pairs;
	
	my %inference_firstAllele;
	my %inference_secondAllele; 
	my %firstAllele_called;
	my %secondAllele_called;		
	my %bestguess_exist_alleles;
	
	my %seen_PGF_positions;
	my %pgf_genotypes_labels;
	
	my %taken_alignment_positions;
	my %taken_haplotype_positions;
	
	my $iteration = 1; # this is just one iteration
	foreach my $fRef ($haplo1_fields_aref, $haplo2_fields_aref)
	{
		my @fields = @$fRef;
		
		die Dumper($#fields, $#loci_in_graph) unless($#fields == $#loci_in_graph);
		
		my %pgf_genotypes;
		
		my $nextHaplo = 1;
		my %prepared_this_PGF_index;
		for(my $pI = 0; $pI <= $#fields; $pI++)
		{
			my $locusID = $loci_in_graph[$pI];
			my @locusID_fields = split(/_/, $locusID);
			die unless($#locusID_fields == 2);

			my $pgf_position = $pgf_start + $locusID_fields[2];
			unless($prepared_this_PGF_index{$pgf_position})
			{
				unless($taken_alignment_positions{$pgf_position})
				{	
					$taken_alignment_positions{$pgf_position} = [];
				}
				push(@{$taken_alignment_positions{$pgf_position}}, []);
				
				unless($taken_haplotype_positions{$pgf_position})
				{	
					$taken_haplotype_positions{$pgf_position} = [];
				}
				push(@{$taken_haplotype_positions{$pgf_position}}, []);		
				
				$prepared_this_PGF_index{$pgf_position} = 1;
			}
		}
		
		for(my $pI = 0; $pI <= $#fields; $pI++)
		{
			my $locusID = $loci_in_graph[$pI];
			my @locusID_fields = split(/_/, $locusID);
			die unless($#locusID_fields == 2);

			my $pgf_position = $pgf_start + $locusID_fields[2];
			
			if($locusID_fields[0] =~ /^rs/)
			{
				$pgf_genotypes_labels{$pgf_position}{$locusID_fields[0]}++;
			}
			
			my $char_from_inference = $fields[$pI];
			
			$pgf_genotypes{$pgf_position} .= $char_from_inference;
			
			$seen_PGF_positions{$pgf_position}++;
			
			die unless($#{$taken_alignment_positions{$pgf_position}} > -1);
			push(@{$taken_alignment_positions{$pgf_position}[$#{$taken_alignment_positions{$pgf_position}}]}, $pI);
			
			if($individual_haplotypes_positions)
			{
				my $haplotype_position;
				if($fRef == $haplo1_fields_aref)
				{
					$haplotype_position = $individual_haplotypes_positions->[0][$pI];
				}
				else
				{
					$haplotype_position = $individual_haplotypes_positions->[1][$pI];				
				}

				push(@{$taken_haplotype_positions{$pgf_position}[$#{$taken_haplotype_positions{$pgf_position}}]}, $haplotype_position) if(defined $haplotype_position);
			}
		}
		
		foreach my $pgfPos (keys %pgf_genotypes)
		{
			push(@{$inference{$individualID}{$pgfPos}[$iteration-1]}, $pgf_genotypes{$pgfPos});
		}
	}
	close(INFERENCE);
	
	# print the alignment positions that this VCF refers to
	
	my $output_alignment_positions = $output_file.'.positionsInRefToAlignment';
	
	open(POSITIONS, '>', $output_alignment_positions) or die "Cannot open $output_alignment_positions";
	print POSITIONS join(' ', 'reference_position', 'alignment_positions_h1', 'alignment_positions_h2'), "\n";
	my @k = sort {$a <=> $b} keys %taken_alignment_positions;
	foreach my $k (@k)
	{
		print POSITIONS join(' ', $k, join(';', @{$taken_alignment_positions{$k}[0]}), join(';', @{$taken_alignment_positions{$k}[1]}), ), "\n";
	}
	close(POSITIONS);
	
	if($individual_haplotypes_positions)
	{
		$output_alignment_positions = $output_file.'.positionsInRefToHaplotypes';
		
		open(POSITIONS, '>', $output_alignment_positions) or die "Cannot open $output_alignment_positions";
		print POSITIONS join(' ', 'reference_position', 'haplotype_positions_h1', 'haplotype_positions_h2'), "\n";
		my @k = sort {$a <=> $b} keys %taken_haplotype_positions;
		foreach my $k (@k)
		{
			print POSITIONS join(' ', $k, join(';', @{$taken_haplotype_positions{$k}[0]}), join(';', @{$taken_haplotype_positions{$k}[1]}), ), "\n";
		}
		close(POSITIONS);	
	}

		
	foreach my $individualID (keys %inference)
	{
		foreach my $pgfPos (keys %seen_PGF_positions)
		{
			foreach my $inference_combination (@{$inference{$individualID}{$pgfPos}})
			{
				warn "Empty field in inference combinations!" unless ($inference_combination);
				my @alleles = sort @{$inference_combination};
				
				my $alleles = join('_', @alleles);
				$inference_seen_pairs{$pgfPos}{$alleles} = 1;

				$inference_PP{$individualID}{$pgfPos}{$alleles}++;
				$inference_PP{$individualID}{$pgfPos}{'ALL'}++;			
				
				my %tmp_alleles = map {$_ => 1} @alleles;
			
				foreach my $allele (keys %tmp_alleles)
				{
					$inference_firstAllele{$individualID}{$pgfPos}{$allele}++;
				}
				$inference_firstAllele{$individualID}{$pgfPos}{'ALL'}++;		

			}
		}
	}
		
	foreach my $individualID (keys %inference)
	{
		foreach my $pgfPos (keys %seen_PGF_positions)
		{					
			my $call = maxAllele($inference_firstAllele{$individualID}{$pgfPos});
			$firstAllele_called{$individualID}{$pgfPos} = $call;
			$firstAllele_called{$individualID}{$pgfPos}[2] = $firstAllele_called{$individualID}{$pgfPos}[1];
		}
	}
			
	foreach my $individualID (keys %inference)
	{
		foreach my $pgfPos (keys %seen_PGF_positions)
		{		
			foreach my $inference_combination (@{$inference{$individualID}{$pgfPos}})
			{
				warn "Empty field in inference combinations!" unless ($inference_combination);
				my @alleles = sort @{$inference_combination};
				my $firstAllele = $firstAllele_called{$individualID}{$pgfPos}[0];
				if($firstAllele ~~ @alleles)
				{
					my $gotFirst = 0;
					my $otherAllele;
					for(@alleles)
					{
						if(! $gotFirst)
						{
							if($_ eq $firstAllele)
							{
								$gotFirst = 1;
								next;
							}
						}
						$otherAllele = $_;
					}
					die unless ($gotFirst);
					
					$inference_secondAllele{$individualID}{$pgfPos}{$otherAllele}++;
					$inference_secondAllele{$individualID}{$pgfPos}{'ALL'}++;	
				}
			}
		}
	}
			
	foreach my $individualID (keys %inference)
	{
		foreach my $pgfPos (keys %seen_PGF_positions)
		{		
			my $call = maxAllele($inference_secondAllele{$individualID}{$pgfPos});
			$secondAllele_called{$individualID}{$pgfPos} = $call;
			$secondAllele_called{$individualID}{$pgfPos}[2] = $secondAllele_called{$individualID}{$pgfPos}[1] * $firstAllele_called{$individualID}{$pgfPos}[2];
		}	
	}		
		
	foreach my $individualID (keys %inference_firstAllele)
	{	
		foreach my $pgfPos (keys %seen_PGF_positions)
		{	
			$bestguess_exist_alleles{$pgfPos}{$firstAllele_called{$individualID}{$pgfPos}[0]}++;
			$bestguess_exist_alleles{$pgfPos}{$secondAllele_called{$individualID}{$pgfPos}[0]}++;
		}
	}
		
		
	open(VCF, '>', $output_file) or die "Cannot open $output_file";
	my $VCF_header = qq(##fileformat=VCFv4.0
##fileDate=27/04/12
##phasing=none, though some calls involve phasing clustered variants
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=COV,Number=2,Type=Integer,Description="Number of reads on ref and alt alleles">
##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">
##FORMAT=<ID=SITE_CONF,Number=1,Type=Float,Description="Probabilitic site classification confidence. Difference in log likelihood of most likely and next most likely model (models are variant, repeat and error)">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=SNP_FROM_COMPLEX,Description="SNP called from a cluster of phased SNPs or complex SNP/indel , split out for easier comparison with other SNP call sets">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INV_INDEL,Description="Inversion+indel - this script overcalls these, so worth checking">
##ALT=<ID=DEL_INV,Description="Deletion + Inversion">
##ALT=<ID=INS_INV,Description="Insertion + Inversion">
##ALT=<ID=PH_SNPS,Description="Phased SNPs">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##FILTER=<ID=MAPQ,Description="5prime flank maps to reference with mapping quality below 40">);
	my @header_fields_VCF = qw/CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/;
	
	my @individualIDs = sort (keys %inference_firstAllele);
	print VCF $VCF_header, "\n";
	print VCF '#'.join("\t", @header_fields_VCF, @individualIDs), "\n";
	
	my $last_pgfPos_with_IGV_H1;
	my $last_pgfPos_with_IGV_H2;
	
	my @output_variants = (keys %seen_PGF_positions);
	foreach my $pgfPos (sort {$a <=> $b} @output_variants)
	{
		my $variantID = '.';
		my $rsID;
		
		if(scalar( keys %{$pgf_genotypes_labels{$pgfPos}}) == 1)
		{
			$rsID = ((keys %{$pgf_genotypes_labels{$pgfPos}})[0]);
			die Dumper('$pgfPos != $identifier_to_position{$rsID}', $rsID, $pgfPos, $identifier_to_position{$rsID}) unless($pgfPos == $identifier_to_position{$rsID});
			$variantID = $rsID;
		}
		
		my @output_fields;
		push(@output_fields, 6);
		push(@output_fields, $pgfPos);
		push(@output_fields, $variantID);
		
		my $reference_charPos_PGF = $pgfPos - $pgf_start;
		die "???" if($reference_charPos_PGF < 0);
		die "???" if($reference_charPos_PGF > (length($pgf_ref_seq) - 1));
		my $reference_char = substr($pgf_ref_seq, $reference_charPos_PGF, 1);

		if($rsID)
		{
			die Dumper($., $pgfPos, $rsID, $identifier_to_refalleles{$rsID}, $reference_char) unless($identifier_to_refalleles{$rsID} eq $reference_char);
		}
		
		my $refAllele = $reference_char;
		push(@output_fields, $refAllele);
		
		die unless($bestguess_exist_alleles{$pgfPos});
		delete $bestguess_exist_alleles{$pgfPos}{$refAllele};
		my @nonref_alleles = sort keys %{$bestguess_exist_alleles{$pgfPos}};
		
		my %allele_2_num = ($refAllele => 0);
		my $allele_counter = 1;
		for(@nonref_alleles)
		{
			$allele_2_num{$nonref_alleles[$allele_counter-1]} = $allele_counter;
			$allele_counter++
		}
		if(scalar(@nonref_alleles) == 0)
		{
			push(@output_fields, '.');
		}
		else
		{
			push(@output_fields, join(',', @nonref_alleles));
		}
		push(@output_fields, '.');
		
		# this is a hack!
		die if($#individualIDs > 0);
		if(not ($secondAllele_called{$individualIDs[0]}{$pgfPos}[2] > 0.7))
		{
			push(@output_fields, 'FAIL');
		}
		else
		{
			push(@output_fields, 'PASS');
		}
		
		my $real_pgfPos = $pgfPos - $pgf_start;
		
		my @positions_in_string_haplotypes;
		my @positions_in_IGV_haplotypes;
		my @preDeletion_IGV_positions;
		
		if($individual_haplotypes_positions)
		{
			if($pgfPos == 31009300)
			{
				# die Dumper($taken_haplotype_positions{$pgfPos}[0]);
			}
			@positions_in_string_haplotypes = (join('-', @{$taken_haplotype_positions{$pgfPos}[0]}), join('-', @{$taken_haplotype_positions{$pgfPos}[1]}));
			@positions_in_IGV_haplotypes = (join('-', map {$_ + $remapping_flank_length + 1} @{$taken_haplotype_positions{$pgfPos}[0]}), join('-', map {$_ + $remapping_flank_length + 1} @{$taken_haplotype_positions{$pgfPos}[1]}));
			
			if(scalar(@{$taken_haplotype_positions{$pgfPos}[0]}) == 0)
			{
				die unless($last_pgfPos_with_IGV_H1);
				$preDeletion_IGV_positions[0] = join('-', map {$_ + $remapping_flank_length + 1} @{$taken_haplotype_positions{$last_pgfPos_with_IGV_H1}[0]});				
			}
			else
			{
				$last_pgfPos_with_IGV_H1 = $pgfPos;
				$preDeletion_IGV_positions[0] = '';
			}
			
			if(scalar(@{$taken_haplotype_positions{$pgfPos}[1]}) == 0)
			{
				die unless($last_pgfPos_with_IGV_H2);
				$preDeletion_IGV_positions[1] = join('-', map {$_ + $remapping_flank_length + 1} @{$taken_haplotype_positions{$last_pgfPos_with_IGV_H2}[1]});				
			}
			else
			{
				$last_pgfPos_with_IGV_H2 = $pgfPos;
				$preDeletion_IGV_positions[1] = '';
			}			
		}
		else
		{		
			@positions_in_string_haplotypes = ('', '');
			@positions_in_IGV_haplotypes = ('', '');
			@preDeletion_IGV_positions = ('', '');
		}
		
		my @positions_in_alignments = (join('-', @{$taken_alignment_positions{$pgfPos}[0]}), join('-', @{$taken_alignment_positions{$pgfPos}[1]}));
		
		# warnings seem to be normal if reading Viterbi haplotypes
		push(@output_fields, 'SVTYPE=SNP;SVLEN=0;'.qq(POS_IN_PGF=${real_pgfPos};H1_ALIGNMENT=${positions_in_alignments[0]};H2_ALIGNMENT=${positions_in_alignments[1]};H1_STRING=${positions_in_string_haplotypes[0]};H2_STRING=${positions_in_string_haplotypes[1]};H1_IGV=${positions_in_IGV_haplotypes[0]};H2_IGV=${positions_in_IGV_haplotypes[1]};H1_PREDELETION_IGV=${preDeletion_IGV_positions[0]};H2_PREDELETION_IGV=${preDeletion_IGV_positions[1]};));
		push(@output_fields, 'GT');
		
		foreach my $indivID (@individualIDs)
		{
			(exists $allele_2_num{$firstAllele_called{$indivID}{$pgfPos}[0]}) or die;
			(exists $allele_2_num{$secondAllele_called{$indivID}{$pgfPos}[0]}) or die;					
			my $gt = join('/', $allele_2_num{$firstAllele_called{$indivID}{$pgfPos}[0]}, $allele_2_num{$secondAllele_called{$indivID}{$pgfPos}[0]});
			push(@output_fields, $gt);
		}
		
		print VCF join("\t", @output_fields), "\n";
	}
	
	close(VCF);
}

