

package MHC::Utils::HLAtypeinference; 

use strict;
use warnings;

use Exporter;
use autodie;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw( 
              averageFractionAlignmentOK 
              compatibleAlleles_individual
              compatibleStringAlleles
              compatibleStringAlleles_individual 
              compatibleStringAlleles_noFlip
              exon_folder_files 
              inferReadLength
              list_comparison
              load_coverages_from_pileup
              load_pileup
              min_avg_max              
              print_which_exons
              read_exon_sequences  
              twoClusterAlleles
              twoValidationAlleles_2_proper_names
            );



sub compatibleStringAlleles
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my ($OK1, $NOTOK1, $MISSING1) = compatibleStringAlleles_noFlip($alleles_validation, $alleles_inference);
	
	my $alleles_validation_flipped = [reverse(@$alleles_validation)];
	my ($OK2, $NOTOK2, $MISSING2) = compatibleStringAlleles_noFlip($alleles_validation_flipped, $alleles_inference);
	
	# die unless($MISSING1 == $MISSING2);
	
	if(($NOTOK1 < $NOTOK2))
	{
		return ($OK1, $NOTOK1, $MISSING1);
	}
	elsif($NOTOK1 == $NOTOK2)
	{
		return (($OK1 > $OK2) ? ($OK1, $NOTOK1, $MISSING1) : ($OK2, $NOTOK2, $MISSING2));		
	}
	else
	{
		return ($OK2, $NOTOK2, $MISSING2);
	}
	
	
}


sub compatibleStringAlleles_noFlip
{
	my $alleles_validation = shift;
	my $alleles_inference = shift;

	die unless($#{$alleles_inference} == 1);
	die unless($#{$alleles_validation} == 1);

	my $alleles_OK = 0;
	my $alleles_NOTOK = 0;
	my $alleles_MISSING = 0;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{	
		my ($OK, $NOTOK, $MISSING) = (compatibleStringAlleles_individual(undef, $alleles_validation->[$aI], $alleles_inference->[$aI]));
		$alleles_OK += $OK;
		$alleles_NOTOK += $NOTOK;
		$alleles_MISSING += $MISSING;
	}
		
	return ($alleles_OK, $alleles_NOTOK, $alleles_MISSING);
}



sub compatibleStringAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;	
	my $allele_inference = shift;

	if(($allele_validation =~ /\?/) || ($allele_inference =~ /\?/))
	{
		return (0, 0, 1);
	}
	else
	{
		if($allele_validation eq $allele_inference)
		{
			return (1, 0, 0);
		}
		else
		{
			return (0, 1, 0);
		}
	}
}

sub compatibleAlleles_individual
{
	my $locus = shift;
	my $allele_validation = shift;
	my $allele_inference = shift;
	

	# die Dumper($allele_validation, $allele_inference);
	
	my @components_allele_validation = split(/;/, $allele_validation);
	if(scalar(@components_allele_validation) > 1)
	{
		my $found_compatible = 0;
		foreach my $validation_allele (@components_allele_validation)
		{
			$found_compatible += compatibleAlleles_individual($locus, $validation_allele, $allele_inference);
			last if($found_compatible);
		}
		if($found_compatible)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}	
	
	my $allele_validation_original = $allele_validation;
		
	if($locus eq 'C')
	{
		# print Dumper([$locus, $allele_validation, $allele_inference, scalar(@components_allele_validation)]);
	}
	
	die "Weird allele $locus $allele_validation" unless(length($allele_validation) >= 4);
	die unless($allele_inference =~ /:/);

	
	#die unless($allele_validation =~ /^\d\d/);

	my $components_validation;
	if($allele_validation !~ /\:/)
	{
		die if($allele_validation =~ /\//);
		$allele_validation = substr($allele_validation, 0, 2).':'.substr($allele_validation, 2);
		$components_validation = 2;
		$allele_validation =~ s/g//i;
		
	}
	else
	{
		my @_components = split(/\:/, $allele_validation);
		$components_validation = scalar(@_components);
		die "Weird allele code $allele_validation" if($allele_validation =~ /g/i);
		die if($allele_validation =~ /D/i);
			
	}
	die unless(defined $components_validation);
	die unless($components_validation >= 2);

	
	my $true_allele = $allele_validation;

		
	my $inferred_alleles = $allele_inference;
	my @inferred_alleles = split(/;/, $inferred_alleles);
	
	# if inferred allele has fewer components than validation allele: add 0s, mismatch if validation allele is not 0 at this position
	# => conservative
	# if inferred allele has more components than validation allele, truncate inferred allele
	# e.g. inferred: 02:03:01, validation allele: 03:01, truncated inferred allele: 02:03
	# if inferred and validation allele have different 4-digit groups, we always generate a mismatch
	# if inferred and validation allele have the same 4-digit group, we generate a match
	# e.g. inferred 02:03:01, validation: 02:03
	# is this a problem? No. Consider two scenarios:
	# - validation allele typed to complete 4-digit (amino acid sequence, all exons) resolution: this is what we want
	# - validation allele typed to SBT G-groups: this case is not a problem.
	# 		Case 1: specified as G, i.e. 02:03G
	#			This is a bad example, because 02:03G is not a valid G allele - needs to be six digits, e.g.
	#			02:03:01G. This is fine as long as our output is also at G group level (i.e. no
	#			higher resolution than exon 2 / exons 2,3), as 02:03:01 will be one of our output alleles
	#			if and only if we infer 02:03:01G as a result.
	# 			Gs are confusing and we will deactivate them now.
	#			(Or enable proper translation)
	#  		Case 2: specified as explicit allele list, separated with slashes - we now assume that 02:03 is 
	#   			one such explicit member of a list.
	#			This case does not exist. If 02:03:01 exists as a named allele, other alleles like 02:03:02 must also exist,
	#			and 02:03 would not be an individual member of any G group (because all members with identical
	#			4-digit group would be differentiated by their 6-digit groups).

	
	@inferred_alleles = map {
		# die "Can't parse allele $_" unless($_ =~ /^\s*\w+\*(\d\d\d?)\:(\d\d\d?N?).*$/);
		# my $allele = $1.':'.$2;
		# $allele
		
		die "Can't parse allele $_" unless($_ =~ /^\s*(\w+)\*([\d\:N]+Q?)L?S?$/);
		
		my $allele = $2;
		
		my @components_allele = split(/:/, $allele);
		my @components_allele_rightLength;
		
		for(my $i = 0; $i < $components_validation; $i++)
		{
			if($i <= $#components_allele)
			{
				push(@components_allele_rightLength, $components_allele[$i]);
			}
			else
			{
				push(@components_allele_rightLength, '00');
			}
		}
		$allele = join(':', @components_allele_rightLength);
		
	} @inferred_alleles;
	
	my %inferred = map {$_ => 1} @inferred_alleles;
	
	if($inferred{$true_allele})
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 1), "\n" if ($locus eq 'DQB1');	
		
		return 1;  
	}		
	else
	{
		# print Dumper($allele_validation_original, $allele_inference, \@inferred_alleles, 0), "\n"  if ($locus eq 'DQB1');	
	
		return 0;
	}
} 

sub print_which_exons
{
	my $locus = shift;
	if((length($locus) == 1) or ($locus =~ /HLA\w/))
	{
		return (2, 3);
	}
	else
	{
		return (2);
	}
}	



sub load_pileup
{
	my $r_href = shift;
	my $file = shift;
	my $indivID = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line, -1);
		die unless(($#f == 3) or ($#f == 2));
		if($#f == 3)
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = $f[3];		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		}
		else
		{
			my $exon = $f[0];
			my $exonPos = $f[1];
			my $pileUp = '';		
			$r_href->{$indivID}{'HLA'.$locus}[$exon][$exonPos] = $pileUp;
		}
	}	
	close(F);
}




sub load_coverages_from_pileup
{
	my $file = shift;
	
	die unless($file =~ /R\d+_pileup_(\w+).txt$/);
	my $locus = $1;
	
	my %forReturn;
	
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die "Cannot parse pileup coverage line" unless($line =~ /^(\d+)\s(\d+)\s(\d+)/);
		my $exon = $1;
		my $exonPos = $2;
		my $coverage = $3;				
		
		die "Pos twice in $file ? $exon $exonPos" if(defined $forReturn{$exon}{$exonPos});
		$forReturn{$exon}{$exonPos} = $coverage;
	}	
	close(F);
	
	return \%forReturn;
}

sub twoValidationAlleles_2_proper_names
{
	my $alleles_validation = shift;
	my $locus = shift;
	my $exon_sequences = shift;
	
	die unless($#{$alleles_validation} == 1);
	
	my @forReturn;
	
	$locus =~ s/HLA//;
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my @vA = split(/;/, $alleles_validation->[$aI]);
		REFALLELE: foreach my $validation_allele (@vA)
		{
			my $original_validation_allele = $validation_allele;
			
			if($validation_allele =~ /\:/)
			{
				$validation_allele = $locus . '*' . $validation_allele;
				# $validation_allele =~ s/g//;
				# my @components = split(/\:/, $validation_alleles),
				
				# die unless(length($validation_allele) >= 4);
				
				# $validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
				# $validation_allele = $locus . '*' . $validation_allele;
			
			}
			else
			{
				$validation_allele =~ s/g//;
				die unless(length($validation_allele) >= 4);
				
				$validation_allele = substr($validation_allele, 0, 2) . ':' . substr($validation_allele, 2);
				$validation_allele = $locus . '*' . $validation_allele;
			}
			my @validation_alleles = ($validation_allele);
			
			my @history_extensions = ($validation_allele);
			my $extensions = 0;
			my $notFoundFirstRound = 0;
			while( not exists $exon_sequences->{$validation_allele} )
			{
				$extensions++;
			
				$validation_allele .= ':01';
				push(@history_extensions, $validation_allele);
				if($extensions > 3)
				{
					$notFoundFirstRound = 1;
					last;
				}
				
			}
			
			if($notFoundFirstRound)
			{
				my $extensions = 0;
				while(scalar(grep {exists $exon_sequences->{$_}} @validation_alleles) == 0)
				{
					$extensions++;
					my $viMax = $#validation_alleles;
					for(my $vI = 0; $vI <= $viMax; $vI++)
					{
						my $a1 = $validation_alleles[$vI];
						$validation_alleles[$vI] .= ':01';
						push(@validation_alleles, $a1.':02');
					}
					
					$validation_allele .= ':01';
					if($extensions > 3)
					{
						# print Dumper([grep {$_ =~ /DQB1/} keys %$exon_sequences]);
						if($validation_allele eq $vA[$#vA])
						{
							push(@forReturn, '?');								
							warn Dumper("Can't identify (II) exon alleles for $alleles_validation->[$aI]", \@validation_alleles, $locus, $alleles_validation, $extensions, [(keys %$exon_sequences)[0 .. 10]], $validation_allele, $original_validation_allele, \@history_extensions);							
							next REFALLELE;
						}
						else
						{
							next REFALLELE;
						}
					}
				}		
				
				my @foundAlleles = grep {exists $exon_sequences->{$_}} @validation_alleles;
				die unless(scalar(@foundAlleles) > 0);
				push(@forReturn, $foundAlleles[0]);	
				last REFALLELE;				
				
			}
			else
			{
				push(@forReturn, $validation_allele);
				last REFALLELE;
			}
		}
	}
			
	return @forReturn;
}

sub list_comparison
{
	my $list1_aref = shift;
	my $list2_aref = shift;
	
	my %l1 = map {$_ => 1} @$list1_aref;
	my %l2 = map {$_ => 1} @$list2_aref;
	
	unless(scalar(keys %l1) == scalar(@$list1_aref))
	{
		die "List 1 non-unique elements";
	}
	
	unless(scalar(keys %l2) == scalar(@$list2_aref))
	{
		die "List 2 non-unique elements";
	}	
	
	my %combined = map {$_ => 1} (@$list1_aref, @$list2_aref);
	
	my @l1_exclusive;
	my @l2_exclusive;
	my $n_shared = 0;
	foreach my $e (keys %combined)
	{
		if($l1{$e} and $l2{$e})
		{
			$n_shared++;
		}
		elsif($l1{$e})
		{
			push(@l1_exclusive, $e);
		}
		elsif($l2{$e})
		{
			push(@l2_exclusive, $e);
		}
		else
		{
			die;
		}
				
	}
	
	die unless(($n_shared + scalar(@l1_exclusive) + scalar(@l2_exclusive)) == scalar(keys %combined));
	
	return ($n_shared, \@l1_exclusive, \@l2_exclusive);
}

sub twoClusterAlleles
{
	my $alleles_inference = shift;
	die Dumper("Don't have two inferred alleles", $alleles_inference) unless($#{$alleles_inference} == 1);
	
	my @forReturn;
	
	for(my $aI = 0; $aI < 2; $aI++)
	{
		my $inferred_alleles = $alleles_inference->[$aI];
		my @inferred_alleles = split(/;/, $inferred_alleles);
	
		$inferred_alleles[0] =~ s/^\s+//;
		
		push(@forReturn, $inferred_alleles[0]);
	}
	
	return @forReturn;
}


sub averageFractionAlignmentOK
{
	my $dir = shift;
	my $f = $dir . '/summaryStatistics.txt';
	die "File $f not there" unless(-e $f);
	
	my $fractionOK;
	open(STATISTICS, '<', $f) or die;
	while(<STATISTICS>)
	{
		my $line = $_;
		chomp($line);
		last if($line =~ /unpaired/);
		if($line =~ /Alignment pairs, average fraction alignment OK:\s+([\d\.]+)/)
		{
			die if(defined $fractionOK);
			$fractionOK = $1;
		}
	}
	close(STATISTICS);	
	
	die "Could not find fraction OK entry" unless(defined $fractionOK);
	return $fractionOK;
}


sub inferReadLength
{
	my $file = shift;
	open(ALIGNED, '<', $file) or die;
	my $firstLine = <ALIGNED>;
	chomp($firstLine);
	my @have_lengths;
	while(<ALIGNED>)
	{
		my $line = $_;
		die unless($line =~ /Aligned pair /);
		my @p1 = map {scalar(<ALIGNED>)} (1 .. 9);
		my @p2 = map {scalar(<ALIGNED>)} (1 .. 9);
		die unless($p1[0] =~ /Read /);
		die unless($p2[0] =~ /Read /);
		my $sequence_1 = $p1[7];
		my $sequence_2 = $p1[7];
		chomp($sequence_1);
		chomp($sequence_2);
		$sequence_1 =~ s/\s+//g;
		$sequence_2 =~ s/\s+//g;
		push(@have_lengths, length($sequence_1));
		push(@have_lengths, length($sequence_2));
		last if(scalar(@have_lengths) >= 20);
	}
	close(ALIGNED);	
	
	my $avg = int(sum(@have_lengths)/scalar(@have_lengths));
	return 2*$avg;
}




sub read_exon_sequences
{
	my $file = shift;
	my %r;
	open(F, '<', $file) or die "Cannot open $file";
	my $firstLine = <F>;
	chomp($firstLine);
	my @header_fields = split(/ /, $firstLine);
	die unless($header_fields[0] eq 'IndividualID');
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_fields = split(/ /, $line);
		my $allele = shift(@line_fields);
		$r{$allele} = join('', @line_fields);
	}
	close(F);
	
	return \%r;
}



sub find_exon_file { 
    my ($locus, $exon, $exon_folder) = @_;  

	$locus =~ s/HLA//;
	
	opendir(my $dh, $exon_folder) or die;
	my @exon_folder_files = readdir($dh);
	closedir($dh);
	
	@exon_folder_files = grep {$_ =~ /^${locus}_\d+_exon_${exon}\.txt$/} @exon_folder_files;
	if($#exon_folder_files != 0)
	{
		die Dumper("Can't find exon file for $locus // $exon", @exon_folder_files);
	}
	else
	{
		return $exon_folder.'/'.$exon_folder_files[0];
	}	
}





sub min_avg_max
{
	my @v = @_;
	@v = sort {$a <=> $b} @v;
	if($#v >= 1)
	{
		die unless($v[0] <= $v[1]);
	}
	
	if(scalar(@v) == 0)
	{
		return ('', '', '');
	}
	if(scalar(@v) == 1)
	{
		return($v[0], $v[0], '');
	}
	
	my $min = $v[0];
	my $max = $v[$#v];
	  
	my $sum = 0;
	foreach my $vE (@v)
	{
		$sum += $vE;
	}
	
	my $avg = '';
	if(scalar(@v) > 0)
	{
		$avg = $sum / scalar(@v);
	}
	
	return ($min, $avg, $max);
}
