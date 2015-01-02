#!/usr/bin/perl

my @samples = ('I1_AA02O9Q_A1', 'I1_AA02O9Q_A3');

foreach my $sample (@samples)
{
	my %origins;
	foreach my $oneTwo (1, 2)
	{
		my $fn = '../tmp/hla/'.$sample.'/reads.p.n_'.$oneTwo;
		my $fh;
		open($fh, '<', $fn) or die "Cannot open $fn";
		while(<$fh>)
		{
			my $line = $_;
			chomp($line);
			die unless(substr($line, 0, 1) eq '@');
			
			if($line =~ /FROM:(.+?):(\d+)/)
			{
				my $chromosome = $1;
				my $position = $2;
				$origins{$chromosome}++;
			}
			
			getlines($fh, 3);		
		}
		close(F);
	}
	
	print "Sample $sample\n";
	my @origins = sort {$origins{$b} <=> $origins{$a}} keys %origins;
	my $S = 0;
	for(@origins){$S += $origins{$_};}
	foreach my $origin (@origins)
	{
		my $perc = sprintf("%.2f", $origins{$origin} / $S * 100).'%';
		print "\t", $origin, "\t", $perc, "\t", $origins{$origin}, "\n";
	}	
}

sub getlines
{
	my $fh = shift;
	my $n = shift;
	my @forReturn;
	for(my $i = 0; $i < $n; $i++)
	{
		my $l = <$fh>;
		die unless(defined $l);
		push(@forReturn, $l);
	}
	return @forReturn;
}