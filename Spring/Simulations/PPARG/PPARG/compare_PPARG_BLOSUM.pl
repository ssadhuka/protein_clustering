#!/usr/bin/perl -w

open(IN,"PPARG-scores.txt") || die;
chomp($header = <IN>);
while(<IN>) {
	chomp;
	($r,$m,$s) = split /\t/, $_;
	$r =~ s/\d+//;
	$a2a2s{$r}{$m} += $s; 
	$a2a2s{$m}{$r} += $s; 
	$a2a2n{$r}{$m} +=  1;
	$a2a2n{$m}{$r} +=  1; 
}

foreach $a1 (keys %a2a2s) {
	foreach $a2 (keys %{$a2a2s{$a1}}) {
		$a2a2s{$a1}{$a2} /= $a2a2n{$a1}{$a2};
	}
}

open(IN,"BLOSUM.txt") || die;
chomp($header = <IN>);
@A = split /\t/, $header;
shift @A;
while(<IN>) {
	chomp;
	@scores = split /\t/,$_;
	$A = shift @scores;
	for($i=0; $i<@A; $i++) {
		$b2b2s{$A}{$A[$i]} = $scores[$i];
	}
}

for($i=0; $i<@A; $i++) {
	$a1 = $A[$i];
	for($j=0; $j<@A; $j++) {
		next if $i == $j;
		$a2 = $A[$j];
		$s1 = $a2a2s{$a2}{$a1};
		$s2 = $b2b2s{$a2}{$a1};
		print "$s1\t$s2\n";
	}
}
