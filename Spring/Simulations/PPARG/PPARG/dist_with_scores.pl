#!/usr/bin/perl -w
use Statistics::Basic qw(:all nofill);

#read dists
open(IN,"PPARG-AA-dists.txt") || die;
chomp($header = <IN>);
@header = split /\t/, $header;
shift @header;
$i = 0; foreach(@header) { $r2i{$_} = $i++ } 
$i = 0;
while(<IN>) {
	chomp;
	@d=split /\t/, $_;
	shift @d;
	push @{$dist[$i++]}, @d;
}

#read scores
open(IN,"PPARG-scores.txt") || die;
chomp($header = <IN>);
while(<IN>) {
	chomp;
	($r,$i,$s) = split /\t/, $_;
	next if not exists $r2i{$r};
	$score[$r2i{$r}] += $s;
	$count[$r2i{$r}] +=  1;
}
for($i=0; $i<@score; $i++) { $score[$i]/=$count[$i] }

die @score+0,"\t",@header+0,"\n" if @score != @header;

for($dist = 4; $dist < -15; $dist+=0.1) {
	@myscore = ();
	@prscore = ();
	for($i=0; $i<@score; $i++) {
		$score = 0;
		$count = 0;
		for($j = 0; $j<@score; $j++) {
			next if $i == $j;
			if($dist[$i][$j] < $dist) {
				$score += $score[$j];
				$count += 1;
			}
		}
		next if $count == 0;
		$score /= $count;
		push @myscore, $score[$i];
		push @prscore, $score;
	}
	next if @myscore == 0;
	$rmsd = 0;
	for($i=0; $i<@myscore; $i++) { $rmsd += ($myscore[$i]-$prscore[$i])**2 }
	$rmsd = sqrt($rmsd/@myscore);
	$v1 = vector(@myscore);
	$v2 = vector(@prscore);
	$corr = sprintf "%.5f", corr($v1,$v2);
	$dist = sprintf "%.2f", $dist;
	print "$dist\t$corr\t$rmsd\n";
}

print "\n\n";

for($t = 0.1; $t < -2; $t+=0.01) {
	@myscore = ();
	@prscore = ();
	for($i=0; $i<@score; $i++) {
		$score = 0;
		$count = 0;
		for($j = 0; $j<@score; $j++) {
			next if $i == $j;
			$weight = exp(-$dist[$i][$j]**$t);
			$score += $weight*$score[$j];
			$count += $weight;
		}
		next if $count == 0;
		$score /= $count;
		push @myscore, $score[$i];
		push @prscore, $score;
	}
	next if @myscore == 0;
	$rmsd = 0;
	for($i=0; $i<@myscore; $i++) { $rmsd += ($myscore[$i]-$prscore[$i])**2 }
	$rmsd = sqrt($rmsd/@myscore);
	$v1 = vector(@myscore);
	$v2 = vector(@prscore);
	$corr = sprintf "%.5f", corr($v1,$v2);
	$t    = sprintf "%.2f", $t;
	print "$t\t$corr\t$rmsd\n";
}

print "\n\n";

for($dist = 5; $dist < 30; $dist+=1) {
	for($t = 0.1; $t < 1; $t+=0.1) {
		@myscore = ();
		@prscore = ();
		for($i=0; $i<@score; $i++) {
			$score = 0;
			$count = 0;
			for($j = 0; $j<@score; $j++) {
				next if $i == $j;
				if($dist[$i][$j] < $dist) {
					$weight = exp(-$dist[$i][$j]**$t);
					$score += $weight*$score[$j];
					$count += $weight;
				}
			}
			next if $count == 0;
			$score /= $count;
			push @myscore, $score[$i];
			push @prscore, $score;
		}
		next if @myscore == 0;
		$rmsd = 0;
		for($i=0; $i<@myscore; $i++) { $rmsd += ($myscore[$i]-$prscore[$i])**2 }
		$rmsd = sqrt($rmsd/@myscore);
		$v1 = vector(@myscore);
		$v2 = vector(@prscore);
		$corr = sprintf "%.5f", corr($v1,$v2);
		$dist = sprintf "%.2f", $dist;
		$t    = sprintf "%.2f", $t;
		print "$dist\t$t\t$corr\t$rmsd\n";
	}
	print "\n";
}
