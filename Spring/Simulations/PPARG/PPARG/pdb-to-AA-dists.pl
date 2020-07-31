#!/usr/bin/perl -w

if(@ARGV==0) { system "$0 PPARG-3dzu.pdb D 28" ; exit }

%AAA2A=('ALA' => 'A', 'VAL' => 'V', 'LEU' => 'L', 'ILE' => 'I', 'PRO' => 'P', 'TRP' => 'W', 'PHE' => 'F', 'MET' => 'M', 'GLY' => 'G', 'SER' => 'S', 
        'THR' => 'T', 'TYR' => 'Y', 'CYS' => 'C', 'ASN' => 'N', 'GLN' => 'Q', 'LYS' => 'K', 'ARG' => 'R', 'HIS' => 'H', 'ASP' => 'D', 'GLU' => 'E');

die "Usage: $0 pdb chain\n" if @ARGV < 2;

$myC   = $ARGV[1];
$SHIFT = @ARGV > 2 ? $ARGV[2] : 0;

open(IN,$ARGV[0]) || die;
while(<IN>) { 
	next if !/^ATOM/ && !/^HETATM/;
	chomp($line=$_);
	if(/^ATOM/)      { $id=substr $line, 21, 1 }
	elsif(/^HETATM/) { $id=substr $line, 21, 5 }
	next if $myC ne $id;
	push @r,        substr $line, 17, 4;
	push @n, $SHIFT+substr $line, 22, 4; 
	push @x,      0+substr $line, 30, 8;
	push @y,      0+substr $line, 38, 8;
	push @z,      0+substr $line, 46, 8;
}

foreach(@r) { $_=~s/ //g  }
foreach(@n) { $_=~s/ //g  }

for($j=$i=0;$i<@n;$i++) { 
	if(!exists $n{$n[$i]}) { 
		die  $r[$i] if !exists $AAA2A{$r[$i]};
		push @myseq, $AAA2A{$r[$i]};
		push @mynum, $n[$i];
		push @resX,  $x[$i];
		push @resY,  $y[$i];
		push @resZ,  $z[$i];
		push @resN,  1;
		$n{$n[$i]} = $j++;
	}
	else {
		$k = $n{$n[$i]};
		$resX[$k] += $x[$i];
		$resY[$k] += $y[$i];
		$resZ[$k] += $z[$i];
		$resN[$k] += 1;
	}
}

for($i = 0; $i<@resN; $i++) {
	$resX[$i] /= $resN[$i];
	$resY[$i] /= $resN[$i];
	$resZ[$i] /= $resN[$i];
}

print "Res";
for($i=0;$i<@resN;$i++) { print "\t$myseq[$i]$mynum[$i]" } print "\n";
for($i=0;$i<@resN;$i++) {
	$x=$resX[$i];
	$y=$resY[$i];
	$z=$resZ[$i];
	print $myseq[$i].$mynum[$i];
	for($j=0;$j<@resN;$j++) {
		$dist=sqrt(($x-$resX[$j])*($x-$resX[$j])+($y-$resY[$j])*($y-$resY[$j])+($z-$resZ[$j])*($z-$resZ[$j]));
		print sprintf "\t%.2e",$dist;
	}
	print "\n";
}

