#!/usr/bin/perl -w

%AAA2A=('ALA' => 'A', 'VAL' => 'V', 'LEU' => 'L', 'ILE' => 'I', 'PRO' => 'P', 'TRP' => 'W', 'PHE' => 'F', 'MET' => 'M', 'GLY' => 'G', 'SER' => 'S', 
        'THR' => 'T', 'TYR' => 'Y', 'CYS' => 'C', 'ASN' => 'N', 'GLN' => 'Q', 'LYS' => 'K', 'ARG' => 'R', 'HIS' => 'H', 'ASP' => 'D', 'GLU' => 'E');

$myC   = "D";
$SHIFT =  28;

#rewrite PDB file
open(IN,"PPARG-3dzu.pdb")         || die;
open(OT,">PPARG-renum-3dzuD.pdb") || die;
$indx = 0;
while(<IN>) { 
	next if !/^ATOM/ && !/^HETATM/;
	chomp($line=$_);
	if(/^ATOM/)      { $id=substr $line, 21, 1 }
	elsif(/^HETATM/) { $id=substr $line, 21, 5 }
	next if $myC ne $id;
	$aaa = substr $line, 17, 4;
	$n   = $SHIFT+substr $line, 22, 4;
	$aaa =~s/ //g; 
	$n   =~s/ //g;
	$res = "$AAA2A{$aaa}$n";
	if(not exists $done{$res}) {
		$r2i{$res}  = $indx++;
		$done{$res} = 1;
		push @r, $res;
	}
	print OT substr $line, 0, 22;
	print OT sprintf "%4s", $n;
	print OT substr $line, 26;
	print OT "\n";
}

#rewrite score file
$header = "";
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

open(OT,">PPARG-betalike.txt") || die;
foreach($i=0; $i<@r; $i++) {
	print OT "p.$r[$i]X\t$score[$i]\t0\t0_$i\_X_X\tPPARG\n";
}
