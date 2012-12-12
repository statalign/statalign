#!/usr/bin/perl -w

my $nseqs = 0;
my @ids;
my $id;

my $pir = shift;

my (@pdb,@start,@end,@chain,@code,%seq);

open(PIR, "< $pir");
my $n=-1;
while(<PIR>) {
    
    #next if (/^>/);

    if (/structureX:(\w{4}):\s+(\d+)\s*:([A-Z\s]):\s*(\d+)\s*:/) {
	
	$n++;
	$pdb[$n] = $1;
	$start[$n] = $2;
	$chain[$n] = $3;
	my $lower_chain = $chain[$n];
	$lower_chain =~ tr/[A-Z]/[a-z]/;
	$code[$n] = $pdb[$n].$lower_chain;
	$end[$n] = $4;
	
	`wget http://www.rcsb.org/pdb/files/$pdb[$n].pdb`;
	if ($chain[$n] && $chain[$n] ne " ") {
	    print "./selecta.pl miss[]b[A]c[$chain[$n]]coor[]r[$start[$n]-$end[$n]]o[$code[$n].coor] $pdb[$n].pdb\n";
	    `./selecta.pl miss[]b[A]c[$chain[$n]]coor[]r[$start[$n]-$end[$n]]o[$code[$n].coor] $pdb[$n].pdb`;
	}
	else {
	    print "./selecta.pl miss[]b[A]coor[]r[$start[$n]-$end[$n]]o[$code[$n].coor] $pdb[$n].pdb\n";
	    `./selecta.pl miss[]b[A]coor[]r[$start[$n]-$end[$n]]o[$code[$n].coor] $pdb[$n].pdb`;
	}
	`rm $pdb[$n].pdb`;
	next;
    }
    else if (/>(\w+)/) {
	$id = $1;
	push(@ids,$id);
	$nseqs++;
	next;
    }
    chomp;
    $seq{$id} .= $_; 
}
close PIR;
