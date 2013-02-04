#!/usr/bin/perl -w

my $fasta = shift;
open(FASTA, "< $fasta") or die "$!: Could not open file $fasta.\n";
my @tmp = split(/\./,$fasta);
my $pir = $tmp[0].".pir";
open(PIR, "> $pir");
while(<FASTA>) {

    my @range;
    my @name;
    my $pdb_code;
    my $chain;
    if (/^>(.+)$/) {
	@name = split(/\//,$1);

	if (@name > 1) {
	    @range = split(/-/,$name[1]);
	}
	if (length($name[0]) == 7) {
	    $pdb_code = substr($name[0],1,4);
	    $chain = substr($name[0],5,1);
	    $chain =~ s/_/ /;
	    $chain =~ tr/[a-z]/[A-Z]/;
	}
	
	print PIR ">$name[0]\n";
	print PIR "structureX:$pdb_code: $range[0] :$chain: $range[1] : : : : : \n";
	next;
    }
    print PIR;
}
close FASTA;
close PIR;
