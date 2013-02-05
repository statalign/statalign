#!/usr/bin/perl -w
my $MAX_SET_SIZE = 30;

my $pir = shift;
my $MAX_PERCENT_ID = 25;
if (@ARGV) {
    $MAX_PERCENT_ID = shift;
}
my $MAX_ID = $MAX_PERCENT_ID/100;
my $structalign_dir = '~/private/workspace/statalign';

my $nseqs = 0;
my $n;
my @ids;
my $id;
my (@pdb,%start,%end,%chain,@code,%seq);

my $aliname = $pir; 
$aliname =~ s/\.pir//;
open(PIR, "< $pir");
while(<PIR>) {

    $n = $nseqs - 1;
    if (/structureX:(\w{4}):\s+(\d+)\s*:([A-Z\s]):\s*(\d+)\s*:/) {
	#print;
	$pdb[$n] = $1;
	$chain[$n] = $3;
	my $lower_chain = $chain[$n];
	$lower_chain =~ tr/[A-Z]/[a-z]/;
	$code[$n] = $pdb[$n];
	if ($lower_chain && $lower_chain ne " ") {
	    $code[$n] .= $lower_chain;
	}
	while (exists $start{$code[$n]}) {
	    $code[$n] .= "x";
	}
	$start{$code[$n]} = $2;
	$end{$code[$n]} = $4;
	push(@ids,$code[$n]);
	next;
    }
    elsif (/^>/) {
	$nseqs++;
	next;
    }
    chomp;
    s/\*//;
    $seq{$code[$n]} .= $_; 
}
close PIR;

#print "Number of sequences = $nseqs.\n";

my @average_ids = compute_average_identity(\@ids);

my @selected_set;
my @remaining_set = sort {$average_ids[$a] <=> $average_ids[$b] } (0..$nseqs-1);
#foreach (@remaining_set) { print; print " ";} print "\n";

SET_SIZE: if (@selected_set < $MAX_SET_SIZE) {
    #foreach (@remaining_set) { print; print " ";} print "\n";
    REMAINING: for (my $i=0; $i<@remaining_set; $i++) {
	#print $ids[$i]."\n";
	if (@selected_set) {
	    foreach $id (@selected_set) {
		my $seq_id = seq_id($seq{$ids[$remaining_set[$i]]},$seq{$id});
		#print "...".$id."\t".$seq_id."\n";
		if ($seq_id > $MAX_ID) {
		    next REMAINING;
		}
	    }
	}
	push(@selected_set,$ids[$remaining_set[$i]]);
	splice(@remaining_set,$i,1);
	goto SET_SIZE;
    }  
}

#print_identity_matrix(\@selected_set);
print_alignment(\@selected_set);

sub seq_id {
    
    my $seq1 = $_[0];
    my $seq2 = $_[1];

    if (length($seq1) != length($seq2)) {
	die "Cannot compute sequence identity between two sequences of different length.\n length(seq1) = ".length($seq1).", length(seq2) = ".length($seq2)."\n";
    }
    my $identity = 0;
    for (my $i=0; $i<length($seq1); $i++) {
	 $identity += (substr($seq1,$i,1) eq substr($seq2,$i,1));
    }
    return($identity/length($seq1));
}

sub compute_identity_matrix {

    my $id_list = $_[0];
    my @identity; 

    my $n = @{$id_list};
    for (my $i=0; $i<$n; $i++) {
	$identity[$i] = [(0) x $n];
    }
    for (my $i=0; $i<($n-1); $i++) {
	$$identity[$i][$i] = 1;
	for (my $j=$i+1; $j<$n; $j++) {
	    $$identity[$i][$j] = seq_id($seq{$id_list->[$i]},$seq{$id_list->[$j]});
	    $$identity[$j][$i] = $$identity[$i][$j];
	}
    }
    $$identity[$n-1][$n-1] = 1;
    return @identity;
}

sub compute_average_identity {

    my $id_list = $_[0];
    my @average_id; 

    my $n = @{$id_list};
    for (my $i=0; $i<$n; $i++) {
	$average_id[$i] = 0;
	for (my $j=0; $j<$n; $j++) {
	    next if ($i==$j);
	    $average_id[$i] += seq_id($seq{$id_list->[$i]},$seq{$id_list->[$j]})/$n
	}
    }
    return @average_id;
}

sub print_identity_matrix {

    my $id_list = $_[0];
    my @identity = compute_identity_matrix($id_list);
    my $n = @{$id_list};
    for (my $i=0; $i<$n; $i++) {
	printf "%8s", $id_list->[$i];
	for (my $j=0; $j<$n; $j++) {
	    printf "%8.3f", $$identity[$i][$j] if ($$identity[$i][$j]);
	}
	print "\n";
    }
}

sub print_alignment {

    my $id_list = $_[0];
    my $n = @{$id_list};
    my $filename = "$aliname"."_$MAX_PERCENT_ID";
    open(FASTA, "> $filename.fasta");
    open(COOR, "> $filename.coor");
    for (my $i=0; $i<$n; $i++) {
	#print $i."\n";
	my $chain = "";
	if (length($id_list->[$i]) == 5) {
	    $chain = substr($id_list->[$i],4,1);
	    $chain =~ tr/[a-z]/[A-Z]/;
	}
	my $coorfile = $id_list->[$i].".coor";

	my $pdb_id = substr($id_list->[$i],0,4);
	print STDERR $pdb_id."\n";
	`wget http://www.rcsb.org/pdb/files/$pdb_id.pdb 2> /dev/null`;
	my $selecta_string = "$structalign_dir/scripts/selecta.pl miss[]b[A]coor[]o[$coorfile]a[]";
	if ($start{$id_list->[$i]}) {
	    $selecta_string .= "r[$start{$id_list->[$i]}-$end{$id_list->[$i]}]";
	}
	else {
	    my $sequence = $seq{$id_list->[$i]};
	    $sequence =~ s/-//g;
	    open(SEQ, "> tmp.fasta");
	    print SEQ $sequence."\n";
	    close SEQ;
	    $selecta_string .= "seq[tmp.fasta]";
	}
	if ($chain && $chain ne " ") {
	    $selecta_string .= "c[$chain]";
	}
	my $errors = `$selecta_string $pdb_id.pdb`;
	`rm $pdb_id.pdb`;
	if ($errors) {
	    print STDERR $errors;
	    print STDERR "Skipping this sequence.\n";
	}
	else {
	    print FASTA ">$id_list->[$i]\n";	
	    print COOR ">$id_list->[$i]\n";	
	    print FASTA $seq{$id_list->[$i]}."\n";
	    open(THIS_COOR, "< $coorfile");
	    my $this_coor = "";
	    while(<THIS_COOR>) {
		$this_coor .= $_;
	    }
	    close THIS_COOR;
	    print COOR $this_coor;
	}
	`rm $coorfile`;
    }
}
