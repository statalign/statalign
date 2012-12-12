#!/usr/bin/perl -w
my $MAX_SET_SIZE = 10;

my $pir = shift;
my $MAX_PERCENT_ID = 25;
if (@ARGV) {
    $MAX_PERCENT_ID = shift;
}
my $MAX_ID = $MAX_PERCENT_ID/100;
my $structalign_dir = "~/workspace/structalign/";

my $nseqs = 0;
my $n;
my @ids;
my $id;
my (@pdb,@start,@end,%chain,@code,%seq);

my $aliname = $pir; 
$aliname =~ s/\.pir//;
open(PIR, "< $pir");
while(<PIR>) {

    $n = $nseqs - 1;
    if (/structureX:(\w{4}):\s+(\d+)\s*:([A-Z\s]):\s*(\d+)\s*:/) {
	
	$pdb[$n] = $1;
	$chain[$n] = $3;
	my $lower_chain = $chain[$n];
	$lower_chain =~ tr/[A-Z]/[a-z]/;
	$code[$n] = $pdb[$n];
	if ($lower_chain && $lower_chain ne " ") {
	    $code[$n] .= $lower_chain;
	}
	$start{$code[$n]} = $2;
	$end{$code[$n]} = $4;
	push(@ids,$code[$n]);
	next;
    }
    elsif (/>P1;(\w+)/) {
	$nseqs++;
	next;
    }
    chomp;
    s/\*//;
    $seq{$code[$n]} .= $_; 
}
close PIR;

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
	die "Cannot compute sequence identity between two sequences of different length.\n";
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
	    $average_id[$i] += seq_id($seq{$id_list->[$i]},$seq{$id_list->[$j]})/$n;
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
	my $chain = "";
	if (length($id_list->[$i]) == 5) {
	    $chain = substr($id_list->[$i],4,1);
	    $chain =~ tr/[a-z]/[A-Z]/;
	}
	print FASTA ">$id_list->[$i]\n";	
	print COOR ">$id_list->[$i]\n";	
	print FASTA $seq{$id_list->[$i]}."\n";

	my $coorfile = $id_list->[$i].".coor";

	my $pdb_id = substr($id_list->[$i],0,4);
	`wget http://www.rcsb.org/pdb/files/$pdb_id.pdb 2> /dev/null`;
	if ($chain{$id_list->[$i]} && $chain{$id_list->[$i]} ne " ") {
	    `$structalign_dir/scripts/selecta.pl miss[]b[A]c[$chain{$id_list->[$i]}]]coor[]r[$start{$id_list->[$i]}-$end{$id_list->[$i]}]o[$coorfile] $pdb_id.pdb`;
	}
	else {
	    `$structalign_dir/scripts/selecta.pl miss[]b[A]coor[]r[$start{$id_list->[$i]}-$end{$id_list->[$i]}]o[$coorfile] $pdb_id.pdb`;
	}
	`rm $pdb_id.pdb`;

	open(THIS_COOR, "< $coorfile");
	my $this_coor = "";
	while(<THIS_COOR>) {
	    $this_coor .= $_;
	}
	close THIS_COOR;
	`rm $coorfile`;
	print COOR $this_coor;
    }
}
