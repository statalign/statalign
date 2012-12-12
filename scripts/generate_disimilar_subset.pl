#!/usr/bin/perl -w
my $MAX_SET_SIZE = 10;
my $MAX_ID = 0.35;

my $nseqs = 0;
my @ids;
my $id;
while(<>) {

    if (/>(\w+)/) {
	$id = $1;
	push(@ids,$id);
	$nseqs++;
	next;
    }
    chomp;
    $seq{$id} .= $_; 
}

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
    for (my $i=0; $i<$n; $i++) {
	print ">$id_list->[$i]\n";
	print $seq{$id_list->[$i]}."\n";
    }
}
