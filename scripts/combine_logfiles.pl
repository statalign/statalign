#!/usr/bin/perl -w

my @samples = ();

sub read_logfile {
    my $log = $_[0];
    open(LOG, "< $log"); 
    while(<LOG>) {
	if (/^Sample (\d+)\tAlignment:\t[^->\w]*([->\w]+)$/) { 
	    $samples[$1] .= $_;
	}
    }
    close LOG;
}

foreach $logfile (@ARGV) {
    read_logfile($logfile);
}

my $combined_file = $ARGV[0];
$combined_file =~ s/chain(\d+)/combined/;
open(COMBINED, "> $combined_file");
foreach $sample (@samples) {
    last unless (defined $sample);
    print COMBINED "$sample";
}
close COMBINED;
