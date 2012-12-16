#!/usr/bin/perl -w

my $CALPHA = 0;
my $ALTLOC = "";
my $CHAIN = "";
my $NOHETATM = 1; # by default filter out heteroatoms
my $ATOMONLY = 0;
my $COOR_ONLY = 0;
my $NOMISSING = 0; # Flags if there are any missing residues
my $MODEL = 0;
my $range = 0;
my @range = [1,100];

my $options = shift @ARGV;

#print STDERR $options."\n";

my $output = "";

if ($options =~ /a\[\]/) {

    $CALPHA++;
}
if ($options =~ /miss\[\]/) {

    $NOMISSING++;
}
if ($options =~ /coor\[\]/) {

    $COOR_ONLY++;
    $ATOMONLY++;
    $NOHETATM = 0;
}
if ($options =~ /b\[(\w)\]/) {

    $ALTLOC = $1;
    # The [A] option includes the first observed location for each residue 
    # with multiple locations.
    # [B] will print only the second, but only if the ALTLOC column contains 
    # a "B", otherwise the first AND the second will be printed.
}
if ($options =~ /c\[(\w)\]/) {

    $CHAIN = $1;
}
if ($options =~ /h\[(\d)\]/) {

    $NOHETATM = $1;
}
if ($options =~ /at\[\]/) {

    $ATOMONLY++;
    $NOHETATM = 0;
}
if ($options =~ /mod\[(\d+)\]/) {

    $MODEL = $1;
}
if ($options =~ /o\[(\S+?)\]/) {

    $output = $1;
}
if ($options =~ /r\[(\d+)-(\d+)\]/) {

    $range++;
    $range[0] = $1;
    $range[1] = $2;
}

#print STDERR $output."\n";

my $pdb = shift @ARGV;
unless($output) {
    $output .= $pdb;
    $output =~ s/\.pdb/\.mod/;
}
my $atoms_printed = 0;
my $res_printed = 1; # will be 1 even if actually 0 are printed, for purposes of counting.
my $prev_res = $range[0]-1; # keeps track of duplicate residues
my $prev_atom = 0; # keeps track of duplicate residues
my $curr_alt_loc = "A";
my $model = 0;

my $line_num = 0;

open(PDB, "< $pdb");
open(OUTPUT, "> $output");
# if ($COOR_ONLY) {
#     my $code = $pdb;
#     $code =~ s/\.pdb//;
#     print OUTPUT ">$code\n";
# }
while (<PDB>) {

    my $line = $_;

    my ($atom_n, $atom_name, $alt_loc, $res_name, $chain, $res_n, $ins, $x, $y, $z, $occupancy, $b_factors, $footnote, $seg_id, $element, $charge);

    if ($line =~ /^MODEL\s+(\d+)/) { 

	$model = $1; 
	$line_num = 0;
    }

    # if no MODEL lines in PDB file, then $model = 0 and $MODEL = 0 always

    if ($line =~ /^ATOM/) {

	#print "1:\n";
	#print;
	$line_num++;
	
	if (($MODEL == 0)&&($model > 0)) { 

	    $MODEL = $model;
	    #print STDERR "Only one model present in $pdb therefore selection of $MODEL ignored.\n";
	}
	next unless ($model == $MODEL); 

	$atom_n = substr($line, 6, 5);
	$atom_name = substr($line, 12, 4);
	$alt_loc = substr($line, 16, 1);
	$res_name = substr($line, 17, 3);
	$chain = substr($line, 21, 1);
	$res_n = substr($line, 22, 4);
	$ins = substr($line, 26, 1);
	$x = substr($line, 30, 8);
        $y = substr($line, 38, 8);
	$z = substr($line, 46, 8);

	if (length ($line) > 54) {

	    $occupancy = substr($line, 54, 6);
	    
	    if (length ($line) > 60) {

		$b_factors = substr($line, 60, 6);
		
		if (length ($line) > 66) {

		    $footnote = substr($line, 66, 6);
		    
		    	if (length ($line) > 72) {
			    
			    $seg_id = substr($line, 72, 4);

			    	if (length ($line) > 76) {

				    $element = substr($line, 76, 2);
				    
				    	if (length ($line) > 78) {

					    $charge = substr($line, 78, 2);
					}
				}
			}
		}
	    }
	}
	#print "2:\n";
	#print;
	if ($CHAIN) {

	    next unless (($chain eq $CHAIN)||($chain eq " "));
	}
	#print "3:\n";
	#print;
	if ($CALPHA) {

	    next unless ($atom_name =~ /CA/);
	}
	#print "4:\n";
	#print;
	if ($ALTLOC) {

	    if ($atom_n != $prev_atom) {	

		$curr_alt_loc = "A";
	    }
	    else {

		$curr_alt_loc++;
		$alt_loc = $curr_alt_loc;
	    }
	    if ($alt_loc =~ /\S/) {

		next unless ($alt_loc eq $ALTLOC);
	    }
	}
	#print "5:\n";
	#print;
	if ($range) {

	    next unless (($res_n >= $range[0])&&($res_n <= $range[1]));
	}
	#print "6:\n";
	#print;
	if ($NOMISSING) {

	    if ((($res_n - $prev_res)>1)||(($line_num == 1)&&($res_n>$prev_res))) {
		print "Missing residues detected. Gap between $prev_res and $res_n\n";
	    }
	}
	#print "7:\n";
	#print;
	if (($res_n - $prev_res)==1) {

	    $res_printed++;
	}

	$prev_atom = $atom_n;
	$prev_res = $res_n;
	
	if ($COOR_ONLY) {
	    print OUTPUT ($x+0.0)."\t".($y+0.0)."\t".($z+0.0)."\n";
	    #print ($x+0.0)."\t".($y+0.0)."\t".($z+0.0)."\n";
	}
	else {
	    print OUTPUT $line;
	}
	$atoms_printed++;

    }
    else {

	next if ($ATOMONLY);
	next unless ($model == $MODEL); 
	if ($NOHETATM) {

	    #print STDERR "HETATM skipped\n";
	    next if ($line =~ /^HETATM/);
	}
	print OUTPUT $line;
    }
}
close PDB;
close OUTPUT;

if ($range && ($prev_res != $range[1])) {

    print "Expected residue range not observed. Chain terminated at residue $prev_res instead of $range[1].\n";
}

#print STDERR "Number of atoms printed = ".$atoms_printed."\n";
print "Number of residues printed =               ".$res_printed."\n\n";
