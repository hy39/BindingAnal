#!/usr/bin/perl
#Read fasta sequences and save each sequence into a separate seq
#Jun 24, 2011
#Hsiang-Yu Yuan

my $filename = $ARGV[0];
open FILE, "$filename" or die $!;
#my @lines = <FILE>;

my $prev_filename;
my $sequence;

while (my $line = <FILE>) {
    if($line =~ m/>(\w+)/) {
	my $fa = $1;
        #print "$fa\n";
        #print $sequence;
        if(defined($prev_filename)) {
            print "$prev_filename\n";
            print "$sequence\n";
 
	    open OUTFILE, ">seq/$prev_filename";
            print OUTFILE $sequence;
            close OUTFILE;
	} else {
        }
        $prev_filename = $fa . ".seq";
        $sequence = "";
    } else {
        $sequence = $sequence . $line;
#        print $line;
    }
}
    print "$prev_filename\n";
    print "$sequence\n";

    open OUTFILE, ">seq/$prev_filename";
    print OUTFILE $sequence;
    close OUTFILE;

 
