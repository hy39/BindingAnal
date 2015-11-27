#!/usr/bin/perl
#Read fasta sequences, format the age information and save it.
#Aug 30, 2012
#Hsiang-Yu Yuan

my $filename = $ARGV[0];
my $out_filename = "hm_h1n1_flu_whole_age.csv";

# specify the filename if user didn't give it
if ($filename eq '') {
 $filename = "seq/hm_h3n2/hm_h3n2_flu_ny_any.csv";
}

open FILE, "$filename" or die $!;

print "reading file: $filename\n";

my %labelhash = ();
my @tokens;
my @results;
my $sequence;

my $i=0;
while (my $line = <FILE>) { 
       @tokens = split(/,/, $line);
       #print "$tokens[0] $tokens[9]\n";
       
       $tokens[0] =~ m/"(.*)"/; #gb access 
       my $gbacc = $1;
       $tokens[9] =~ m/"(.*)"/;
       my $age = $1;
       if ($age =~ m/([0-9]*)(m|M)/) {
         $age = 1;
       }
       if ($age =~ m/([0-9]*)(y|Y)/) {
         $age = $1;
       }
       
       $results[$i][0] = $gbacc;
       $results[$i][1] = $age;
       $i++;     
}
   
    
print "write to file: seq/$out_filename";
open OUTFILE, ">seq/$out_filename";   
print OUTFILE "GBACC,AGE\n";
for (my $j = 0; $j < scalar(@results); $j++) { 
       my $out_str = "\"$results[$j][0]\",$results[$j][1]";
       print OUTFILE "$out_str\n";
}
close OUTFILE;

 
