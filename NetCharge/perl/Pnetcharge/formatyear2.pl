#!/usr/bin/perl
#Read fasta sequences, format the age information and save it.
#Aug 30, 2012
#Hsiang-Yu Yuan

my $filename = $ARGV[0];
my $out_filename = "hm_h1n1_whole_any_ay.csv";

# specify the filename if user didn't give it
if ($filename eq '') {
# $filename = "seq/hm_h3n2/hm_h3n2_flu_ny_any.csv";
 $filename = "seq/hm_h3n2_sea/hm_h3n2_sea_any.csv";
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
       
       #extract gb accession
       $tokens[0] =~ m/"(.*)"/; #gb access 
       my $gbacc = $1;
       
       #extract age
       $tokens[9] =~ m/"(.*)"/;
       my $age = $1;
       if ($age =~ m/([0-9]*)(d|D)/) {
         $age = 1;
       }
       if ($age =~ m/([0-9]*)(m|M)/) {
         $age = 1;
       }
       if ($age =~ m/([0-9]*)(y|Y)/) {
         $age = $1;
       }
       if ($age =~ m/N.*/) {
         next;
       }
 
       #extract year
       $tokens[7] =~m/"(.*)"/;
       my $year = $1;
       if ($year =~ m/(\d\d\d\d\/\d\d\/\d\d)/) {
         $year = $1;
       } else {
         next;
       }
      
 

       $results[$i][0] = $gbacc;
       $results[$i][1] = $age;
       $results[$i][2] = $year;
       $i++;     
}
   
    
print "write to file: seq\/$out_filename";
open OUTFILE, ">seq\/$out_filename";   
print OUTFILE "GBACC,AGE,ISODATE\n";
for (my $j = 0; $j < scalar(@results); $j++) { 
       my $out_str = "$results[$j][0],$results[$j][1],$results[$j][2]";
       print OUTFILE "$out_str\n";
}
close OUTFILE;

 
