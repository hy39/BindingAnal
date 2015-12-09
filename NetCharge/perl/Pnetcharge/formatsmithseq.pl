#!/usr/bin/perl
#Read fasta sequences, format the age information and save it.
#Aug 30, 2012
#Hsiang-Yu Yuan

my $filename = $ARGV[0];
my $out_filename = "hm_h3n2_smith_seq.csv";

# specify the filename if user didn't give it
if ($filename eq '') {
 $filename = "seq/smith/smith_list.txt";
}

open FILE, "$filename" or die $!;

print "reading file: $filename\n";

my %labelhash = ();
my @tokens;
my @results;
my @removed_set;
my $sequence;

my $i=0;
while (my $line = <FILE>) { 
     
       if ($line =~ /^$/) {
         next;
       }

       @tokens = split(/\s/, $line);
       #print "$tokens[0] $tokens[9]\n";
       
       #print "$tokens[0]\n";
       $tokens[0] =~ m/(.*)/; #gb access 
       my $antigenic_cluster = $1;
       #print "$antigenic_cluster\n";      

       if ($tokens[1] =~ m/.*\/.*\/(\d\d)/) {
         $strain_name = $tokens[1];
         $isotime = $1;
       }
 
       #if ($isotime =~ m/(\d\d\d\d)\/(\d\d)\/(\d\d)/) {
       #    print "$1\n";
       #    print "$2\n";
       #    print "$3\n";
       #}

       $results[$i][0] = $antigenic_cluster;
       $results[$i][1] = $strain_name;
       $results[$i][2] = $isotime;       
       $i++;     
}
   

    
print "write to file: seq/$out_filename";
open OUTFILE, ">seq/$out_filename";   
print OUTFILE "CLUSTER,STRAIN,ISODATE\n";
for (my $j = 0; $j < scalar(@results); $j++) { 
       my $out_str = "$results[$j][0],$results[$j][1],$results[$j][2]";
       print OUTFILE "$out_str\n";
}
close OUTFILE;

 
