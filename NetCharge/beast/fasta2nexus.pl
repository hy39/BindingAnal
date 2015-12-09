#!/usr/bin/perl
#change fasta to phylip
#Date:May 30, 2009
#Author:Hsiang-yu Yuan


#Set a Variable

# Run a while loop

#First: read a list of filenames
$count = 1;
$n1=$ARGV[0];


@files = <$n1/*.fa>;

foreach (@files) {
    $file = $_;

#    $str = "B_1041_aligned_subsetyears_30_" . "$index" . ".fa";
#    print "$str";
    $n1=$ARGV[0];
    $n2=$ARGV[1];
#    print "$n1 $n2";

    $infile = "$file";
    #print $infile;
    $outfile = "$file" . '.nexus';
    #print "$outfile\n";
    #print("java -cp readseq.jar run format=Phylip output=$outfile  $infile\n");
    system("java -cp readseq.jar run format=nexus inform=fasta output=$outfile  $infile");
    $count++;
}
