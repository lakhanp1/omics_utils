#!/usr/bin/perl
use strict;
use warnings;

open IN1,"trimmed_Tn-13_R1.fastq";
open OUT1,"> TA_trimmed_Tn-13_R1.txt";
open OUT2,"> num_trimmed_Tn-13_R1.txt";
my $total_num=0;
my $ta_num=0;
while(<IN1>){
       my $line1;
       if(/^@/){
           $line1=$_;
           $total_num++;
       }
       my $line2=<IN1>;
       my $line3=<IN1>;
       my $line4=<IN1>;
       if($line2=~/^TA/){
           $ta_num++;
           print OUT1 $line1,$line2,$line3,$line4;
       }
}
print OUT2 "sequence number:",$total_num,"\n","TA sequence number:",$ta_num,"\n";

