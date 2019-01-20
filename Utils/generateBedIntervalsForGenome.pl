#!/usr/bin/perl -w

use strict;
use warnings;

#This script generates the BED file of intervals with given interval length for all the chromosomes. Input is a tab separated file with Chromosome name in first column and its size in second column

my $usage = "#This script generates the BED file of intervals with given interval length for all the chromosomes. Input is a tab separated file with Chromosome name in first column and its size in second column
USAGE: generateBedIntervalsForGenome.pl <genome.size> <interval length>
";

if($#ARGV != 1){
	print STDERR "Please provide chromosome size file and bed interval length\n";
	die $usage;
}

my $size = 1000;

if($ARGV[1] =~m/^\d+$/){
	$size = $ARGV[1];
}
else{
	print STDERR "Please provide an integer as bed interval length\n";
	die $usage;
}


open(my $fh, $ARGV[0]) or die "Cannot open file $ARGV[0]: $!";
while(<$fh>){
	chomp;
	my ($chr, $len) = split(/\t/, $_);
	
	my $regionsFile = $chr."_regions.bed";
	open(my $out, '>',$regionsFile) or die "Cannot create $regionsFile: $!";
	
	for(my $i=1; $i<=$len; $i+=$size){
		if($i+$size-1 > $len){
			print $out join("\t", $chr, $i, $len, "$chr:$i-$len"),"\n";
		}
		else{
			print $out join("\t", $chr, $i, $i+$size-1, join("", "$chr:$i-",$i+$size-1)),"\n";
		}
	}
	
	close($out);
	
}



close($fh);