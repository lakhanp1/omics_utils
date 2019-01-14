#!/usr/bin/perl -w


#This script returns the coordinates on whole karyotype such that the points are placed equidistant on the whole circle.

use strict;
use warnings;

#This is karyotype hash
my %data = (1 => 249250621,
2 => 243199373,
3 => 198022430,
4 => 191154276,
5 => 180915260,
6 => 171115067,
7 => 159138663,
8 => 146364022,
9 => 141213431,
10 => 135534747,
11 => 135006516,
12 => 133851895,
13 => 115169878,
14 => 107349540,
15 => 102531392,
16 => 90354753,
17 => 81195210,
18 => 78077248,
19 => 59128983,
20 => 63025520,
21 => 48129895,
22 => 51304566,
'X' => 155270560,
'Y' => 59373566);

my @chr = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);


my $dist = 119000000;		#Change this parameter as per your need

my $remaining = 0;
my @newPos = ();

for(my $i=0; $i<@chr; $i++){
	my $total = $remaining + $data{$chr[$i]};
	
	if($total < $dist){
		#print "$remaining < $dist: NEXT\n";
		$remaining = $total;
		next;
	}
	
	my $newPt = 0;
	
	while(($total - $dist) > 0){
		# print "*CHR $chr[$i]\t$remaining*\n";
		$newPt = $newPt + $dist - $remaining;
		
		$remaining = 0;
		$total = $data{$chr[$i]} - $newPt;
		print join("\t", $chr[$i], $newPt),"\n";
	}
	
	$remaining = $total;
}