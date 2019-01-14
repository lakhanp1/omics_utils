#!/usr/bin/perl -w

use strict;
use warnings;

my $rawFile = "rawData.txt";

my $heatmap = "heatmap.txt";
my $textFile = "heatmapText.txt";

open(my $fh, $rawFile) or die "Cannot open file $rawFile: $!";
open(my $mapOut, ">", $heatmap) or die "Cannot create $heatmap: $!";
open(my $textOut, ">", $textFile) or die "Cannot create $textFile: $!";

my $header = 0;
my @hmIds = ();

while(my $line = <$fh>){
	chomp $line;
	if(!$header){
		#header line
		@hmIds = map{/^(\w)(\w)(\d)/? [$1, $2, $3, $_] : [$_]}split(/\t/, $line);
		$header = 1;
		#map{print "@{$_}\n";}@hmIds;
		next;
	}
	
	my $id = undef;
	my $i =0;
	my $end = 1600;
	my $start = 0;
	foreach(split(/\t/, $line)){
		if($i == 0){
			$id = $_;
		}
		elsif($i>1){
			if($hmIds[$i]->[2] == 1){
				$start = 0;
				$end = 1999;
			}
			print $mapOut $id," ",$start," ",$end," ",$_," id=",$i-2,"\n";
			print $textOut $id," ",$start," ",$end," ",$hmIds[$i]->[1].$hmIds[$i]->[2],"\n";
			$start = $end + 1;
			$end+=2000;
		}
		
		$i++;
	}
	
}