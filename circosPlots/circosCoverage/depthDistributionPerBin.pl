#!/usr/bin/perl -w

use strict;
use warnings;


# This script takes input stream from 'bedtools intersect -loj -a <regions.bed> -b <meanCoverageOut.bed>'
# Process each line and calculate the fraction of exons in each bin with coverage = 0/1-5/6-10/11<


my %freq;

while(<>){
	chomp;
	my @tmp=split(/\t/, $_);
	my $key=$tmp[0].":".$tmp[1]."-".$tmp[2];
	
	if(!exists $freq{$key}){
		$freq{$key}->{'0'}=0;
		$freq{$key}->{'1-5'}=0;
		$freq{$key}->{'5-10'}=0;
		$freq{$key}->{'11'}=0;
		$freq{$key}->{'sum'}=0;
	}
	
	$freq{$key}->{'sum'}++;
	
	if($tmp[-1] < 1){
		$freq{$key}->{'0'}++;
	}elsif($tmp[-1] <= 5){
		$freq{$key}->{'1-5'}++;
	}elsif($tmp[-1] <= 10){
		$freq{$key}->{'5-10'}++;
	}else{
		$freq{$key}->{'11'}++;
	}
}

foreach(sort{my @aa=$a=~m/(.*):(\d+)-(\d+)/; my @bb=$b=~m/(.*):(\d+)-(\d+)/; $aa[0] cmp $bb[0] || $aa[1] <=> $bb[1]}keys %freq){
	#print "$_\t$freq{$_}->{'sum'}\n";
	print join(" ", split(/[:-]/, $_))," ";
	print sprintf("%.2f", $freq{$_}->{'0'}/$freq{$_}->{'sum'}).",";
	print sprintf("%.2f", $freq{$_}->{'1-5'}/$freq{$_}->{'sum'}).",";
	print sprintf("%.2f", $freq{$_}->{'5-10'}/$freq{$_}->{'sum'}).",";
	print sprintf("%.2f", $freq{$_}->{'11'}/$freq{$_}->{'sum'}),"\n";
}

