#!/usr/bin/perl -w

use strict;
use warnings;


my $usage = 'This script parse the \'bedtools genomecov -bga -split \' output and returns only those regions which has read.
USAGE: perl getCoveredRegionsAsBed.pl <bedtools output>
';

if(@ARGV != 1){
	print "ERROR: Please provide the bedtools coverage output file\n";
	die $usage;
}

my $newBed = 0;
my ($chr, $start, $end) = ('', undef, undef);

open(my $fh, $ARGV[0]) or die "Cannot open file $ARGV[0]\n";

while(<$fh>){
	chomp $_;
	my @tmp = split(/\t/, $_);
	

	if($tmp[3] != 0 && $newBed == 0){
		$start = $tmp[1];
		$chr = $tmp[0];
		$end = undef;
		$newBed = 1;
	}
	elsif($tmp[3] == 0 && $newBed == 1 && $chr eq $tmp[0]){
		print join("\t", $chr, $start, $end),"\n";
		($chr, $start, $end) = ('', undef, undef);
		$newBed = 0;
	}
	elsif($chr ne $tmp[0] && $newBed == 1){
		print join("\t", $chr, $start, $end),"\n";
		($chr, $start, $end) = ('', undef, undef);
		$newBed = 0;
	}
	
	$end = $tmp[2];
}

close($fh);

