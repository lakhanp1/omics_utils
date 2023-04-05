#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

#bedtools coverage gives output for a given bam file/STDIN of bam file. (eg. 'samtools view -b -F 4 <bamFile> NC_000077.6:101488264-101553955'). But it gives flanking regions also. This script extract the coverage for the specific region from bedtools coverage output.

my %options;
GetOptions(\%options, 'region=s') or die("Error in command line arguments\n");

if(!defined $options{'region'}){
	die "Please provide chromosomal region in format CHR_NAME:START-END\n";
}

my ($chr, $start, $end) = (undef, undef, undef);

if($options{'region'}=~m/(\S+):(\d+)-(\d+)/){
	($chr, $start, $end) = ($1, $2, $3);
}
else{
	die "Genome region format is wrong. Please provide chromosomal region in format CHR_NAME:START-END\n";
}

my $cmp = $start;
my $inRegion = 0;

while(my $i=<>){
	chomp $i;
	my @tmp = split(/\s+/, $i);
	
	#only the required chromosomes
	if($i=~m/^$chr/){
		
		if($tmp[2] > $start && $inRegion==0){
			$inRegion = 1;
			print join("\t", $tmp[0],$start,$tmp[2],$tmp[3]),"\n";
			next;
		}
		
		if($tmp[2] >= $end){
			$inRegion = 0;
			print join("\t", $tmp[0],$tmp[1],$end,$tmp[3]),"\n";
			last;
		}
		
		if($inRegion==1){
			print $i,"\n";
		}
	}
	
}