#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = '
This scripts downloads the SRA files. It accepts a file having SRA IDs on each line and downloads the SRA from NCBI ftp directory.
Usage: perl sraDownload.pl <Comma_Separated_SRA_IDs>';


if($#ARGV != 0){
	print "No valid SRA ids given.\n";
	die $usage;
}

my @sras = split(/,/, $ARGV[0]);

foreach my $sra(@sras){
	print "Downloading SRA $sra from: ";
	my $link = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/";
	$sra=~m/((([SED]RR)\d{3})\d+)/;
	$link.=$3."/".$2."/".$1."/$1.sra";
	print $link,"\n";
	
	system("wget $link");
	if($? != 0){
		print "Status: SRA $sra was not downloaded\n";
	}
	else{
		print "Status: SRA $sra download completed successfully\n";
	}
}