#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

#This script gets the 4th column value of 'samtools mpileup', start and end region. This value is used to calculate median and average coverage.
#NC_000069.6     41203795        N       34      aaaaAAAaaaAaAAAaAAAaaaAAAaaAAAaa^]a^]a  BBF<F<FFFFFFFFFFFFkFFFFkF/FFFBFFFF
#NC_000069.6     41203430        N       3       ttT     BFF
#NC_000069.6     41203438        N       1       T       F

#input coming for above samtools mpileup output will be: 34 \n 3 \n 

my %options;
GetOptions(\%options, 'start=i', 'end=i') or die("Error in command line arguments\n");

if(!defined $options{'start'} || !defined $options{'end'}){
	die "Please provide start and end: $!";
}


my %covFreq=();

my %stats = (
	'medianCov' => 0,
	'averageCov' => 0,
	'length' => 0,
	'sum' => 0,
	'bases' => 0);


while(my $i=<>){
	chomp $i;
	$stats{'bases'}++;
	$covFreq{$i}++;
	$stats{'sum'}+=$i;
}

$stats{'length'} = $options{'end'}-$options{'start'};
$covFreq{0} = $stats{'length'}-$stats{'bases'};

my $midPoint = 0;

if($stats{'length'}%2 == 0){
	$midPoint = $stats{'length'}/2
}
else{
	$midPoint = ($stats{'length'}+1)/2
}

#get the median coverage
my $temp = 0;
foreach(sort{$a <=> $b} keys %covFreq){
	$temp+=$covFreq{$_};
	# print "$_\t$covFreq{$_}\ttemp=$temp\tmid=$midPoint\n";
	if($temp>=$midPoint){
		$stats{'medianCov'} = $_;
		last;
	}
}

$stats{'averageCov'} = sprintf("%.3f", $stats{'sum'}/$stats{'length'});

print "Non zero coverage bases = $stats{'bases'}\n";
print "MedianCoverage = $stats{'medianCov'}\n";
print "AverageCoverage = $stats{'averageCov'}\n";
