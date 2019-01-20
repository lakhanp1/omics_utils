use strict;
use warnings;

use Getopt::Long;

my $usage = 'This script splits the each read into overlapping reads of given LENGTH starting from 5\' end. 
If you are not getting any reads in output files, the read ID pattern is different than standard Illumina pattern. Please change the regex in the IF statement inside WHILE loop.
USAGE: perl fastqReadSpliterOverlappingSE.pl --fq <INPUT FASTQ FILE> --length <LENGTH OF NEW READS> --overlap <OVERLAP BETWEEN TWO CONSECUTIVE READS>';

my %options;
GetOptions(\%options, 'fq=s', 'length=i', 'overlap=i') or die("Error in command line arguments\n");


if(!defined $options{'fq'}){
	print  "ERROR: Please provide input fastq file\n";
	die $usage,"\n";
}
if(!defined $options{'length'}){
	print "Please provide length into which each read to split from 5' and 3' end\n";
	die $usage,"\n";
}
if(!defined $options{'overlap'}){
	print "Please provide required overlapp between consecutive reads\n";
	die $usage,"\n";
}


open(my $fh, $options{'fq'}) || die "cannot open $options{'fq'}\n";

$options{'fq'}=~m/(.*)\.fastq.*/;

my $r1 = $1.'_'.$options{'length'}.'_overlapping.fastq.gz';

open(my $oR1, "|gzip >$r1") or die "Cannot create file $r1: $!";

my $i = 1;
while(1){
	if(eof($fh)){
		last;
	}
	
	my $line = <$fh>;
	
	if($line=~m/^(@[\w-]+:\d+:[\w-]+(:\d+){4})/){
		#read ID
		# print $line;
		
		my $id = $1.':'.$i;
		
		#read seq
		my $seq = <$fh>;
		chomp $seq;
				
		#comment line
		my $comment = <$fh>;
		chomp $comment;
		
		#qual scores
		my $qual = <$fh>;
		chomp $qual;
		
		my $start = 0;
		while(($start+$options{'length'}) <= length $seq){
		
			print $oR1 $id.':'.$start,"\n";
			print $oR1 substr($seq, $start, $options{'length'}),"\n";
			print $oR1 $comment,"\n";
			print $oR1 substr($qual, $start, $options{'length'}),"\n";
			
			$start += $options{'overlap'};
		}
		
		if(($start+$options{'length'}) < (length $seq)-($options{'length'}/2)){
			print $oR1 $id.':',(length($seq) - $options{'length'}),"\n";
			print $oR1 substr($seq, length($seq) - $options{'length'}),"\n";
			print $oR1 $comment,"\n";
			print $oR1 substr($qual, length($seq) - $options{'length'}),"\n";
		}
	}
	
	$i++;
}

close($oR1);

