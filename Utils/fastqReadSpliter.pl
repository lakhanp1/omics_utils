use strict;
use warnings;

use Getopt::Long;

my $usage = 'This script splits the each read into 2 mates: forward read from 5\' end and reverse read from 3\' end. The R2 is reversed to make the orientation of mates -->	<--. You can comment the reverse statement inside script if you want --->	---> mates.
If you are not getting any reads in output files, the read ID pattern is different than standard Illumina pattern. Please change the regex in the IF statement inside WHILE loop.
USAGE: perl fastqReadSpliter.pl --fq <INPUT FASTQ FILE> --length <LENGTH OF NEW READS>';

my %options;
GetOptions(\%options, 'fq=s', 'length=i') or die("Error in command line arguments\n");

if(!defined $options{'fq'}){
	print  "ERROR: Please provide input fastq file\n";
	die $usage,"\n";

}
if(!defined $options{'length'}){
	print "Please provide length into which each read to split from 5' and 3' end\n";
	die $usage,"\n";
}



open(my $fh, $options{'fq'}) || die "cannot open $options{'fq'}\n";

$options{'fq'}=~m/(.*)\.fastq.*/;
my $r1 = $1.'_'.$options{'length'}.'_R1.fastq.gz';
my $r2 = $1.'_'.$options{'length'}.'_R2.fastq.gz';

open(my $oR1, "|gzip >$r1") or die "Cannot create file $r1: $!";
open(my $oR2, "|gzip >$r2") or die "Cannot create file $r2: $!";

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
		print $oR1 $id."/1\n";
		print $oR2 $id."/2\n";
		
		#read seq
		my $seq = <$fh>;
		chomp $seq;
		
		print $oR1 substr($seq, 0, $options{'length'}),"\n";
		# print $oR2 substr(reverse($seq), 0, $options{'length'}),"\n";
		print $oR2 substr($seq, length($seq) - $options{'length'}),"\n";
		
		#3rd line
		$seq = <$fh>;
		chomp $seq;
		print $oR1 $seq,"\n";
		print $oR2 $seq,"\n";
		
		#qual scores
		$seq = <$fh>;
		chomp $seq;
		
		print $oR1 substr($seq, 0, $options{'length'}),"\n";
		# print $oR2 substr(reverse($seq), 0, $options{'length'}),"\n";
		print $oR2 substr($seq, length($seq) - $options{'length'}),"\n";
	}
	
	$i++;
}

close($oR1);
close($oR2);

