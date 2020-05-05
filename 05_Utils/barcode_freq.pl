use strict;
use warnings;

my $usage = 'This script extract the most 10 frequent 5\' sequences counts from a fastq file.
USAGE: zcat barcode_file.txt | perl barcode_freq.pl
printf "sample\tATCGCA | perl barcode_freq.pl"
';

while(my $line = <STDIN>){
	chomp $line;
	my ($sp, $bc) = split(/\t/, $line);
	my $bcLen = length($sp =~s/^.*_($bc[ATGC]*)_.*/$1/r) + 1;

	# print "$sp\t$bc\t$bcLen\n";
	my $i=3;
	my %freq=();
	
	open(my $fastq, "zcat $sp*.fastq.gz | ") or die;
	
	while(<$fastq>){
		
		# print $_;
		# last;

		if($i % 4 == 0){
			$freq{substr($_, 0, $bcLen)}++;
		}
		
		$i++;
		
	}
	
	my $count = 0;
	foreach my $ad(sort{$freq{$b} <=> $freq{$a}}keys %freq){
		print "$bc\t$ad\t$freq{$ad}\n";
		$count++;
		if($count >= 10){
			last;
		}
	}
	
	# print "\n";
}
