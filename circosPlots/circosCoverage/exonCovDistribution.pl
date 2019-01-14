use strict;
use warnings;

use Statistics::Basic qw(median mean);
use Getopt::Long;

#This script parse the output of command: samtools bedcov <BED> <BAM> | perl -F"\t" -ane 'chomp @F; print join("\t", @F,sprintf("%.3f", $F[-1]/($F[2]-$F[1]))),"\n"' > <NAME>_exonCov.bed
#Calculate the %of exons with coverage greater than x coverage value.
# If want to skip a specific chromosome (e.g. chr Y in case of female sample), provide the --skipChr argument.

my %options;
GetOptions(\%options, 'skipChr=s') or die("Error in command line arguments\n");

my $skip = 0;
if(defined $options{'skipChr'} ){
	chomp $options{'skipChr'};
	$skip = 1;
}


my %freq=();
my $sum=0;
my $i=0;

my @bin = (0, 1, 3, 5, map{10*$_}1..20);
my @vals = ();
my $noCov = 0;

while(<>){
	next if($skip && /^($options{'skipChr'})\t/);
	
	my @F = split(/\t/, $_);
	push(@vals, $F[4]);

	$i++;
	$sum+=$F[4];
	if($F[4]==0){
		$noCov++;
	}
	else{
		foreach my $j(@bin){
			$F[4] > $j ? $freq{$j}++ : last;
		}
	}
}

my $median = median(@vals);

foreach(sort {$a <=> $b} keys %freq){
	print sprintf("%.2f", $freq{$_}*100/$i),"\t";
}
print sprintf("%.2f", $sum/$i),"\t",$median;

print "\n";
