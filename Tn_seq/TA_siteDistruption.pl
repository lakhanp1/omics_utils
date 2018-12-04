#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;
use Math::Round qw(nearest);

my $usage = "This script calculates the frequency table for percentage of TA sites disrupted plot
USAGE: perl TA_siteDistruption.pl <TA_counts.bed file 1>
";


if($#ARGV == -1){
	print STDERR "Please provide atleast one bedtools coverage output file generated for the TA site location BED file\n";
	die $usage;
}

#empty %geneDistruptionFreq to store freq data for each sample
my %geneDistruptionFreq = ();
my %TA_insertionStats = ();

foreach(map{$_*5}(0..20)){
	$geneDistruptionFreq{$_} = 0;
}

my $i = 0;
my %TA_count = ();
my @bedFields = qw(chr start end name count);


my $file = shift @ARGV;
my $sampleName = &fileparse($file);

my $TA_statsFile = $sampleName."_disruptionFreq.stats";
my $genewiseTA_statsFile = $sampleName."_genewiseIns.stats";

open(my $o1, '>', $TA_statsFile) or die "Cannot open file $TA_statsFile: $!";
open(my $o2, '>', $genewiseTA_statsFile) or die "Cannot open file $genewiseTA_statsFile: $!";

open(my $fh, $file) or die "Cannot open file $file: $!";

my ($totalGenomicTA, $totalGenomicInsertions) = (0, 0);
my ($totalGeneTA, $totalGeneInsertion) = (0, 0);
my ($totalGenes, $totalGenesDisrupted) = (0, 0);
my $geneTA_count = 0;
my $gene = '';
my $ins = 0;

print $o2 join("\t", qw(Gene total_TA_sites TA_sites_disrupted Fraction)),"\n";

while(<$fh>){
	@TA_count{@bedFields} = split(/\t/, $_);
	
	#TA sites outside gene
	if($TA_count{'name'}=~m/_(up|down)_\d+$/){
		$totalGenomicTA++;
		if($TA_count{'count'} > 0){
			$totalGenomicInsertions++;
		}
	}
	#For TA sites withing gene		
	elsif($TA_count{'name'}=~m/TA_\d+_(.*)_0$/){		
		$totalGenomicTA++;
		$totalGeneTA++;
		
		#calculate the fraction of TA sites distrupted for a gene
		if($gene ne $1){
			if($gene ne ''){				
				#Calculate the fraction and print to stats file
				my $fraction = sprintf("%.2f", $ins*100/$geneTA_count);
				print $o2 join("\t", $gene, $geneTA_count, $ins, $fraction),"\n";
				
				#nearest multiple of 5
				$fraction = &nearest(5, $fraction);
				$geneDistruptionFreq{$fraction}++;
				
			}
			
			$gene = $1;
			$ins = 0;
			$geneTA_count = 0;
		}
		
		$geneTA_count++;
		if($TA_count{'count'} > 0){
			$ins++;
			$totalGenomicInsertions++;
			$totalGeneInsertion++;
		}		
	}		
}

#Run for the last TA site in the file
if($gene ne ''){
	my $fraction = sprintf("%.2f", $ins*100/$geneTA_count);
	print $o2 join("\t", $gene, $geneTA_count, $ins, $fraction),"\n";
	$fraction = &nearest(5, $fraction);
	$geneDistruptionFreq{$fraction}++;
}

$TA_insertionStats{'totalGenomicTA'} = $totalGenomicTA;
$TA_insertionStats{'totalGenomicInsertions'} = $totalGenomicInsertions;
$TA_insertionStats{'totalGeneTA'} = $totalGeneTA;
$TA_insertionStats{'totalGeneInsertion'} = $totalGeneInsertion;




print $o1 '%disruption',"\t",$sampleName,"\n";

foreach(sort {$a <=> $b} keys %geneDistruptionFreq){
	print $o1 join("\t", $_, $geneDistruptionFreq{$_}),"\n";
}

foreach(sort keys %TA_insertionStats){
	print $o1 join("\t", $_, $TA_insertionStats{$_}),"\n";
}


close($fh);
close($o1);
close($o2);


