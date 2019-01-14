#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use DataFilePrepare qw(karyotypeData histogramAndLinks setParameters getQueryLength %circosFiles setChrUnit);

my %options;

GetOptions(\%options, 'reference=s', 'bam=s', 'coverageScript=s', 'skipChr=s', 'max=i', 'help|h') or die("Error in command line arguments\n");

#perl -I /home/lakhanp/scripts/ /home/lakhanp/scripts/genomeCovCircos.pl --reference /home/lakhanp/database/Mouse/mm_m38/reference/GCF_000001635.24_GRCm38.p4_genomicChromosomes.fna --bam /home/lakhanp/Analysis/CoreProjects/5_MiaoKaiSC/mapping/SC11/SC11_recalibrated.bam --max 30  --coverageScript /home/lakhanp/scripts/samtoolsAvgCoverage.pl


if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
}


if(!$options{'reference'}){
	print STDERR "Error: Please provide the reference genome file in FASTA format\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'coverageScript'}){
	print STDERR "Error: Please provide the coverageScript file path\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'bam'}){
	print STDERR "Error: Please provide the alignment file in BAM format\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(defined $options{'max'}){
	print "NOTE: Will lump all positions in the genome having feature coverage greater than or equal to -max into the -max histogram bin.\n";
}


foreach(keys %options){
	print "$_\t\t$options{$_}\n";
}

print "\n";

my $dataStorePath = "circos/data";
my $confStorePath = "circos/etc";
my $bedFile = $dataStorePath.'/chrBins.bed';
my $avgCov = $dataStorePath.'/averageCov.txt';
my $medianCov = $dataStorePath.'/medianCov.txt';
my $uncoveredFraction = $dataStorePath.'/uncoveredFraction.txt';
my $readCount = $dataStorePath.'/readCount.txt';
my $textFile = $dataStorePath.'/text.txt';

my $histInterval = 10000;			#this is determined internally based on the chromosome length

#circos runtime parameters
my %config;
#1) The chromosomes_unit value is used as a unit (suffix "u") to shorten values in other parts of the configuration file.
$config{chromosomes_units} = 100000;


#create karyotype file
my $idMap = &karyotypeData($options{'reference'}, $dataStorePath);

if(defined $options{'skipChr'}){
	foreach(split(/,/, $options{'skipChr'})){
		print "Not plotting chromosome: $_\n";
		delete $idMap->{$_};
	}
}

#decide the 'chromosomes_unit', 'multiplier' and 'histInterval' parameters based on genome size so that the ticks and tick lables will be appropriate
($config{chromosomes_units}, $config{multiplier}) = &setChrUnit();
$histInterval = $config{chromosomes_units} * 2;

print "Histogram interval = $histInterval\n";
print "Chromosome units = $config{chromosomes_units}\n";
print "Tick lable multiplier = $config{multiplier}\n";

my ($maxAvgCov, $maxMedCov, $maxCount) = (0,0,0);

print "Creating the bins for plotting histograms...\n";
open(my $bOut, ">", $bedFile) or die "Cant create file ",$bedFile,": $!";
open(my $avgOut, ">", $avgCov) or die "Cant create file ",$avgCov,": $!";
open(my $medOut, ">", $medianCov) or die "Cant create file ",$medianCov,": $!";
open(my $uncovOut, ">", $uncoveredFraction) or die "Cant create file ",$uncoveredFraction,": $!";
open(my $cntOut, ">", $readCount) or die "Cant create file ",$readCount,": $!";
open(my $txtOut, ">", $textFile) or die "Cant create file ",$textFile,": $!";

foreach(sort keys %{$idMap}){
	for(my $i = 0; $i < $idMap->{$_}->[1]; $i += $histInterval){
		#print $_," => $idMap->{$_}->[1]\t",$i,"\t",$i+$histInterval-1,"\t",0,"\n";
		if($i+$histInterval-1 > $idMap->{$_}->[1]){
			
			print $bOut join("\t", $_, $i, $idMap->{$_}->[1]),"\n";
			my ($count, $uncoveredFraction, $med, $avg) = &samtoolsAvgCoverage($_, $i, $idMap->{$_}->[1]);
			
			#Print the max median value if median is more than defined limit
			if(defined $options{'max'} && $med > $options{'max'}){
				print $medOut join("\t", $_, $i, $idMap->{$_}->[1], $options{'max'}),"\n";
			}
			else{
				print $medOut join("\t", $_, $i, $idMap->{$_}->[1], $med),"\n";
			}
			
			#Print the max avg value if average is more than defined limit
			if(defined $options{'max'} && $avg > $options{'max'}){
				print $avgOut join("\t", $_, $i, $idMap->{$_}->[1], $options{'max'}),"\n";
			}
			else{
				print $avgOut join("\t", $_, $i, $idMap->{$_}->[1], $avg),"\n";
			}
			
			
			
			print $uncovOut join("\t", $_, $i, $idMap->{$_}->[1], $uncoveredFraction),"\n";
			print $cntOut join("\t", $_, $i, $idMap->{$_}->[1], $count),"\n";
		}
		else{
			
			print $bOut join("\t", $_, $i, $i+$histInterval-1),"\n";
			my ($count, $uncoveredFraction, $med, $avg) = &samtoolsAvgCoverage($_, $i, $i+$histInterval-1);
			
			if(defined $options{'max'} && $med > $options{'max'}){
				print $medOut join("\t", $_, $i, $i+$histInterval-1, $options{'max'}),"\n";
				print $avgOut join("\t", $_, $i, $i+$histInterval-1, $options{'max'}),"\n";
			}
			else{
				print $medOut join("\t", $_, $i, $i+$histInterval-1, $med),"\n";
				print $avgOut join("\t", $_, $i, $i+$histInterval-1, $avg),"\n";
			}
			
			
			print $uncovOut join("\t", $_, $i, $i+$histInterval-1, $uncoveredFraction),"\n";
			print $cntOut join("\t", $_, $i, $i+$histInterval-1, $count),"\n";
		}		
	}
}

print $txtOut join("  ",'NC_000073.6', '0', '145441459', "Max_Fraction_Covered=1", 'rad0=0.40r,rad1=0.53r,color=vdorange'),"\n";
print $txtOut join("  ",'NC_000069.6', '0', '160039680', "Max_Mean_coverage=$maxAvgCov", 'rad0=0.70r,rad1=0.83r,color=vdgreen'),"\n";
print $txtOut join("  ",'NC_000067.6', '0', '195471971', "Max_Median_coverage=$maxMedCov", 'rad0=0.40r,rad1=0.53r,color=vdorange'),"\n";
print $txtOut join("  ",'NC_000071.6', '0', '151834684', "Max_Read_Freq=$maxCount", 'rad0=0.55r,rad1=0.68r,color=vdblue'),"\n";

print "Max average coverage = $maxAvgCov\n";
print "Max median coverage = $maxMedCov\n";
print "Max read count = $maxCount\n";
print "Done...\n\n";




close($bOut);
close($avgOut);
close($medOut);
close($uncovOut);
close($cntOut);
close($txtOut);














sub samtoolsAvgCoverage{
	# samtools mpileup -r NC_000087.7:10000000-11999999 --ff 4 SC11bt2Mapped.bam | awk '{ print $4 }' | perl /home/lakhanp/scripts/samtoolsAvgCoverage.pl --start 10000000 --end 11999999
	#NC_000069.6     41203795        N       34      aaaaAAAaaaAaAAAaAAAaaaAAAaaAAAaa^]a^]a  BBF<F<FFFFFFFFFFFFkFFFFkF/FFFBFFFF
	#NC_000069.6     41203430        N       3       ttT     BFF
	#NC_000069.6     41203438        N       1       T       F
	
	my $command = "samtools mpileup -r $_[0]:$_[1]-$_[2] $options{'bam'} | awk '{ print \$4 }' | perl $options{'coverageScript'} --start $_[1] --end $_[2]
	echo \"PIPESTATUS:\${PIPESTATUS[@]}\"";
	
	my $res = `$command`;
	
	if($res !~ m/PIPESTATUS:0 0 0/){
		print $res;
		die "Could not calculate average coverage with command $command: $!";
	}
	
	my ($count, $uncoveredFraction, $med, $avg) = (undef, undef, undef, undef);
	
	#Number of reads in region
	$count=`samtools view -c -F 4 $options{'bam'} $_[0]:$_[1]-$_[2]`;
	
	chomp $count;
	
	if($res=~m/Non zero coverage bases = (\d+)/){
		$uncoveredFraction = sprintf("%.3f", ($_[2] - $_[1] + 1 - $1)*100/($_[2] - $_[1] + 1));
	}
	if($res=~m/MedianCoverage = (\d+)/){
		$med=$1;
	}
	if($res=~m/AverageCoverage = (\d+\.\d+)/){
		$avg=$1;
	}
	
	#set the max values
	if($med > $maxMedCov){
		$maxMedCov = $med;
	}
	if($avg > $maxAvgCov){
		$maxAvgCov = $avg;
	}
	if($count > $maxCount){
		$maxCount = $count;
	}
	
	
	if(defined $count && defined $uncoveredFraction && defined $med && defined $avg){
		return($count, $uncoveredFraction, $med, $avg);
	}
	else{
		die "Could not calculate average coverage with command $command: $!";
	}
	
}






__END__


=head1 NAME


=head1 SYNOPSIS

perl -I <Code dir Path> /complete/path/to/genomeCovCircos.pl --reference <reference file> --bam <bam file>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

Something

=head1 OPTIONS

=over 30

=item B<--help>

Show this scripts help information.

=item B<--reference>

[STR] Reference FASTA sequence

=item B<--coverageScript>

[STR] Complete path to samtoolsAvgCoverage.pl

=item B<--bam>

[STR] BAM file 

=item B<--max>

[INT] Will lump all positions in the genome having feature coverage greater than or equal to -max into the -max histogram bin.

=item B<--skipChr>

[STR] Chromosome to skip in coverage calculation

=back


=cut



