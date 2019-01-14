#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;


my %options;

GetOptions(\%options, 'reference=s', 'bam=s', 'coverageScript=s', 'bed=s', 'skipChr=s', 'help|h') or die("Error in command line arguments\n");

#perl -I /home/lakhanp/scripts/ /home/lakhanp/scripts/genomeCovCircos.pl --reference /home/lakhanp/database/C_albicans/SC5314_A21/reference/C_albicans_SC5314_A21_current_chromosomes.fasta --bam $PWD/Sample1.bam --coverageScript /home/lakhanp/scripts/samtoolsAvgCoverage.pl




if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 1, VERBOSE => 2});
	exit 1;
}


# if(!$options{'reference'}){
	# print STDERR "Error: Please provide the reference genome file in FASTA format\n";
	# &pod2usage({EXIT => 2, VERBOSE => 0});
# }

if(!$options{'coverageScript'}){
	print STDERR "Error: Please provide the coverageScript file path\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'bam'}){
	print STDERR "Error: Please provide the alignment file in BAM format\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'bed'}){
	print STDERR "Error: Please provide the regions in BED format\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}



my $dataStorePath = ".";
my $covFile = $dataStorePath.'/chrBins.bed';

my ($maxAvgCov, $maxMedCov, $maxCount) = (0,0,0);


open(my $fh, $options{'bed'}) or die "Cannot open file ",$options{'bed'},"$!\n";


while(<$fh>){
	chomp $_;
	if(/^\w\S*\s+\d+\s+\d+/){
		my @tmp = split(/\t/, $_);
		my ($count, $uncoveredFraction, $med, $avg) = &samtoolsAvgCoverage(@tmp[0..2]);
		
		print join("\t", $_, $tmp[2]-$tmp[1], $count, $uncoveredFraction, $avg, $med),"\n";
		
	}
}



close($fh);




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
		die "Could not calculate coverage with command $command: $!";
	}
	
}








__END__


=head1 NAME


=head1 SYNOPSIS

perl -I <Code dir Path> /complete/path/to/genomeCovCircos.pl --bam <bam file> --coverageScript <path/to/samtoolsAvgCoverage.pl> --bed <exon regions in BED format>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script parse the SORTED exon.bed file and calculates the number of mapped reads, uncovered fraction, median and average coverage for each region in BED file. 

=head1 OPTIONS

=over 30

=item B<--help>

Show this scripts help information.

=back


=cut



