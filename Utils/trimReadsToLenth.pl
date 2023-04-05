#!/usr/bin/perl -w

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Basename;


my %options;
my $isPaired = 1;

GetOptions(\%options, 'size=i', '1=s', '2=s', 'help|h') or die("Error in command line arguments\n");


if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
	#exit 1;
}


if(!$options{'size'}){
	print STDERR "Error: Please provide the barcode file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


if(!$options{'1'}){
	print STDERR "Error: Please provide the Read_1 file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


if(!$options{'2'}){
	print  "Warning: Read 2 file not provided. Assuming single end data\n";
	$isPaired = 0;
}


my $isGz = 0;
if($options{1}=~m/gz$/){
	$isGz=1;
}


my $r1Trimmed = basename($options{1});
$r1Trimmed = 'trimmed_'.$r1Trimmed;
open(my $out1, "|gzip >$r1Trimmed") or die "Cannot create file $r1Trimmed: $!";

#Read fastq data
open(my $fh1, $isGz ? "gzip -dc $options{1} |" : $options{1}) or die "Cannot open file $options{1}: $!";
my $fh2 = undef;

my $r2Trimmed = '';
my $out2;

if($isPaired){
	$r2Trimmed = basename($options{2});
	$r2Trimmed = 'trimmed_'.$r2Trimmed;
	open($out2, "|gzip >$r2Trimmed") or die "Cannot create file $r2Trimmed: $!";
	open($fh2, $isGz ? "gzip -dc $options{2} |" : $options{2}) or die "Cannot open file $options{2}: $!";
}


my ($p1, $p2, $outBar);
my ($trimmedReads, $filteredReads) = (0,0);
my $absLen = abs($options{'size'});

if($isPaired){
	#for paired end data
	while(1){
		if(eof($fh1)){
			last;
		}
		
		#read line1: headers
		$p1 = <$fh1>;
		$p2 = <$fh2>;
		
		if($p1=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N)/ && $p2=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N)/){
			#@SIM:1:FCX:1:15:6329:1045 1:N:0:2
			
			#read line 2: sequence
			my $sq1 = <$fh1>;
			my $sq2 = <$fh2>;
			if(length($sq1) > $absLen && length($sq2) > $absLen){
				$trimmedReads++;
				
				$p1 .= substr($sq1, -1*$absLen);
				$p2 .= substr($sq2, -1*$absLen);
				
				#read line 3: +
				$p1 .= <$fh1>;
				$p2 .= <$fh2>;
				
				#read line 4: qual
				#trim the qual of 5_prime barcode
				$p1 .= substr(<$fh1>, -1*$absLen);
				$p2 .= substr(<$fh2>, -1*$absLen);
				
				#print the reads
				print $out1 $p1;
				print $out2 $p2;
			}
			else{
				$filteredReads++;
			}
	
		}
	}
}
else{
	#for single end
	while(1){
		if(eof($fh1)){
			last;
		}
		
		#read line1: headers
		$p1 = <$fh1>;
		
		if($p1=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N)/){
			#@SIM:1:FCX:1:15:6329:1045 1:N:0:2
			
			#read line 2: sequence
			my $sq1 = <$fh1>;
			if(length($sq1) > $absLen){
				$trimmedReads++;
				
				$p1 .= substr($sq1, -1*$absLen);
				
				#read line 3: +
				$p1 .= <$fh1>;
				
				#read line 4: qual
				#trim the qual of 5_prime barcode
				$p1 .= substr(<$fh1>, -1*$absLen);
				
				#print the reads
				print $out1 $p1;
			}
			else{
				$filteredReads++;
			}			
		}
	}
}




print "Total reads: ",$trimmedReads+$filteredReads,"\n";
print "Trimmed reads: ",$trimmedReads,"\n";
print "Reads cannot be trimmed because of length smaller than $absLen: ",$filteredReads,"\n";











__END__


=head1 NAME


=head1 SYNOPSIS

perl trimReadsToLenth.pl --size <Trim length> -1 <Mate1> -2 <Mate2>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script trims the reads to a given size from 3' end.


=head1 OPTIONS

=over 30

=item B<--size>

[STR] Desired length of the read to trim to

=item B<--1>

[STR] Forward read file

=item B<--2>

[STR] Reverse read file (Optional)

=item B<--help>

Show this scripts help information.

=back


=cut


