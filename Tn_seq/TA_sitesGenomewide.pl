#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my %options;

GetOptions(\%options, 'fasta=s', 'help|h') or die("Error in command line arguments\n");

if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
}

# if(!$options{'fasta'}){
	# print STDERR "Error: Please provide the FASTQ file for the analysis\n";
	# &pod2usage({EXIT => 2, VERBOSE => 0});
# }


my $id = '';
my ($position, $i, $site) = (0, 0, 0);
my $TA_end = 0;


while(my $line = <>){
	chomp $line;
	
	#FASTA header line
	if($line =~ m/^>(\S+)/){
		$id = $1;
		$TA_end = 0;
		$position = 0;
	}
	elsif($line =~ m/^\w+$/){
		#If the TA site is split over two lines i.e. T at the end and A at the begining of next sequence line
		if($TA_end){
			if($line =~ m/^A/){
				$site = $position;
				print join("\t", $id, $site, $site+1, 'TA_'.$i),"\n";
				$i++;
			}
			$TA_end = 0;		#Cannot merge these nested if loops. $TA_end has to be set to 0 for this line because T should be present only in last line end.
		}
		
		#Normal TA site position finding
		while($line=~m/TA/g){
			$site = $position + pos($line) - 1;
			print join("\t", $id, $site, $site+1, 'TA_'.$i),"\n";
			$i++;
		}
		
		#If the current sequence line ends with T, set $TA_end flag to check for A at the begining of next line
		if($line =~ m/T$/){
			$TA_end = 1;
		}
		
		#Update the base position
		$position += length($line);
	}

}



__END__


=head1 NAME


=head1 SYNOPSIS

perl TA_sitesGenomewide.pl <reference/FASTA/file/path>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

TA_sitesGenomewide.pl: This script finds all the TA sites in the genome and report them as a BED region


=head1 OPTIONS

=over 30

=item B<--help>

Show this scripts help information.

=back


=cut






