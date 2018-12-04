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

if(!$options{'fasta'}){
	print STDERR "Error: Please provide the FASTQ file for the analysis\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


while(<>){
	chomp;
	next if /(^#|^\s*$)/;
	
	my ($chr, $start, $end, $gene) = split(/\t/, $_);
	
	my $fa = `samtools faidx $options{'fasta'} $chr:$start-$end`;
	
	if($? != 0){
		print STDERR "samtools faidx failed for region: $chr:$start-$end\n";
		die;
	}
	else{
		&TA_sites($fa, $chr, $start, $gene);
	}
}


# my ($chr, $start, $end, $gene) = qw(NC_004603.1 371 805 gene);

# my $fa = ">NC_004603.1:371-805
# TTAAAAGCGCTCGATATTCTCTTTTAACCATGCTTCTGCAGCATCTTCTGGGACGTCATG
# AGCCAGTACATCAATCAACAAGCAATCTGCAATCGGCGTTGCGCCAATATCTTCTAACAA
# GTTGTAAGCGTGCTTACCGGCAGCGCAGAACGTATCGTAGCTCGAATCACCAATCGCAAT
# GACTGCAAACTTCACGTCGGCCATCTTAGGTGGCGTATTTTGCAATGCTGCGATGAAAGG
# TTGAATATTGTCTGGATATTCCCCCGCACCGTGGGTGGAGGTGATCACAAGCCAAGTCCC
# TTGGTTATCGATCTCGTCTAATGACGGTTGGTTGTGGATCGTGGTTTCAAATCCTTGTTC
# AACGAGCAGATCGCTCAGGTGATCACCGACGTATTCCGCACCACCAAGAGTGCTACCTGT
# AATAATGTGGATCAT
# ";
# &TA_sites($fa, $chr, $start, $gene);

sub TA_sites{
	my $id = '';
	my $seq = '';
	#my ($fa, $chr, $start, $end) = @_;
	
	foreach my $line(split(/\n/, $_[0])){
		chomp $line;
		if($line =~ m/^>/){
			$id = $line;
		}
		elsif($line =~ m/^\w+$/){
			$seq .= $line;
		}
	}
	
	my $i=1;
	while($seq=~m/TA/g){
		my $site = pos($seq)+$_[2]-2;
		print join("\t", $_[1], $site, $site+1, $_[-1].'_TA'.$i),"\n";
		$i++;
	}
	
}




__END__


=head1 NAME


=head1 SYNOPSIS

perl TA_sitesPerGene.pl --fasta <reference/FASTA/file/path>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

TA_sitesPerGene.pl


=head1 OPTIONS

=over 30

=item B<--fasta>

[STR] Reference FASTA file path

=item B<--help>

Show this scripts help information.

=back


=cut






