#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Excel::Writer::XLSX;

my %diffFiles;
my $i=0;
my @header=();

#perl /home/lakhanp/scripts/utils/splitCuffdiffSamples.pl --geneDiff gene_exp_annotations.diff --qval 0.05 --log2 1

my %options;

GetOptions(\%options, 'geneDiff=s', 'qval=f', 'log2=f', 'help|h') or die("Error in command line arguments\n");

if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
}


if(!$options{'geneDiff'} || !$options{'qval'} || !$options{'log2'}){
	print STDERR "Error: Please provide all the arguments\n";
	&pod2usage({EXIT => 2, VERBOSE => 1});
}


open(my $fh, $options{'geneDiff'}) or die "Cannot open file $options{geneDiff}: $!";


while(<$fh>){
	chomp;
	if($i==0){
		@header = split(/\t/, $_);
		$i++;
		next;
	}

	my @temp = split(/\t/, $_);
	my $pair = $temp[4].'_vs_'.$temp[5];
	if(!exists $diffFiles{$pair}){
		
		$diffFiles{$pair}->{'workbook'} = Excel::Writer::XLSX->new($pair.".xlsx");
		
		#create new Excel file and create 4 worksheets in that
		foreach my $fil(qw(all up down filtered)){
			$diffFiles{$pair}->{$fil}->{'worksheet'} = $diffFiles{$pair}->{'workbook'}->add_worksheet($fil);
			$diffFiles{$pair}->{$fil}->{'row'} = 0;
			$diffFiles{$pair}->{$fil}->{'worksheet'}->write_row($diffFiles{$pair}->{$fil}->{'row'}, 0, \@header);
			$diffFiles{$pair}->{$fil}->{'row'}++;
		}
				
	}
	
	# write to all worksheet
	$diffFiles{$pair}->{'all'}->{'worksheet'}->write_row($diffFiles{$pair}->{'all'}->{'row'}, 0, \@temp);
	$diffFiles{$pair}->{'all'}->{'row'}++;
	
	
	#print {$diffFiles{$temp[4].'_vs_'.$temp[5]}} $_;
	#filter the rows with status !OK and q_value > cutoff
	if($temp[6] ne 'OK' || $temp[12] > $options{'qval'}){
		$diffFiles{$pair}->{'filtered'}->{'worksheet'}->write_row($diffFiles{$pair}->{'filtered'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'filtered'}->{'row'}++;
		next;
	}
	
	#write the -/+inf rows to down/up sheets respectively
	if($temp[9] eq 'inf'){
		$diffFiles{$pair}->{'up'}->{'worksheet'}->write_row($diffFiles{$pair}->{'up'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'up'}->{'row'}++;
		next;
	}
	elsif($temp[9] eq '-inf'){
		$diffFiles{$pair}->{'down'}->{'worksheet'}->write_row($diffFiles{$pair}->{'down'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'down'}->{'row'}++;
		next;
	}
	
	#up regulated
	if($temp[9] >= $options{'log2'}){
		$diffFiles{$pair}->{'up'}->{'worksheet'}->write_row($diffFiles{$pair}->{'up'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'up'}->{'row'}++;
	}
	#down regulated
	elsif($temp[9]*-1 >= $options{'log2'}){
		$diffFiles{$pair}->{'down'}->{'worksheet'}->write_row($diffFiles{$pair}->{'down'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'down'}->{'row'}++;
	}
	#rows where qval is < threshold but log2 threshold is failed
	else{
		$diffFiles{$pair}->{'filtered'}->{'worksheet'}->write_row($diffFiles{$pair}->{'filtered'}->{'row'}, 0, \@temp);
		$diffFiles{$pair}->{'filtered'}->{'row'}++;
		next;
	}
	
}


close($fh);







__END__


=head1 NAME


=head1 SYNOPSIS

This script parse the gene_exp.diff file and prepares different MS-Excel files for pairwise sample comparison. Inside the Excel file, it creates three sheets: all data, up regulated and down regulated genes.

USAGE: perl splitCuffdiffSamples.pl --geneDiff <PATH_TO_GENES_DIFF> --log2 <LOG2_CUTOFF> --qval <Q_VALUE_CUTOFF>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

Split the gene_exp.diff file into pairwise samples.

=head1 OPTIONS

=over 30

=item B<--help>

Show this scripts help information.

=item B<--geneDiff>

[STR] Complete path to gene_exp.diff output file from cuffdiff

=item B<--log2>

[FLOAT] log2 fold change cutoff

=item B<--qval>

[FLOAT] q-value cutoff

=back


=cut




