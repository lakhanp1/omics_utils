#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = "This script read the annotation/gene description file and writes the description for each gene in cuffdiff output as a last column.
USAGE: perl addAnnotationToCuffdiffOut.pl <cuffdiff_out_file> <gene_description_file> <Out_file>\n";

if(@ARGV != 3){
	die $usage;
}


#store the gene description information in data structure
open(my $descFh, $ARGV[1]) or die "Cannot open file $ARGV[1]: $!";

my %data = ();
my %idRef = ();

my $i = 1;
while(<$descFh>){
	if(/^WBGene/){
		chomp;
		my @tmp = split(/\t/, $_);
		
		if($tmp[3] eq 'none available'){
			if($tmp[4] eq 'none available'){
				if($tmp[5] eq 'none available'){
					if($tmp[6] eq 'none available'){
						$data{$i} = $tmp[7];
					}
					else{
						$data{$i} = $tmp[6];
					}
				}
				else{
					$data{$i} = $tmp[5];
				}
			}
			else{
				$data{$i} = $tmp[4];
			}
		}
		else{
			$data{$i} = $tmp[3];
		}
		
		$idRef{$tmp[0]} = $i;
		$idRef{$tmp[1]} = $i;
		$idRef{$tmp[2]} = $i;
		
		$i++;
	}
	
	
}

open(my $out, '>',$ARGV[2]) or die "Cannot open file $ARGV[2]: $!";
open(my $fh, $ARGV[0]) or die "Cannot open file $ARGV[0]: $!";

while(<$fh>){
	chomp $_;
	
	if(/^(Gene|Pseudogene):(\S+):?/){
		print $out $_,"\t";
		my $gene = $2;
		
		if(defined $idRef{$gene}){
			print  $out $data{$idRef{$gene}},"\n";
		}
		else{
			print $out "Not known\n";
		}
	}
	elsif(/^test_id/){
		print $out $_,"\tGene_annotation\n";
	}
	else{
		print $out $_,"\tNot known\n";
	}
}



close($descFh);
close($fh);
close($out);
