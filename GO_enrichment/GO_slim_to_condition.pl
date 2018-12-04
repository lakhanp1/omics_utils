#!/usr/bin/perl -w

use strict;
use warnings;
use Spreadsheet::XLSX;

my $usage = 'This script parse the Genes -> GO slim master file. The GO slim master file has to be created from Amigo GO slimmer tool.
For the given excel file with differential list of gene in Up and Down worksheet, this script calculates the number of genes falling in each GO slim category.
Usage: perl GO_slim_to_condition.pl <Genes_to_GO_slim_master_file> <Diff_genes_excelx_file> <gene_id_mapping_file>
';

if(@ARGV != 3){
	print STDERR $usage;
	die;
}


my $goSlim_masterFile = $ARGV[0];
my $diffGeneXlsxFile = $ARGV[1];
my $ensemblToZfinFile = $ARGV[2];	#This file has ENSEMBL ID to ZFIN ID mapping for genes

my @goSlimFilelds = qw(GO_ID Term Count Genes GO_category);
my %geneToSlim = ();

#Read the output of Amigo GO Slimmer tool for all the genes and store the GO slim information for each gene
print STDERR "Read the output of Amigo GO Slimmer tool for all the genes and store the GO slim information for each gene\n";

open(my $fh1, $goSlim_masterFile) or die "Cannot open file $goSlim_masterFile: $!";

my %goSlim = ();
while(<$fh1>){
	chomp;
	@goSlim{@goSlimFilelds} = split(/\t/, $_);
	
	$goSlim{'GO_category'} =~ s/^\s*(\S+)\s*/$1/;
	my $GO = join('~', $goSlim{'GO_ID'}, $goSlim{'Term'});
	
	foreach my $zfinGene(split(/\s+/, $goSlim{'Genes'})){
		$zfinGene=~s/^\s*ZFIN:(\S+)\s*/$1/;
		
		# {ZFin gene id} -> {P|C|F} -> {GO_ID~GO_slim term}
		push(@{$geneToSlim{$zfinGene}->{$goSlim{'GO_category'}}}, $GO);
	}
}


# foreach my $gene(keys %geneToSlim){
	# foreach my $cat(keys %{$geneToSlim{$gene}}){
		# print join("\t", $gene, $cat, @{$geneToSlim{$gene}->{$cat}}),"\n";
	# }
# }



print STDERR "Done\n";

# Read the ENSEMBL gene ID -> ZFIN gene ID map file and store the data
print STDERR "Read the ENSEMBL gene ID -> ZFIN gene ID map file and store the data\n";
open(my $fh2, $ensemblToZfinFile) or die "Cannot open file $ensemblToZfinFile: $!";

my %ensemblToZfin = ();

while(<$fh2>){
	chomp;
	next if /^#/;
	my @temp = split(/\t/, $_);
	
	$ensemblToZfin{$temp[0]} = $temp[2];
}

print STDERR "Done\n";

# Read the differential genes excel file
# for the genes in Up and Down worksheets, assign the GO slim 
# make a GO slim count table for Up regulated and Down regulated genes
# {GO} -> [{name}, {up_count}, {down_count}]

print STDERR "Read differentially expressed gene list XLSX file and make GO Slim count table\n";

my %goSlimTable = ();

my $excel = Spreadsheet::XLSX->new($diffGeneXlsxFile);


my $i=0;

for my $sheet (@{$excel -> {Worksheet}}) {
		
	if($sheet->{Name} eq 'all' || $sheet->{Name} eq 'filtered'){
		next;
	}
	
	# printf("Sheet: %s\n", $sheet->{Name});
	
	my $DEG_class = $sheet->{Name};
	$sheet -> {MaxRow} ||= $sheet -> {MinRow};
	# print $sheet -> {MaxRow},"\t",$sheet -> {MinRow},"\t",$sheet->{Name},"\n";
	
	for my $row ( $sheet -> {MinRow} .. $sheet -> {MaxRow} ){
		$sheet -> {MaxCol} ||= $sheet -> {MinCol};
		
		my $cell = $sheet -> {Cells}[$row][1];
		next unless $cell;
		
		# print "Value = ",$cell->value(),"\t",$DEG_class,"\n";
		# get the ZFIN gene ID for current
		my $ensemblId = $cell->value();
		my $zfinId = exists $ensemblToZfin{$ensemblId} ? $ensemblToZfin{$ensemblId} : '-';
		
		$i++;
		# print "$i\t$ensemblId\t$zfinId\t$DEG_class\n";
		
		# %goSlimTable = {GO_category(P|C|F)} -> {GO_slim_term} -> {up|down} -> count
		if(defined $geneToSlim{$zfinId}){
			
			# print "$i\t$ensemblId\t$zfinId\t$DEG_class\n";
			
			foreach my $category(qw(P C F)){	
				
				if(exists $geneToSlim{$zfinId}->{$category}){
					# print "Found: $ensemblId\t$zfinId\t$category\t$DEG_class\n";
					
					foreach my $GO_term(@{$geneToSlim{$zfinId}->{$category}}){
						
						if(!exists $goSlimTable{$category}->{$GO_term}->{'up'}){
							$goSlimTable{$category}->{$GO_term}->{'up'} = 0;
						}
						if(!exists $goSlimTable{$category}->{$GO_term}->{'down'}){
							$goSlimTable{$category}->{$GO_term}->{'down'} = 0;
						}
						
						$goSlimTable{$category}->{$GO_term}->{$DEG_class}++;
						# print join("\t", $ensemblId, $zfinId, $category, $GO_term),"\n";
					}
				}
			}
		}
	}
}

# print "Number of rows read: $i\n";

print "#",join("\t", qw(GO_slim_term up_regulated_genes down_regulated_genes category)),"\n";
foreach my $category(keys %goSlimTable){
	foreach my $GO_slim(sort keys %{$goSlimTable{$category}}){
		print join("\t", $GO_slim, $goSlimTable{$category}->{$GO_slim}->{'up'}, $goSlimTable{$category}->{$GO_slim}->{'down'}, $category),"\n";
	}
}






close($fh1);
close($fh2);

