#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $usage = 'This script parses the nucleotide genbank records and extracts 
the CGD IDs and GeneIds for each gene in gp file.
USAGE: perl parseGp.pl <genBank file>
';


if(!defined $ARGV[0] || !-e $ARGV[0]){
	print STDERR "ERROR: Input gp file not found\n\n";
	print STDERR $usage;
	die;
}


my $stream = Bio::SeqIO->new(-file   => $ARGV[0], -format => 'GenBank');

print join("\t", qw(accession GI cgdId geneId locusTag)),"\n";

while( my $seqObject = $stream->next_seq() ) {
	my $accession = $seqObject->accession_number();
	my $gi = $seqObject->primary_id();

	
	my ($cgdId, $geneId, $locus) = ('', '', '');
	my @locusTags = ();
	
	for my $featObject ($seqObject->get_SeqFeatures()) {
		
		## get CGD cross-reference
		if($featObject->primary_tag() eq 'CDS'){

			if($featObject->has_tag('db_xref')){
				
				foreach my $tagVal($featObject->get_tag_values('db_xref')){
					if($tagVal =~m/CGD:(\w+)/){
						$cgdId .= $1.';';
					}
				}
				
			}
			
		}
		## get GeneId and locus_tag 
		elsif($featObject->primary_tag() eq 'gene'){
			if($featObject->has_tag('locus_tag')){
				@locusTags = $featObject->get_tag_values('locus_tag');
			}
			
			if($featObject->has_tag('db_xref')){
				foreach my $tagVal($featObject->get_tag_values('db_xref')){
					if($tagVal =~m/GeneID:(\w+)/){
						$geneId .= $1.';';
					}
				}
			}
		}
		else{
			next;
		}
	}
	
	if($#locusTags >= 0){
		$locus = join(';', @locusTags);
	}
	
	chop $cgdId;
	chop $geneId;
	
	print join("\t", $accession, $gi, $cgdId, $geneId, $locus),"\n";
	
	# last;
}














