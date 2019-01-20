#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $usage = 'This script parses the protein genbank records and extracts 
the CGD IDs and GeneIds for each protein in gp file.
USAGE: perl parseGb.pl <genPept file>
';


if(!defined $ARGV[0] || !-e $ARGV[0]){
	print STDERR "ERROR: Input gp file not found\n\n";
	print STDERR $usage;
	die;
}


my $stream = Bio::SeqIO->new(-file   => $ARGV[0], -format => 'GenBank');

print join("\t", qw(accession GI dbId geneId locusTag)),"\n";

while( my $seqObject = $stream->next_seq() ) {
	my $accession = $seqObject->accession_number();
	my $gi = $seqObject->primary_id();

	
	my ($cgdId, $geneId, $locus) = ('', '', '');
	my @locusTags = ();
	
	for my $featObject ($seqObject->get_SeqFeatures()) {
		
		## get information only from CDS tag
		if($featObject->primary_tag() eq 'CDS'){
		
			if($featObject->has_tag('locus_tag')){
				@locusTags = $featObject->get_tag_values('locus_tag');
			} 
			
			if($featObject->has_tag('db_xref')){
				
				foreach my $tagVal($featObject->get_tag_values('db_xref')){
					if($tagVal =~m/SGD:(\w+)/){
						$cgdId .= $1.';';
					}
					elsif($tagVal =~m/GeneID:(\w+)/){
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














