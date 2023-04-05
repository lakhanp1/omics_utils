#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = "This script removes the tRNA and its associated gene, exons from GFF file
perl filter_tRNA_GFF.pl <GFF>";


open(my $fh, $ARGV[0]) or die "Cannot open file $ARGV[0]: $!";

my @record = ();

my ($id, $pk, $count, $filterCount) = ('', '', 0, 0);
my $filter = 0;
my $firstRecordFlag = 1;

while(<$fh>){
	if(/^#/){
		print $_;
		next;
	}
	
	if(/\tID=([^;]+);/){
		$id=$1;
		
		if(/\ttRNA\t/){
			$filter = 1;
			$filterCount++;
		}
		
		#First entry in the record i.e. gene row
		if(!/;Parent=/){
			$pk = $id;
			$count++;
			#process the last saved record
			if(!$firstRecordFlag){
				if(!$filter){
					print @record;
				}
			}
			else{
				$firstRecordFlag = 0;
			}
			
			#empty the @record to start saving new record
			@record = ();
			$filter = 0;

			
		}
		
		push(@record, $_);

	}
	else{
		print "New pattern found:\n",$_;
		die;
	}
	
}

if(!$filter){
	print @record;
}

print STDERR "Total features: $count\n";
print STDERR "Total filtered features: $filterCount\n";
