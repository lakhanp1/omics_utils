#!/usr/env perl

use strict;
use warnings;
use File::Copy;

my $usage = "Usage: perl merge_fastqs.pl work/dir/path raw_data_file\n";
	
if($#ARGV != 1){
	print STDERR "Error: wrong arguments\n$usage";
	exit 1;
}

my $workDir = $ARGV[0];

open(my $fh, $ARGV[1]) or die "File not found: $!\n";

my $firstLine = 1;
my @colNames = ();
my $colCount = undef;
my %rawFiles = ();
my %fqCopied = ();

while(my $line = <$fh>){
	chomp $line;
	next if $line=~m/^$/;
	
	my @data = split(/\t/, $line);
	my %fqInfo = ();

	if($firstLine){
		@colNames = @data;
		$colCount = $#colNames;
		%fqInfo = map {$_ => 1} @colNames;
		# print join(" || ", @colNames), "\n";
		$firstLine = 0;
		
		if(!exists($fqInfo{sampleId}) || !exists($fqInfo{sampleName}) || !exists($fqInfo{R1}) || !exists($fqInfo{outDir})){
			die "Error: Missing one or more columns of sampleId/sampleName/R1/outDir\n";
		}
		
		next;
	}

	if($#data != $colCount){
		die "Error: Value count does not match in row\n";
	}
	
	@fqInfo{@colNames} = @data;
	# print join(" || ", @fqInfo{qw(sampleId R1 outDir)}), "\n";
	
	my $outPath = "$workDir/$fqInfo{outDir}";
	if(! -d $outPath){
		mkdir($outPath);
	}
	

	if(! -e $fqInfo{R1}){
		die "Error: File not found: $fqInfo{R1}\n";
	}
	else{
		## a raw file is not expected twice.
		if(exists($rawFiles{$fqInfo{R1}})){
			die "Error: $fqInfo{R1} being used repeatedly...\n";
		}
		
		
		$rawFiles{$fqInfo{R1}}++;
		
		my $newFqFile = "$outPath/$fqInfo{sampleName}_R1.fastq.gz";
		
		if(-e $newFqFile){
			## make sure that file is being appended to recently copied file
			if(!exists($fqCopied{$newFqFile})){
				die "Error: Older version of file unrelated to current processing already exists: $newFqFile\n";
			}
			
			# append
			my $command = "zcat $fqInfo{R1} | gzip -c >> $newFqFile";
			system($command) == 0 or die "Error: Command failed: $command: $!\n";
			print "Append: $fqInfo{R1} $newFqFile\n";
		}
		else{
			# copy
			copy($fqInfo{R1}, $newFqFile) or die "Copying failed: $!";
			print "Copy: $fqInfo{R1} $newFqFile\n";
			$fqCopied{$newFqFile}++;
		}
		
	}
	
	
}

close($fh);


