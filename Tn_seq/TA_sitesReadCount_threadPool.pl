#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use threads;
use Thread::Queue;

my %options;
$options{'threads'} = 2;

GetOptions(\%options, 'sites=s', 'bam=s', 'threads=s', 'help|h') or die("Error in command line arguments\n");

if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
}

if(!$options{'sites'}){
	print STDERR " Error: Please provide the TA sites BED file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'bam'}){
	print STDERR " Error: Please provide sorted and indexed alignment BAM file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}



my @bedFields = qw(chr start end rest);
my @samFields = qw(qname flag chr pos mapq cigar rnext pnext tlen seq qual rest);

open(my $out, '>', 'result.bed') or $!;

# #Read the BED file and enqueue the TA sites in the processing queue
# open(my $fh, $options{'sites'}) or die "Cannot open file $options{sites}: $!";
	
# my $i = 1;
# while(eof($fh)){
	# my %parallel = ();
	
	# #Run in parallel for number of threads
	# foreach(1..$options{'threads'}){
		# if(!eof($fh)){
			# my $line = <$fh>;
			# chomp $line;
			# $parallel{$_}->{'bed'} = $line;
			# $parallel{$_}->{'thread'} = threads->create(\&getTAReadCount, $parallel{$_}->{'bed'});
			
			# $i++;
		# }
		# else{
			# last;
		# }
	# }
	
	# foreach(sort{$a <=> $b}keys %parallel){
		# $parallel{$_}->{'threadReturn'} = $parallel{$_}->{'thread'}->join();
		# my $result = join("\t", $parallel{$_}->{'bed'}, $parallel{$_}->{'threadReturn'});
		# print $result,"\n";
	# }
	
	
	# # if($i >= 50){
		# # last;
	# # }
# }

# close($fh);


#Link to understand the Perl thread pool template: http://www.perlmonks.org/?node_id=735923

my $qWork = new Thread::Queue();		#work queue
my $qResults = new Thread::Queue();		#result process queue

## Get the work items from the TA sites BED file
## and queue them up for the workers
my $builder=threads->create(\&buildQueue, $qWork);

sleep(10);

print STDERR localtime." Enqueuing to WORK thread started...\n";

## Create the pool of workers
my @pool = map{
	threads->create(\&worker, $qWork)
}1..$options{'threads'};
print STDERR localtime." All worker threads thread started...\n";


## Process the results as they become available
## until all the workers say they are finished.
# my $printer = threads->create(\&resultPrint, $qResults);
foreach(1..$options{'threads'}){
	my $i = 0;
	while(my $result = $qResults->dequeue()){
		print $out $result,"\n";
		$i++;
	}
	print STDERR localtime." $_ Printed $i records before encountering undef $_\n";
}

#Waiting for our threads to finish.
$builder->join();
$_->join() for @pool;
# $printer->join();

print STDERR localtime." ***All threads finished\n";

close($out);




# http://www.perlmonks.org/?node_id=1068673
# http://perldoc.perl.org/Thread/Queue.html
# http://www.perlmonks.org/?node_id=735923































sub buildQueue{
	my $qWork = shift @_;
	
	#Read the BED file and enqueue the TA sites in the processing queue
	open(my $fh, $options{'sites'}) or die "Cannot open file $options{sites}: $!";
	
	print STDERR localtime." Thread started enqueuing data into WORK queue...\n";

	my $i = 0;
	while(my $line = <$fh>){
		if($line !~ m/\S+\t\d+\t\d+\t\S+$/){
			next;
		}
		
		chomp $line;

		while(1){
			#Wait to see if our queue has < 20 items...
			if($qWork->pending() < 20){
				$qWork->enqueue($line);
				last;	#This breaks out of the infinite loop
			}
		}
		
		$i++;
	}
	
	print STDERR localtime." Enqueued total $i records in the WORK queue...\n";
	
	## Tell the workers there are no more work items
	# foreach(1..$options{'threads'}){
		# $qWork->enqueue(undef);
		# print STDERR localtime." Thread $_ enqueued undef in WORK queue...\n";
	# }
	$qWork->end();
}


sub worker{
	my $tid = threads->tid();
	my $qWork = shift @_;
	print STDERR localtime." Thread $tid started dequeuing and processing data from WORK queue...\n";
	
	my $i = 0;
	while(my $TA_site = $qWork->dequeue()){
		#print STDERR "Site:$tid: $TA_site\t";		#**************
		#Calculate the number of reads at TA site
		my $readCount = &getTAReadCount($TA_site);
		#print STDERR "Site:$tid: $TA_site\t$readCount\n";		#**************
		my $result = join("\t", $TA_site, $readCount);
		# print STDOUT $result,"\n";
		# print STDERR $result,"\n";

		#enqueue the reqult to the result queue
		$qResults->enqueue($result);
		$i++;
	}
	print STDERR localtime." Thread $tid processed total $i records...\n";
	## Signal this thread is finished
	$qResults->enqueue(undef);
	print STDERR localtime." Thread $tid enqueued undef in RESULT queue...\n";

}


sub resultPrint{
	my $qResults = shift @_;
	foreach(1..$options{'threads'}){
		print STDERR localtime." Thread $_ started dequeuing RESULT queue for printing...\n";
		while(my $result = $qResults->dequeue()){
			print $result,"\n";
		}
		print STDERR localtime." Thread $_ finished RESULT printing queue...\n";
	}
}


#Calculate the number of reads at a TA site
sub getTAReadCount{
	my %sam = ();
	my %bed = ();
	
	#Make the region string for 'samtools view ' command
	@bed{@bedFields} = split(/\t/, $_[0], 4);
	
	#Region width to consider the reads as a part of TA site insertion
	my $region = $bed{'chr'}.':'.($bed{'start'}-1).'-'.($bed{'end'}+1);
	my $TA_start = $bed{'start'}-1;
	my $TA_end = $bed{'end'}+1;
	
	# print "\n$region\n";
	my $readCount = 0;
	my $strand = 0;

	#foreach BAM record, decide if the read is valid TA site insertion proof
	foreach my $line(`samtools view $options{'bam'} $region`){
		#print $_;
		@sam{@samFields} = split(/\t/, $line, 12);
		
		#Check the read alignment strand
		if($sam{'flag'} & 16){
			$strand = '-';
		}
		else{
			$strand = '+';
		}
		
		#Find the end position of alignment on the reference sequence using CIGAR string
		my $refLen = 0;
		foreach($sam{'cigar'}=~m/(\d+\D)/g){
			$_=~m/(\d+)(\D)/;
			if($2 eq 'M' || $2 eq 'D'){
				#print $_,"\n";
				$refLen+=$1;
			}
		}
		
		my $endPos = $sam{'pos'} + $refLen - 1;
		
		# print join("\t", $strand, @sam{qw/qname flag chr pos cigar/}, $endPos),"\n";
		
		
		#Check if the current TA site is at the 5' end of the read. If yes, increase the read count for this TA site
		if($strand eq '+'){
			#For positive strand
			if($sam{'pos'} >= $TA_start && $sam{'pos'} <= $TA_end ){
				$readCount++;
			}
		}
		elsif($strand eq '-'){
			#For negative strand
			if($endPos >= $TA_start && $endPos <= $TA_end ){
				$readCount++;
			}
		}
		else{
			#something wrong
			print STDERR "Error at line:\n",$line;
			print STDERR "Alignment END position: $endPos\n";
			die;
		}
		
		$strand = 0;
		
	}
	
	#print STDERR "$_[0]\t$readCount\n";
	return $readCount;
}












__END__


=head1 NAME


=head1 SYNOPSIS

perl TA_sitesReadCount.pl --sites <TA sites bed file> --bam <input BAM file> --threads [number of threads]

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script finds the number of reads from input BAM file at each of the the TA sites from BED file.
Using Perl ithread pool to process the queue of data with constant number of threads.
Link to understand the Perl thread pool template: http://www.perlmonks.org/?node_id=735923


=head1 OPTIONS

=over 30

=item B<--sites>

[STR] BED file with TA sites. Only four columns allowed: (chr start end name)

=item B<--bam>

[STR] Input BAM file. This should be sorted and indexed

=item B<--threads>

[STR] Number of threads

=item B<--help>

Show this scripts help information.

=back


=cut






