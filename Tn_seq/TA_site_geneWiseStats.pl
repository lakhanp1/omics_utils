#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::HTS::Tabix;


my %options;


GetOptions(\%options, 'geneBed=s', 'TA_counts=s', 'threads=s', 'help|h') or die("Error in command line arguments\n");

if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
}

if(!$options{'geneBed'}){
	print STDERR "Error: Please provide the gene regions as BED file format. The name column should have gene name\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}

if(!$options{'TA_counts'}){
	print STDERR "Error: Please provide the genome wide TA counts generated using script TA_sitesReadCount.pl\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}



# Check for the bgzipped file and tabix index
if($options{'TA_counts'} !~ m/gz$/){
	print STDERR $options{'TA_counts'}," file is not bgzipped\n";
	die;
}

if(!-e "$options{'TA_counts'}.tbi"){
	print STDERR "TABIX index not found for the file $options{'TA_counts'}\n";
	die;
}




my @bedFields = qw(chr start end name rest);
my @TAcountFields = qw(chr start end name count);
my @TAstatsFields = qw(TA_siteCount TA_siteDisrupted TA_siteReads centralTA_siteCount centralTA_siteDisrupted centralTA_siteReads);

my $tabix = Bio::DB::HTS::Tabix->new( filename => $options{'TA_counts'} );



print '#',join("\t",qw(chr start end name), @TAstatsFields),"\n";

# read the gene file
open(my $gFh, $options{'geneBed'}) or die "Cannot open file $options{'geneBed'}: $!";

while(my $line = <$gFh>){
	chomp $line;
	my %record = ();
	@record{@bedFields} = split(/\t/, $line, 5);
	my $region = $record{chr}.':'.$record{start}.'-'.$record{end};
	# print $region,"\t", $record{rest},"\n";
	
	# new start and new end for the central 80% region
	my $newStart = $record{start} + int( ($record{end} - $record{start})/10 );
	my $newEnd = $record{end} - int( ($record{end} - $record{start})/10 );
	
	# extract TA count stats
	my $statsHref = &TA_stats($region, $newStart, $newEnd);
	
	print join("\t", @record{qw(chr start end name)}, $statsHref->{TA_siteCount}, $statsHref->{TA_siteDisrupted}, $statsHref->{TA_siteReads}, $statsHref->{centralTA_siteCount}, $statsHref->{centralTA_siteDisrupted}, $statsHref->{centralTA_siteReads}),"\n";
	
}


close($gFh);
$tabix->close;


















sub TA_stats{
	my $centralEnd = pop @_;
	my $centralStart = pop @_;
	
	# extract the TA site counts for this region
	my $iter = $tabix->query($_[0]);
	
	my %TA_stats = ();
	@TA_stats{@TAstatsFields} = (0,0,0,0,0,0);
	
	while ( my $ta = $iter->next() ) {
		# print $ta,"\n";
		my %taSite = ();
		# (chr start end name count)
		@taSite{@TAcountFields} = split(/\t/, $ta);
		
		# (TA_siteCount TA_siteDisrupted TA_siteReads centralTA_siteCount centralTA_siteDisrupted centralTA_siteReads)
		$TA_stats{TA_siteCount}++;
		
		if($taSite{count} > 0){
			$TA_stats{TA_siteReads} += $taSite{count};
			$TA_stats{TA_siteDisrupted}++;			
		}
		
		# for the TA sites lying in central 80% gene region
		if(($taSite{start} >= $centralStart ) && ($taSite{end} <= $centralEnd)){
			$TA_stats{centralTA_siteCount}++;
			
			if($taSite{count} > 0){
				$TA_stats{centralTA_siteReads} += $taSite{count};
				$TA_stats{centralTA_siteDisrupted}++;
			}
		}
		
	}
	
	
	
	return \%TA_stats;
}














__END__


=head1 NAME


=head1 SYNOPSIS

perl TA_site_geneWiseStats.pl --geneBed <gene bed file> --TA_counts <genome wide TA counts> 

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script calculates following statistics related to TA site read count for each gene.
1. Number of TA sites
2. Number of TA reads
3. Number of TA sites in the central 80% of the gene
4. Number of TA reads in central 80% region of the gene

Required tools: htslib, samtools

=head1 OPTIONS

=over 30

=item B<--geneBed>

[STR] gene regions as BED file format. The name column should have gene name (chr start end name)

=item B<--TA_counts>

[STR] genome wide TA counts generated using script TA_sitesReadCount.pl

=item B<--help>

Show this scripts help information.

=back


=cut






