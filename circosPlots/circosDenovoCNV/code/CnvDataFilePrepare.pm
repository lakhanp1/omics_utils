package CnvDataFilePrepare;


use strict;
use warnings;
use File::Path;
use File::Copy; 


require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(cnvGlyphFile);
our @EXPORT = qw();
our $VERSION = 1.1;


#create the glyph data file for circos from BED file
sub cnvGlyphFile{
	my $cnvFile = shift @_;
	my $storeLocation = shift @_;
	my $idMap = shift @_;
	
	#create the circos/data directory to store the histogram.txt and links.txt file
	if(!-d $storeLocation){
		mkpath($storeLocation, 1, 0777) or die "Cannot create path $storeLocation: $!";
	}
	
	my $glyphFile = $storeLocation."/cnvGlyph.txt";
	
	open(my $out, ">", $glyphFile) or die "Cannot crete file $glyphFile: $!";
	open(my $in, $cnvFile) or die "Cannot open file $cnvFile: $!";
	
	my ($longestCnv, $shortestCnv) = (undef, undef);
	while(<$in>){
		if(/(^\w+)\s+(\d+)\s+(\d+)/){
			chomp $_;
			my $cnvSize = $3-$2;
			my $chr = $1;
			#update the longest CNV and shortest CNV
			if(!defined $longestCnv){
				$longestCnv = $cnvSize;
				$shortestCnv = $cnvSize;
			}
			elsif($cnvSize > $longestCnv){
				$longestCnv = $cnvSize;
			}
			elsif($cnvSize < $shortestCnv){
				$shortestCnv = $cnvSize;
			}
			
			$_.="\tlabel_size=".$cnvSize.",color=".$idMap->{$chr}->[0];
			s/\s+/ /g;
			
			print $out $_,"\n";
		}
	}
	
	
	close($in);
	close($out);
	
	return ($glyphFile, $shortestCnv, $longestCnv);
}


1;