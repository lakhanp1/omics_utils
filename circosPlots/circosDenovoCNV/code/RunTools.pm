package RunTools;

use strict;
use warnings;

use File::Copy;
require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw(Circos blast);
our @EXPORT = qw();
our $VERSION = 1.1;


#run CIRCOS executible
sub Circos{
	my $cirExe = shift;
	my $circosConfFile = shift;
	my $kmer = shift;
	my $imageMapFile = shift;
	my $circosEtcPath = shift;
	my $dashboardDir = shift;
	
	my $imageName = $circosEtcPath."/circos_".$kmer.".png";
	my $svgName = "circos_".$kmer.".svg";
	
	if($kmer == -1){
		$imageName = $circosEtcPath."/circos.png";
		$svgName = "circos.svg";
	}
	
	#print $cirExe,"\n",$circosConfFile,"\n";
	print "Running circos with command: $cirExe -conf $circosConfFile -outputfile $imageName -debug_group summary,timer > run.out\n";
	my $result = `$cirExe -conf $circosConfFile -outputfile $imageName -debug_group summary,timer > circosRunLog.out`;
	
	#move the svg file to dashboard_output folder
	move("$circosEtcPath/$svgName", $dashboardDir) or die "Cannot copy $svgName to $dashboardDir: $!";
	move($imageName, $dashboardDir) or die "Cannot copy $imageName to $dashboardDir: $!";
	
	#write svg file name and kmer value to circos_output.result file
	open(my $outAppend, ">>", $imageMapFile) or die "Cannot open file $imageMapFile to append data: $!";
	if($kmer == -1){
		print $outAppend "$svgName\n";
	}
	else{
		print $outAppend "$kmer,$svgName\n";
	}
	close($outAppend);
}



#system("$blastPath/formatdb -i $header -o T -p F");
#system("$blastPath/blastall -p blastn -d $Reference -i $contig -m8 >blastout.txt");
sub blast{
	my $blastDir = shift;
	my $refSeq = shift;
	my $querySeq = shift;
	my $cpus = shift;
	
	$querySeq =~m/(.*)\./;
	my $resultFile = $1.".blastout";
	
	#formatdb
	my $blastDb = `$blastDir/formatdb -i $refSeq -p F -o T`;
	
	#blastall
	my $blastRun = `$blastDir/blastall -p blastn -d $refSeq -i $querySeq -m 8 -a $cpus -o $resultFile > blastall.log`;
	
	return $resultFile;
}



1;