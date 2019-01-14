#!/usr/bin/perl -w

#1) This scripts runs BLAST on the assembled contigs.fa against user selected REFERENCE 
#2) parse the BLAST output to get the frequency of contig coverage on different intervals on chromosomes AND
#3) get the information if a contig or some of its part has mapped to two regions on different chromosomes
#4) output this info in the data file format for circos
#5) on the fly generate the config files for the circos based on current data using already created templates and TemplateToolkit2
#6) run the circos and copy the image to dashboradOutput folder and create the CSV 

#How to run:
#perl circosAutomation.pl --refFile ../workDir/genome/genome.fa --blastPath /LAB3/Sangenix/Install/Ngs_Tools/blast-2.2.25/bin --queryFile ../workDir/genome/yeast_sequence.fasta --circosPath /home/vineet/circos-0.66/bin/circos --confTemplate ../centralData/confTemplates --histCutOffs 70 80 90 --linkThreshold 1000 --linkQueryCov 80 --kmer 23 --proc 8 

#perl circosAutomation.pl --refFile <REFERENCE_FILE_AGAINST_WHICH_BLAST_WILL_BE_RUN> --blastPath <BLAST_PATH_UPTO_BIN> --queryFile ../workDir/genome/yeast_sequence.fasta --circosPath <CIRCOS_EXE_FILE_COMPLETE_PATH> --confTemplate <DIR_PATH_FOR_CIRCOS_CONF_TEMPLATES> --histCutOffs <HISTOGRAM_%_IDENTITY_CUTOFFS_SPACE_SPEARATED> --linkThreshold <INT_VALUE_FOR_LINK_THRESHOLD> --linkQueryCov <QUERY_COVERAGE_VALUE_FOR_LINK> --kmer <UNIQ_IDENTIFIER_TO_ADD_IN_POSTPROCESS_FILE> --proc <NUMBER_OF_THREADS>

use strict;
use warnings;
use Getopt::Long;

use DataFilePrepare qw(karyotypeData histogramAndLinks setParameters getQueryLength %circosFiles setChrUnit);
use UpdateConfigs qw(updateDenovoCircosConf updateTicksConf copyOtherConfs);
use RunTools qw(Circos blast);


#command line inputs
my $workflowPath = ".";
my $refFile = "";
my $queryFile = "";
my $blastPath = "";
my $circosPath = "";
my $confTemplate = "";
my $histCutOffs = [70, 80, 90];			#default
my $linkThreshold = 1000;			#default
my $linkQueryCov = 80;			#default
my $kmer = undef;
my $proc = 8;			#default

my $histCoverage = 51;			#This is hardcoded
my $histInterval = 10000;			#this is determined internally based on the chromosome length



GetOptions(
	"refFile=s" => \$refFile,
	"blastPath=s" => \$blastPath,
	"queryFile=s" => \$queryFile,
	"circosPath=s" => \$circosPath,
	"confTemplate=s" => \$confTemplate,
	'histCutOffs=i@{1,3}' => \$histCutOffs,
	"linkThreshold=i" => \$linkThreshold,
	"linkQueryCov=i" => \$linkQueryCov,
	"kmer=i" => \$kmer,
	"proc=i" => \$proc)
or die "Error in command line arguments";
	


#remove the trailing / from the command line arguments
foreach($blastPath, $confTemplate){
	if($_ =~ m/\/$/){
		chop $_;
	}
}

#print "workflowPath=s\t", $workflowPath,"\n";
#print "refFile=s \t",$refFile,"\n";
#print "blastPath=s\t",$blastPath,"\n";
#print "queryFile=s\t",$queryFile,"\n";
#print "circosPath=s\t",$circosPath,"\n";
#print "confTemplate=s\t",$confTemplate,"\n";
#print "histCutOffs=i{1,3}\t@$histCutOffs\n";
#print "linkThreshold=i\t",$linkThreshold,"\n";
#print "linkQueryCov=i\t",$linkQueryCov,"\n";


=head
	--refFile = [reference file path which is browsed by user and will be used to generate Karyotype and as BLAST DB]
	--blastPath = [BLAST exe path]
	--queryFile = [file path which will be used as query sequence in BLAST, assembled contigs]
	--circosPath = [path for ciccos_exe]
	--confTemplate = [template file path for the circos conf files]
	--histInterval = [bin size of the histogram]		[optional-need to handle default]
	--histCutOffs = [The %identity cutoffs for plotting different histograms]		[optional-need to handle default]
	--linkThreshold = [Threshold query length to be matched on CHRs to establish a link. can be in %]		[optional-need to handle default]
	--linkQueryCov = [Minimum query coverage to establish a link]		[optional-need to handle default]
	--kmer = [kmer length used in the assembly algorithm]
	--proc = [number of processors to use for BLASTALL]
	--pmPath = [path where all *.pm files for circos are kept]
	
	NOT USED
	--histCoverage = [query sequence length to have min % in a bin to add its frequency]		[optional-need to handle default]	
	--workflowPath = [path for current workflow]
=cut










################Line to comment all the code#######################

################Line to comment all the code#######################

my $dataStorePath = $workflowPath."/circos_"."$kmer/data";
my $confStorePath = $workflowPath."/circos_"."$kmer/etc";

#this is a csv file with (kmer value, circos image for that kmer)
my $dashboardDir = "$workflowPath/dashboard_output";
my $imageMapFile = "$dashboardDir/circos_output.result";

if(!-d $dashboardDir){
	mkdir($dashboardDir, 0777) or die "Cannot create directory $dashboardDir: $!";
}


#circos runtime parameters
my %config;
#1) The chromosomes_unit value is used as a unit (suffix "u") to shorten values in other parts of the configuration file.
$config{chromosomes_units} = 100000;

#2) the tick label is derived by multiplying the tick position by 'multiplier' and casting it in 'format'
$config{multiplier} = "1e-5";
$config{linkRadius} = 0.65;


#create karyotype file
my $idMap = &karyotypeData($refFile, $dataStorePath);


#decide the 'chromosomes_unit', 'multiplier' and 'histInterval' parameters based on genome size so that the ticks and tick lables will be appropriate
($config{chromosomes_units}, $config{multiplier}) = &setChrUnit();
$histInterval = $config{chromosomes_units} * 2;

#set the parameters in DataFilePrepare.pm
&setParameters($histCoverage, $histInterval, $histCutOffs, $linkQueryCov, $linkThreshold);


#read the QUERY file and get the length for each query seq. the lengths for each of the query sequences is stored in DataFilePrepare:%queryLen
#DataFilePrepare:%queryLen will be accessed for calculating the query coverage for plotting HISTOGRAM and LINKS in the plot
&getQueryLength($queryFile, $dataStorePath);



#create empty bins for the histogram
#IX => [[chr5	0	24999	0	0	0], [chr5	25000	49999	0	0	0], ....]
my %histogram;
my @histFreq = split(//, "0"x@{$histCutOffs});

foreach(keys %{$idMap}){
	for(my $i = 0; $i < $idMap->{$_}->[1]; $i += $histInterval){
		#print $_," => $idMap->{$_}->[1]\t",$i,"\t",$i+$histInterval-1,"\t",0,"\n";
		if($i+$histInterval-1 > $idMap->{$_}->[1]){
			push(@{$histogram{$_}}, [$_, $i, $idMap->{$_}->[1], @histFreq]);
		}
		else{
			push(@{$histogram{$_}}, [$_, $i, $i+$histInterval-1, @histFreq]);
		}		
	}
}


#foreach(keys %histogram){
#	print $_,"\n";
#	foreach my $bin(@{$histogram{$_}}){
#		print join("\t", @{$bin},"\n");
#	}
#}

#run BLAST against --refFile
print "Starting: Running BLAST against reference file\n";
#my $blastOut = "../workDir/genome/yeast_sequence.blastout";		#for test direct blastout file
my $blastOut = &blast($blastPath, $refFile, $queryFile, $proc);
print "Finished: Running BLAST\n";

#create Histogram and Links file
print "Starting: Processing of BLAST output $blastOut to create LINKS and HISTOGRAM\n";
&histogramAndLinks($blastOut, $dataStorePath, $idMap, \%histogram);
print "Finished: BLAST output $blastOut processing\n";


#update circos.conf file with the appropriate variables and copy it to working_directory/circos/etc folder
print "Starting: Updating circos.conf file using template files\n";
my $circosConfFile = &updateDenovoCircosConf("$confTemplate/circosDenovo.template", $confStorePath, \%config);
print "Finished: Updating circos.conf file\n";


#update ticks.conf file with the appropriate variables and copy it to working_directory/circos/etc folder
print "Starting: Updating ticks.conf\n";
&updateTicksConf("$confTemplate/ticks.template", $confStorePath, \%config);
print "Finished: Updating ticks.conf\n";


#copy other configs to the circos working directory
&copyOtherConfs($confTemplate, $confStorePath);


#run circos
print "Starting: Running circos to generate plot\n";
Circos($circosPath, $circosConfFile, $kmer, $imageMapFile, $confStorePath, $dashboardDir);
print "Finished: Running circos to generate plot\n";
