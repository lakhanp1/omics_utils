#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

use DataFilePrepare qw(karyotypeData setChrUnit %circosFiles);
use CnvDataFilePrepare qw(cnvGlyphFile);
use UpdateConfigs qw(updateCnvCircosConf copyOtherConfs updateTicksConf);
use RunTools qw(Circos);

my $usage = '
This script parse the CNVnator output and generate the data files required for circos plot. 
It creates the circos templates on the fly for current data and runs the circos tool.
USAGE: $SANGENIX_HOME_PERL/bin/perl cnvCircos.pl --cnv <CNVnator_output> --refFile <REFERENCE_FILE_FOR_KARYOTYPE> --circosPath <CIRCOS_EXE_FILE_COMPLETE_PATH> --confTemplate <DIR_PATH_FOR_CIRCOS_CONF_TEMPLATES> --perl $SANGENIX_HOME_PERL';



#perl -I/home/Lakhan/5_cnvCircos/code/ code/cnvCircos.pl  --cnv otherData/CNV.bed --refFile otherData/genome.fa --circosPath /home/Lakhan/circos-0.67/bin/circos --confTemplate conf/


my $workflowPath = ".";
my $cnv = "";
my $refFile = "";
my $circosPath = "";
my $confTemplate = "";

GetOptions(
	"cnv=s" => \$cnv,
	"refFile=s" => \$refFile,
	"circosPath=s" => \$circosPath,
	"confTemplate=s" => \$confTemplate)
or die "Error in command line arguments: $usage";


if(!-e $cnv){
	die "CNVnator output file $cnv does not exists: $!";
}
if(!-e $refFile){
	die "Reference sequence file $refFile does not exists: $!";
}
if(!-e $circosPath){
	die "Circos exe file $circosPath does not exists: $!";	
}
if(!-d $confTemplate){
	die "Circos config template dir $confTemplate: $!";
}

print "CNVnator output file $cnv\n";
print "Reference sequence file $refFile\n";


my $dataStorePath = $workflowPath."/circos/data";
my $confStorePath = $workflowPath."/circos/etc";

my $dashboardDir = "$workflowPath/dashboard_output";
my $imageMapFile = "$dashboardDir/cnvCircos_output.result";

if(!-d $dashboardDir){
	mkdir($dashboardDir, 0777) or die "Cannot create directory $dashboardDir: $!";
}

#circos runtime parameters
my %config;
#1) The chromosomes_unit value is used as a unit (suffix "u") to shorten values in other parts of the configuration file.
$config{chromosomes_units} = 100000;

#2) the tick label is derived by multiplying the tick position by 'multiplier' and casting it in 'format'
$config{multiplier} = "1e-5";

#create karyotype file
my $idMap = &karyotypeData($refFile, $dataStorePath);
print "Generated Karyotype file: $circosFiles{KARYOTYPE}\n";

#decide the 'chromosomes_unit', 'multiplier' parameters based on genome size so that the ticks and tick lables will be appropriate
($config{chromosomes_units}, $config{multiplier}) = &setChrUnit();


foreach(keys %config){
	print "$_\t$config{$_}\n";
}


#process the CNVnator output and prepare the datafile to plot text track with glyph
($config{glyphFile}, $config{shortestCnv}, $config{longestCnv}) = &cnvGlyphFile($cnv, $dataStorePath, $idMap);

print "Generated Glyph data file: $config{glyphFile}\n";
print "Length of shortest CNV: $config{shortestCnv}\n";
print "Length of longest CNV: $config{longestCnv}\n";


#update circos.conf file with the appropriate variables and copy it to working_directory/circos/etc folder
print "Starting: Updating circos.conf file using template files\n";
my $circosConfFile = &updateCnvCircosConf("$confTemplate/circosCnv.template", $confStorePath, \%config);
print "Finished: Updating circos.conf file\n";

#update ticks.conf file with the appropriate variables and copy it to working_directory/circos/etc folder
print "Starting: Updating ticks.conf\n";
&updateTicksConf("$confTemplate/ticks.template", $confStorePath, \%config);
print "Finished: Updating ticks.conf\n";

#copy other configs to the circos working directory
&copyOtherConfs($confTemplate, $confStorePath);


#run circos and copy image to dashboard_output dir
print "Starting: Running circos to generate plot\n";
Circos($circosPath, $circosConfFile, -1, $imageMapFile, $confStorePath, $dashboardDir);
print "Finished: Running circos to generate plot\n";







