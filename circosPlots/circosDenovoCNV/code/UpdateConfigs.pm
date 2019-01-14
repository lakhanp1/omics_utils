package UpdateConfigs;


use strict;
use warnings;

use Template;
use DataFilePrepare qw(%circosFiles);
use File::Path;
use File::Copy;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw(updateDenovoCircosConf updateTicksConf copyOtherConfs updateCnvCircosConf);
our @EXPORT = qw();
our $VERSION = 1.1;

#create circos.conf file for CNVnator dashboard
sub updateCnvCircosConf{
	#parameter to modify from template:
	#1) karyotype
	#2) chromosomes_units
	#3) glyphFile
	#4) shortestCnv
	#5) longestCnv
	my $template = shift;
	my $storeDir = shift;
	my $confHref = shift;

	&createDirPath($storeDir);

	my $confFile = $storeDir."/circos.conf";
	$confHref->{karyotype} = $circosFiles{KARYOTYPE};
	
	#generate the template using TemplateToolkit2
	open(my $outFh, '>:raw', $confFile) or die "Cannot open file $confFile:$!";	
	my $tt2Status = &runTT2($template, $confHref, $outFh);
	
	close($outFh);
	return $confFile;
}


#create circos.conf file from Denovo Assembly dashboard data
sub updateDenovoCircosConf{
	#parameters to be modified
	#1) chromosomes_units
	#2) karyotype 
	#3) histogram files: [file, r0, r1] for each file
	#4) link: file, radius
	my $template = shift;
	my $storeDir = shift;
	my $confHref = shift;

	#craete the path
	&createDirPath($storeDir);
	
	my $confFile = $storeDir."/circos.conf";	
	my $numOfHist = @{$circosFiles{HISTOGRAM}};
	$confHref->{linkFile} = $circosFiles{LINKS};
	$confHref->{karyotype} = $circosFiles{KARYOTYPE};

	my $histTrackDiff = 0.02;
	
	open(my $outFh, '>:raw', $confFile) or die "Cannot open file $confFile:$!";
	
	#print $confHref->{linkRadius},"\n";
	
	#the size on the circos plot which a single histogram track can occupy
	my $histRadiusIncrement = sprintf("%.2f", (0.99 - $confHref->{linkRadius} - ($histTrackDiff * $numOfHist) )/$numOfHist);
	
	my $r0 = $confHref->{linkRadius} + $histTrackDiff;
	my $r1 = $r0 + $histRadiusIncrement;
	
	#print "r0 = $r0\nr1 = $r1\n";
	
	#set the radii r0 and r1 for the histogram plot for differnet %identity cutoffs plots
	for(my $i = 0; $i < @{$circosFiles{HISTOGRAM}}; $i++){
		#print $circosFiles{HISTOGRAM}->[$i],"\n";
		#print "r0 = $r0\nr1 = $r1\n";
		
		$confHref->{hist}->[$i]->{file} = $circosFiles{HISTOGRAM}->[$i];
		$confHref->{hist}->[$i]->{r0} = $r0;
		$confHref->{hist}->[$i]->{r1} = $r1;
		
		$r0 = $r1 + $histTrackDiff;
		$r1 = $r0 + $histRadiusIncrement;
	}
	
	#foreach(@{$confHref->{hist}}){
	#	foreach my $id(sort keys %{$_}){
	#		print "$id\t$_->{$id}\n";
	#	}
	#	print "\n";
	#}
	
	#generate the template using TemplateToolkit2
	my $tt2Status = &runTT2($template, $confHref, $outFh);
	
	#out: a circos.conf file will be created at PATH/circos/etc
	close($outFh);
	return $confFile;
}


#create ticks.conf from ticks.template file
sub updateTicksConf{
	#parameters to be modified
	#1) multiplier
	my $template = shift;
	my $storeDir = shift;
	my $confHref = shift;

	my $confFile = $storeDir."/ticks.conf";
	
	#generate the template using TemplateToolkit2
	open(my $outFh, '>:raw', $confFile) or die "Cannot open file $confFile:$!";
	my $tt2Status = &runTT2($template, $confHref, $outFh);
	#out: a links.conf file will be created at PATH/circos/etc
}


#copy the remaining config files which are not to be modified to the circos/etc
sub copyOtherConfs{
	#Files to be copied
	#1. ideogram.conf
	#2. axes.conf
	#3. backgrounds.conf
	
	my $templatePath = shift;
	my $confPath = shift;
	
	copy("$templatePath/ideogram.conf", "$confPath/ideogram.conf") or die "Cannot copy file $templatePath/ideogram.conf: $!";
	copy("$templatePath/backgrounds.conf", "$confPath/backgrounds.conf") or die "Cannot copy file $templatePath/backgrounds.conf: $!";

	#out: a histogram.conf file will be created at PATH/circos/etc
}

#run TemplateToolkit2 to create a config file from template file
sub runTT2{
	
	my $circosTemplate = Template->new({
		ABSOLUTE => 1,
		TRIM => 1,
	}) or die "Error in TT2: $Template::ERROR";
	
	my $status = $circosTemplate->process(@_) or die $circosTemplate->error();
	
	return $status;
	
}


#craete the directory path
sub createDirPath{
	
	if(!-d $_[0]){
		mkpath($_[0], {mode => 0777, error => \my $error});
		if(@$error){
			for my $diag(@$error){
				my ($folder, $msg) = %$diag;
				if($folder eq ''){
					die "General Error while creating folder: $msg";
				}
				else{
					die "Problem while creating path $folder: $msg";
				}
			}
		}
	}
}



1;