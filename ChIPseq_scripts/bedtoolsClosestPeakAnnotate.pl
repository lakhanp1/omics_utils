#!/usr/bin/env perl -w


use strict;
use warnings;
use Getopt::Long;
use Bio::DB::HTS::Tabix;
use List::Util qw(max min);

my %options;
my $usage = 'USAGE: 
This script parse the "bedtools closest <peaks.bed> <genes_CDS.bed>" output stream. 
For each peak, a nearest target gene/s is/are found. 
Each peak is categorized into one of more of the following categories: 
upstream, overlapStart, overlapEnd, inside, includeFeature, insideOverlapStart, 
insideOverlapEnd, insideOverlapStartOverlapEnd
Peaks targets are also marked pseudo depending on various conditions

IMPORTANT: "bedtools closest" command should be run with following parameters:
bedtools closest -nonamecheck -k 5 -t all -id -D b -a <peaks.broadPeak> -b <genes/transcripts.bed>


perl bedtoolsClosestPeakAnnotate.pl --peakFormat <peakType> --bedFile <tabix index .bed.gz file>

Arguments:
--peakFormat      <STR>       One of narrowPeak, broadPeak, bed
                              Format of the input peak file
                  
--bedFile         <FILE>      Tabix indexed and bgzipped bed file with standard 6 BED columns for
                              gene or transcripts. This file should be same file used in "bedtools
                              closest" command before pipe
				  			
--bindingInGene   <FLAG>      Whether the TF of interest has binding over gene body. This argument 
                              can be set for peaks called on polII ChIP data. Default: FALSE
                  
--bidirectionCut  <INT>       Distance cutoff to decide a peak as bidirectional for two target genes
                              Both the target genes should be withing these bp distance from peak.
                              Default: 500

--tssRegion       <INT INT>   TSS region for deciding confident peak. This argument should include
                              two INT numbers. First number will be used to decide the upstream
                              region from TSS and second number will be used to decide downstream
                              region from TSS. Eg: -200===TSS===100. Default: 200 100
							  
--help                        Print the USAGE
						  
';

$options{bindingInGene} = 0;

GetOptions(\%options, 'bedFile=s', 'peakFormat=s', 'bindingInGene!', 'bidirectionCut=i', 'tssRegion=i{2}', 'help|h') or die("Error in command line arguments\n");

if($options{'help'} || $options{'h'}){
	print STDERR $usage;
	exit 1;
}


## various defaults
$options{'peakFormat'} = lc($options{'peakFormat'});
my $biDirectionCutoff = 500;			## peak distance cutoff to decide a peak as bidirectional for a target
my $bidirectionSkewFraction = 0.2;
my $minTSS_gapForPseudo = 500;			## minimum distance required between two targets' TSS sites to mark one site as pseudo
my $confidentUpstreamCutoff = 200;		## any target gene with distance between TSS and peak < 200bp is confident target
my $insideCutoff = 100;					## distance cutoff to decide that peak is inside the feature

## any peak which is marked as pseudo_upstream for a gene should be within 500bp from TSS.
## Upstream peak with distance from target gene TSS larger than this is definitely false target and no need to consider
my $ignorePseudoUpbeyond = 500;



if(!$options{'bedFile'} || !-e $options{'bedFile'}){
	print STDERR "Error: Please provide the gene CDS file in BED format\n";
	die;
}

if(!$options{'peakFormat'}){
	print STDERR "Error: Please provide the peak format: narrowPeak|broadPeak|bed\n";
	die;
}

if($options{'bidirectionCut'}){
	$biDirectionCutoff = $options{'bidirectionCut'};
}

if($options{'tssRegion'}){
	$confidentUpstreamCutoff = $options{'bidirectionCut'}->[0];
	$insideCutoff = $options{'bidirectionCut'}->[1];
}


print STDERR '# Bidirectional peak cutoff: ', $biDirectionCutoff, "\n";
print STDERR '# BindingInGene: ', $options{bindingInGene}, "\n";
print STDERR '# distance cutoff to decide that peak is inside the feature: ', $insideCutoff, "\n";
print STDERR '# Peak file format: ', $options{'peakFormat'}, "\n";




## field names
my @narrowPeakCol = qw(peakChr peakStart peakEnd peakId peakLength peakStrand peakEnrichment peakPval peakQval peakSummit);

my @colNames = ();
my @outFields = ();

if($options{'peakFormat'} eq "broadpeak"){
	@colNames = @narrowPeakCol[0..8];
	@outFields = @narrowPeakCol[0..3,6..9];
}
elsif($options{'peakFormat'} eq "bed"){
	@colNames = @narrowPeakCol[0..5];
	@outFields = @narrowPeakCol[0..3];
}
elsif($options{'peakFormat'} eq "narrowpeak"){
	@colNames = @narrowPeakCol;
	@outFields = @narrowPeakCol[0..3,6..9];
}
else{
	print STDERR "Error: Please provide the peak format: narrowPeak|broadPeak|bed\n";
	die;
}


push(@colNames, qw(gChr gStart gEnd gName gScore gStrand peakDist));
push(@outFields, qw(gName peakDist summitDist type bidirectional featureCovFrac relativeSummitPos));

my $tabix = Bio::DB::HTS::Tabix->new( filename => $options{'bedFile'} );

print join("\t", @outFields),"\n";


## peak annotation
my @closestSet = ();
my $peakId = '';

while(my $line = <STDIN>){
	chomp $line;
	my @tmp = split(/\t/, $line);
	
	if($#tmp != $#colNames){
		print STDERR "Expected number of fields: ", $#colNames+1,"\n";
		print STDERR "Number of columns expected not found at line:\n", $line,"\n";
		die;
	}
	
	
	if($peakId ne $tmp[3]){
		# print "$peakId != $tmp[3]\n";
		
		## dont process for first line 
		if($#closestSet != -1){
			## collection of all lines for peak_i done; now process it
			getNearest(@closestSet);
			
		}
		
		@closestSet = ();
		$peakId = $tmp[3]
	}
	

	my %data = ();
	@data{@colNames} = @tmp;
	
	## decide the peak summit/center position based on the type of peak file
	my $peakCenter = int(($data{peakEnd} - $data{peakStart})/2);
	
	
	if($options{'peakFormat'} eq "narrowpeak"){
		$data{peakSummit} = $data{peakStart} + $data{peakSummit};
	}
	else{
		## other file types where summit position is not available: BED, broadPeak etc
		$data{peakSummit} = $data{peakStart} + $peakCenter;
	}
	
	
	
	push(@closestSet, \%data);
	
}

## process the last peak
getNearest(@closestSet);


$tabix->close;

























sub getNearest{
	# print "Set:\n";
	if($#_ == -1){
		print STDERR "Error: Target gene set cannot be empty\n";
		exit 1;
	}
	
	my @upstreamTssTarget = ();
	my @tssTarget = ();
	my $tesTarget = undef;
	my @featureInPeakTarget = ();
	my $peakInFeatureTarget = undef;
	my $bidirectionalTarget = undef;
	my @otherTargets = ();
	
	
	my $foundPeakUpTss = 0;
	my $foundPeakNearTSS = 0;
	my $foundPeakNearTES = 0;
	my $foundFeatureInPeak = 0;
	my $foundPeakInFeature = 0;
	my $foundBidirectional = 0;
	
	my $searchBidirectional = 0;
	my $markPseudoTarget = 0;
	my $minDist = undef;
	
	foreach my $target(@_){
		
		my $featureCovFrac = 0.00;
		## upstream feature
		if($target->{peakDist} < 0){
			$target->{type} = "upstream";
		}
		elsif($target->{peakDist} == 0){
			$featureCovFrac = sprintf("%.2f", (min($target->{peakEnd}, $target->{gEnd}) - max($target->{peakStart}, $target->{gStart})) / ($target->{gEnd} - $target->{gStart}));
		}
		
		
		## decide the relative position of peak summit w.r.t. peak to decide the peak tail direction
		## smaller the $summitRelativePos, longer the tail of peak towards gene body
		## IMP: for the target on -ve strand, this value should be 1 - $summitRelativePos
		my $summitRelativePos = sprintf("%.2f", ($target->{peakSummit} - $target->{peakStart}) / ($target->{peakEnd} - $target->{peakStart}));
		
		if($summitRelativePos < 0 || $summitRelativePos > 1){
			print STDERR "#**",join("\t", @{$target}{@outFields}),"*\n";
			print STDERR "Error: Failed while calculating summit tail direction\n";
			print STDERR "Error: Relative summit position should be: 0 < summitRelativePos < 1\n";
			print STDERR "Error: Calculated relative summit position = $summitRelativePos\n";
			die;
		}

		
		$target->{featureCovFrac} = $featureCovFrac;
		$target->{bidirectional} = 0;
		
		## assign the target type
		## strand: +ve
		if($target->{gStrand} eq '+'){
			## distance from summit
			$target->{summitDist} = $target->{peakSummit} - $target->{gStart};
			
			## overlapping feature
			if($target->{peakDist} == 0){
				##
				##          |>=====>=====>======>======>======>|
				##        -------
				##                                         -------
				##                      -------
				##       -------------------------------------------
				##
				
				## start overlap 
				if($target->{peakStart} < $target->{gStart} && $target->{peakEnd} < $target->{gEnd}){
					$target->{type} = "overlapStart";
				}
				## end overlap
				elsif($target->{peakStart} > $target->{gStart} && $target->{peakEnd} > $target->{gEnd}){
					$target->{type} = "overlapEnd";
				}
				## inside gene body
				elsif($target->{peakStart} >= $target->{gStart} && $target->{peakEnd} <= $target->{gEnd}){
					$target->{type} = "inside";
					
					## inside + very near to START: insideOverlapStart
					if(($target->{peakStart} - $target->{gStart}) <= $insideCutoff){
						$target->{type} .= "OverlapStart";
					}
					## inside + very near to END: insideOverlapEnd
					elsif(($target->{gEnd} - $target->{peakEnd}) <= $insideCutoff){
						$target->{type} .= "OverlapEnd";
					}
					
				}
				## include feature
				elsif($target->{peakStart} <= $target->{gStart} && $target->{peakEnd} >= $target->{gEnd}){
					$target->{type} = "includeFeature";
				}
			}
		}
		## strand: -ve
		elsif($target->{gStrand} eq '-'){
			## distance from summit
			$target->{summitDist} = $target->{gEnd} - $target->{peakSummit};
			
			## overlapping feature
			if($target->{peakDist} == 0){
				##
				##          |<=====<=====<======<======<======<|
				##        --------
				##                                         --------
				##                      -------
				##       --------------------------------------------
				##	[2,627,052	2,627,730]		[2,627,458	2,627,730]
				
				## start overlap 
				if($target->{peakStart} > $target->{gStart} && $target->{peakEnd} > $target->{gEnd}){
					$target->{type} = "overlapStart";
				}
				## end overlap 
				elsif($target->{peakStart} < $target->{gStart} && $target->{peakEnd} < $target->{gEnd}){
					$target->{type} = "overlapEnd";
				}
				## inside gene body	
				elsif($target->{peakStart} >= $target->{gStart} && $target->{peakEnd} <= $target->{gEnd}){
					$target->{type} = "inside";
					
					## inside + very near to START: insideOverlapStart
					if(($target->{gEnd} - $target->{peakEnd}) <= $insideCutoff){
						$target->{type} .= "OverlapStart";
					}
					## inside + very near to END: insideOverlapEnd
					elsif(($target->{peakStart} - $target->{gStart}) <= $insideCutoff){
						$target->{type} .= "OverlapEnd";
					}
				}
				## include feature 
				elsif($target->{peakStart} <= $target->{gStart} && $target->{peakEnd} >= $target->{gEnd}){
					$target->{type} = "includeFeature";
				}
			}
		}
		
		
		if(!exists($target->{type})){
			print STDERR "ERROR: Could not assign target type to peak-target pair:\n", join("\t", @{$target}{@colNames[0..3,6..13,15], qw(peakDist summitDist)}),"\n";
			die;
		}
		
		
		## to handle peaks for factor like cclA. cclA has peaks all over the gene body. 
		if($featureCovFrac >= 0.7){
			if($options{bindingInGene} && ($target->{type} eq "overlapStart" || $target->{type} eq "overlapEnd" || $target->{type} =~m/^inside/)){
				$target->{type} = "includeFeature";
			}
			elsif(!$options{bindingInGene} && ($target->{type} eq "overlapStart" || $target->{type} =~m/^inside/ || $target->{type} eq "overlapEnd")){
				$target->{type} = "includeFeature";
			}
			
		}
		
		## for each target type, decide whether to finalize the target or not
		## for overlapping peaks
		if($target->{type} eq "overlapStart" || $target->{type} eq "insideOverlapStart"){
			## overlap near TSS
			push(@tssTarget, $target);
			$foundPeakNearTSS = 1;
			
			## search for bidirectional peaks only if the current found target is of type "overlapStart"
			if($target->{type} eq "overlapStart"){
				$searchBidirectional = 1;
			}
			
			## if a peak near TES was found before this, set it to undef
			if($foundPeakNearTES){
				$foundPeakNearTES = 0;
				$tesTarget = undef;
			}
		}
		## for gene in peak
		elsif($target->{type} eq "includeFeature"){
			push(@featureInPeakTarget, $target);
			$foundFeatureInPeak = 1;

			if($options{bindingInGene}){
				## This is best peak so dont search upstream peaks
				$minDist = $target->{peakDist};
			}
		}
		## for overlap near TES
		elsif($target->{type} eq "overlapEnd" || $target->{type} eq "insideOverlapEnd"){
			## update summit distance with respect to the CDS end
			if($target->{gStrand} eq '+'){
				$target->{summitDist} = $target->{peakSummit} - $target->{gEnd};
			}
			elsif($target->{gStrand} eq '-'){
				$target->{summitDist} = $target->{gStart} - $target->{peakSummit};
			}
			
			## if previous peak at TSS has not been found, then only consider the peak near TES
			if(!$foundPeakNearTSS){
				$tesTarget = $target;
				$foundPeakNearTES = 1;
				
				## exception for the factors like cclA which have binding over gene body. 
				## for genes with known binding in gene body, these set of peaks are true peaks.
				if($options{bindingInGene}){
					## This is good peak even though it is at TES. So upstream peaks wont be searched
					$minDist = $target->{peakDist};
				}
				# # elsif(!$options{bindingInGene}){
					# # # search bidirectional peaks only for the factor which is not expected to bind in gene body
					# # $minDist = $biDirectionCutoff;
				# # }
			}
			
		}
		## for peak in gene
		elsif($target->{type} eq "inside"){ 
			$peakInFeatureTarget = $target;
			$foundPeakInFeature = 1;
			
			## exception for the factors like cclA which have binding over gene body
			## for genes with known binding in gene body, these set of peaks are true peaks.
			## hence search bidirectional peaks only for the factor which does not bind in gene body
			if($options{bindingInGene}){
				## This is good peak even though it is at TES. So upstream peaks wont be searched
				$minDist = $target->{peakDist};
			}
			elsif(!$options{bindingInGene}){
				# other upstream targets will be searched with distance less than this
				$minDist = $biDirectionCutoff;
			}

		}
		elsif($target->{type} eq "insideOverlapStartOverlapEnd"){
			push(@tssTarget, $target);
		}
	
		
		## for the peaks which are inside the gene body, $summitRelativePos w.r.t. peak does not make sense
		## instead, report the relative summit position w.r.t. gene
		if($target->{type} =~m/^inside/){
			$summitRelativePos = sprintf("%.2f", ($target->{peakSummit} - $target->{gStart}) / ($target->{gEnd} - $target->{gStart}));		
		}
		
		
		## if the target gene is on -ve strand, summitRelativePos = 1 - summitRelativePos
		if($target->{gStrand} eq '-'){
			$summitRelativePos = 1 - $summitRelativePos;
		}
		
		$target->{relativeSummitPos} = $summitRelativePos;			
		
		
		## for upstream peaks,
		## if first upstream peak i.e. foundPeakUpTss == 0: add this target to @upstreamTssTarget, set searchBidirectional = 1 for next upstream peaks
		## if foundPeakUpTss == 1 OR searchBidirectional == 1: check if this peak is upstream of other gene too
		if($target->{type} eq "upstream"){
			
			if(!defined($minDist) || abs($minDist) > abs($target->{summitDist})){
				
				## upstream peaks with dist <= $biDirectionCutoff: confident targets 
				if(abs($target->{peakDist}) <= $biDirectionCutoff){
					
					## mark this upstream peak as bidirectional if searchBidirectional = TRUE
					if($searchBidirectional){
						
						if($foundPeakNearTSS && $foundPeakUpTss){
							$target->{bidirectional} = 3;
						}
						elsif($foundPeakNearTSS){
							## for one target, peak is overlapStart|insideOverlapStart and for other, the peak is upstream
							$target->{bidirectional} = 1;
						}
						elsif($foundPeakUpTss){
							## for both the targets, peak is upstream
							$target->{bidirectional} = 2;
						}
						
						## $bidirectionalTarget = $target;
						## need geneBetweenPeakAndTarget($target) filter here to handle the following situation
						##  ------
						##        --------   
						##            ====>===>==  ===>===>===
						##  peak2  peak1   gene1      gene2
						## if gene2 is within $biDirectionCutoff, it will come as upstream and will get classified as bidirectional
						## this can be true for peak1 or peak2
						
						my $noGeneInbetween = &geneBetweenPeakAndTarget($target);
						if($noGeneInbetween){
							push(@upstreamTssTarget, $target);
							$foundBidirectional = 1;
						}
						
					}
					## upstream peaks which are at distance < $confidentUpstreamCutoff (200bp) from TSS of gene1 and include another gene2 withing peak
					elsif($foundFeatureInPeak && !$options{bindingInGene} && (abs($target->{peakDist}) <= $confidentUpstreamCutoff)){
						##  ===<===<===                              ===<===<===
						##                   ===<===<===                             ===<===<===
						##              -------------------                       -------------------
						##                 *                                                     *
						##   gene1       peak   gene2                     gene1     peak    gene2
						## >>>>>>>> Example1 <<<<<<<<<<<<<<<<<<       >>>>>>>> Example2 <<<<<<<<<<<<<<<<<<  
						## * : summit
						## can't resolve this based on summit because most of the times, peaks in such cases are multimodal
						## include both, gene1 and gene2 as targets for peak
						
						push(@upstreamTssTarget, $target);
						$foundPeakUpTss = 1;
						# print STDERR "#**#",join("\t", @{$featureInPeakTarget[0]}{@outFields}),"*\n";
						# print STDERR "#**",join("\t", @{$target}{@outFields}),"*\n";
					}
					## normal upstream peak which is closer to the target
					elsif(!$foundPeakNearTSS && !$foundFeatureInPeak){
						push(@upstreamTssTarget, $target);
						$foundPeakUpTss = 1;
						
						## this is the first upstream peak found. So for next upstream peaks, set searchBidirectional = 1
						$searchBidirectional = 1;
					}
				}
				## other upstream peaks which are farther: check if there is any other gene inbetween
				else{
					## check if there is another gene between peak and current target
					my $noGeneInbetween = &geneBetweenPeakAndTarget($target);
					
					if($noGeneInbetween){
						$minDist = $target->{peakDist};
					}
					
					# print STDERR "##**",join("\t", @{$target}{@outFields}),"*$noGeneInbetween\n";
					
					## final filter: there should not be any other gene in between the peak and its target					
					if($noGeneInbetween && !$foundPeakNearTSS && !$foundFeatureInPeak){
						
						## if already foundPeakUpTss: this is second Upstream target and has to be within $ignorePseudoUpbeyond distance from peak
						if($foundPeakUpTss && abs($target->{peakDist}) < $ignorePseudoUpbeyond){
							push(@upstreamTssTarget, $target);
						}
						elsif(!$foundPeakUpTss){
							push(@upstreamTssTarget, $target);
							$foundPeakUpTss = 1;
						}
						# $searchBidirectional = 1;
					}
				}
			}
		}
		
		# print STDERR "#**",join("\t", @{$target}{@outFields}),"*\n";
	}
	

	#@upstreamTssTarget       $foundPeakUpTss          upstream
	#@tssTarget               $foundPeakNearTSS        overlapStart, insideOverlapStart, insideOverlapStartOverlapEnd
	#$tesTarget               $foundPeakNearTES        overlapEnd, insideOverlapEnd
	#$featureInPeakTarget     $foundFeatureInPeak      includeFeature
	#$peakInFeatureTarget     $foundPeakInFeature      inside
	#$bidirectionalTarget     $foundBidirectional      
	
	# print STDERR "# $peakId targets:\t";
	# print STDERR join("\t", "foundBidirectional=$foundBidirectional", "foundPeakUpTss=$foundPeakUpTss", "foundPeakNearTSS=$foundPeakNearTSS", "foundPeakNearTES=$foundPeakNearTES", "foundFeatureInPeak=$foundFeatureInPeak", "foundPeakInFeature=$foundPeakInFeature"),"\n";

	
	## when the peak has two targets and for both the targets, peak is upstream
	## this is done even for bidirectional peaks
	if($#upstreamTssTarget == 1){
		@upstreamTssTarget = &nearestFromBidirectional(@upstreamTssTarget);
	}
	
	## when the peak has two targets and it is bidirectional:
	## decide the appropriate target based on the closeness of the peak with target
	if($foundBidirectional){
		## one target gene has peak at overlapStart|insideOverlapStart and for other target peak is at upstream 
		if($foundPeakNearTSS){			
			my $distBetween = max($upstreamTssTarget[0]->{gStart}, $tssTarget[0]->{gStart}) - min($upstreamTssTarget[0]->{gEnd}, $tssTarget[0]->{gEnd});
			
			## if the distance between TSS of two genes < 500 : can't decide the other target as pseudo
			## if upstream peak is within $confidentUpstreamCutoff (300bp) of target: do not set it to pseudo
			## if the distance between TSS of two genes > 500 : mark other upstream target as pseudo
			if($distBetween > $minTSS_gapForPseudo && abs($upstreamTssTarget[0]->{peakDist}) > $confidentUpstreamCutoff){
				@upstreamTssTarget = map{&setTargetToPesudo($_)} @upstreamTssTarget;
			}
		}
	}
	
	
	
	## mark pseudo targets: 
	## (overlapEnd, insideOverlapEnd) = pseudo if (upstream|overlapStart|insideOverlapStart|insideOverlapStartOverlapEnd|includeFeature) 
	if($foundPeakNearTSS){
		## mark the tesTarget targets as pseudo targets
		$tesTarget = &setTargetToPesudo($tesTarget);
	}
	
	## remove tesTarget if a foundFeatureInPeak and  a peak near TES is found
	if($foundFeatureInPeak){		
		$foundPeakNearTES = 0;
		$tesTarget = undef;
	}
	
	
	## for the TF which has known binding over gene body. Can also be used for peaks called on polII data
	if($options{bindingInGene}){
		## preference is for "includeFeature". all other targets are pseudo
		if($foundFeatureInPeak){
			## $tesTarget = &setTargetToPesudo($tesTarget);
			# @tssTarget = map{&setTargetToPesudo($_)} @tssTarget;
			@upstreamTssTarget = map{&setTargetToPesudo($_)} @upstreamTssTarget;
			## $bidirectionalTarget = &setTargetToPesudo($bidirectionalTarget);
		}
	}
	## for the TFs with normal binding near promoter region
	elsif(!$options{bindingInGene}){
		## for normal TFs, preference is for binding near TSS
		if($foundPeakNearTSS){
			## mark the "inside" targets as pseudo targets
			$peakInFeatureTarget = &setTargetToPesudo($peakInFeatureTarget);	
		}
		
		if($foundPeakUpTss && abs($upstreamTssTarget[0]->{peakDist}) <= $biDirectionCutoff){
			$tesTarget = &setTargetToPesudo($tesTarget);
			$peakInFeatureTarget = &setTargetToPesudo($peakInFeatureTarget);
		}
	
	}

	

	
	if($foundFeatureInPeak){
		push(@otherTargets, @featureInPeakTarget);
	}
	if($foundPeakNearTES){
		push(@otherTargets, $tesTarget);
	}
	if($foundPeakInFeature){
		push(@otherTargets, $peakInFeatureTarget);
	}
	if($foundPeakUpTss || $foundBidirectional){
		push(@otherTargets, @upstreamTssTarget);
	}
	# if($foundBidirectional){
		# push(@otherTargets, $bidirectionalTarget);
	# }
	
	
	foreach my $t(@tssTarget, @otherTargets){
		print join("\t", @{$t}{@outFields}),"\n";
	}
	
	## if nothing found, print the NA values
	if(!$foundPeakUpTss && !$foundPeakNearTSS && !$foundPeakNearTES && !$foundFeatureInPeak && !$foundPeakInFeature){
		my $t = $_[0];
		@{$t}{qw(gName peakDist summitDist type bidirectional featureCovFrac)} = qw(NA NA NA NA NA NA);
		print join("\t", @{$t}{@outFields}),"\n";
	}
	# print "\n";
	
}






















## for bidirectional targets, select the target which is closest to the peak based on conditions:
## 90% of the peak lies on the side of the target at midpoint of two targets

sub nearestFromBidirectional{
	my ($t1, $t2) = @_;
	
	##     target1                 peak                                 target2
	## ==<=====<=====<===         ----------                     ===>=====>=====>=====>== 
	##                                      |
	##                         center between two targets
	## True target: target1

	
	my ($from, $to) = (undef, undef);
	my ($negative, $positive) = (undef, undef);
	
	if($t1->{gStrand} ne $t2->{gStrand}){
	
		if($t1->{gStrand} eq '+' && $t2->{gStrand} eq '-'){
			$positive = $t1;
			$negative = $t2;
			
			$from = $negative->{gEnd};
			$to = $positive->{gStart};
		}
		elsif($t1->{gStrand} eq '-' && $t2->{gStrand} eq '+'){
			$positive = $t2;
			$negative = $t1;
			
			$from = $negative->{gEnd};
			$to = $positive->{gStart};
		}
		else{
			print STDERR "ERROR: Unexpected strand for bidirectional targets:\n";
			print STDERR 'target1: ', join("\t", @{$t1}{@colNames[0..3], qw(gChr gStart gEnd gName gStrand peakDist type)}),"\n";
			print STDERR 'target2: ', join("\t", @{$t2}{@colNames[0..3], qw(gChr gStart gEnd gName gStrand peakDist type)}),"\n";
			die;
		}
		
		my $center = $from + int(($to - $from)/2);
		
		my $peakLen = $t1->{peakEnd} - $t1->{peakStart};
		my $fraction = $peakLen * $bidirectionSkewFraction;
		my $distBetweenTragets = abs($to - $from);
		
		# print STDERR "from = $from | to = $to | center = $center | 20% = ", ($t1->{peakStart} + $fraction) , " | 80% = ",($t1->{peakEnd} - $fraction),"\n";
		
		## if the distance between the TSS of two bidirectional targets < $minTSS_gapForPseudo: CANNOT decide the pseudo target confidently
		if($distBetweenTragets > $minTSS_gapForPseudo){
			## negative strand target is true; set positive strand target to pseudo
			if(($t1->{peakEnd} - $fraction) <= $center){
				$positive = &setTargetToPesudo($positive);
				# print STDERR "*farthestFromBidirectional*\t",join("\t", @{$positive}{@outFields}),"*\n";
			}
			## positive strand target is true; set negative strand target to pseudo
			elsif(($t1->{peakStart} + $fraction) >= $center){
				$negative = &setTargetToPesudo($negative);
				# print STDERR "*farthestFromBidirectional*\t",join("\t", @{$negative}{@outFields}),"*\n";
			}
		}
		
		return($negative, $positive);
	}
	else{
		return @_;
	}
	
}





sub setTargetToPesudo{
	my $t = shift @_;
	if(defined $t){
		# print "marking pseudo: $t\n";
		$t->{type} = 'pseudo_'.$t->{type};
	}
	
	return($t);
}




sub geneBetweenPeakAndTarget(){
	my $target = shift @_;
	
	## handle such instances:
	##       peak             middle gene             target
	##       --------      |<====<====<====<|     |>====>====>====>|
	##
	##       target                   middle gene         peak
	##     |<====<====<====<|     |>====>====>====>|	 --------
	##
	if($target->{type} ne "upstream"){
		return 1;
	}
	else{
		my $region = $target->{peakChr};
		my ($from, $to) = (undef, undef);
		
		if($target->{gStrand} eq '+'){
			$from = $target->{peakEnd};
			$to = $target->{gStart};
		}
		elsif($target->{gStrand} eq '-'){
			$from = $target->{gEnd};
			$to = $target->{peakStart};
		}
		
		$region .= ':'.$from.'-'.$to;
		# print STDERR "##Region: $region\n";
		
		my $iter = $tabix->query($region);
		
		## check if there is any gene between the region [peak's nearest boundry, target nearest boundry]
		while ( my $n = $iter->next ) {
			my @cdsBed = split(/\t/, $n);
			my $bedLen = $cdsBed[2] - $cdsBed[1];
			
			my $overlap = min($to, $cdsBed[2]) - max($from, $cdsBed[1]);
			my $fractionOverp = $overlap/$bedLen;
			# print STDERR "##$n | $target->{gName} | bedLen = $bedLen | Overlap = $overlap | FractionOverlap = $fractionOverp | min(to=$to, cdsEnd=$cdsBed[2]) - max(from=$from, cdsStart=$cdsBed[1])\n";
			
			if($target->{gName} ne $cdsBed[3] && $fractionOverp >= 0.3){
				# print STDERR '######Gene between peak ', $target->{peakId}, ' and gene ', $target->{gName}, ': ', $cdsBed[3],"\tFraction: $fractionOverp\n";
				return 0;
			}
		}
	}
	
	return 1;
}



