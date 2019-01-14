package DataFilePrepare;


use strict;
use warnings;
use File::Path;
use File::Copy; 


require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(karyotypeData histogramAndLinks setParameters getQueryLength %circosFiles setChrUnit);
our @EXPORT = qw();
our $VERSION = 1.1;


#parameters for plotting Histograms and Links in Circos. These are default parameters. 
our $Hist_Coverage = 51;
our $Hist_Interval_Len = 100000;
our $Hist_CutOffs = [70, 80, 90];
our %Query_Len;
our $Links_Query_Coverage = 80;
our $Link_Threshold = 1000;


#global variables to be exported
our %circosFiles =(
	KARYOTYPE => undef,
	LINKS => undef,
	HISTOGRAM => undef,
	ID_MAP => undef
);





#set the GLOBAL VARIABLES for this package
sub setParameters{
	$Hist_Coverage = shift;
	$Hist_Interval_Len = shift;
	
	$Hist_CutOffs = shift;
	$Hist_CutOffs = [sort{$a <=> $b}@{$Hist_CutOffs}];
	
	$Links_Query_Coverage = shift;
	$Link_Threshold = shift;
	
	#check if the GLOBAL PARAMETER are getting set or not
	&checkParamenter();
}


#check if the paramenter which are defined as GLOBAL VARIABLES in this package are set or not
sub checkParamenter{
	print "Histogram coverage: $Hist_Coverage\n";
	print "Histogram interval length: $Hist_Interval_Len\n";
	print "Different % identity cutoffs for different histograms: ",join("\t", @{$Hist_CutOffs},"\n");
	print "Minimum query coverage for links: $Links_Query_Coverage\n";
	print "Minimum query sequence length to be matched for link: $Link_Threshold\n";
}



#create the karyotype file
sub karyotypeData{
	#usage: &karyotypeData([genome file complete path or karyotype file path], [exec Id folder path]);
	
	my $file = shift;
	my $storeLocation = shift;
	
	#this variables will be returned
	my %seqIdMap;
	my $chromosomeUnit = 10;
	
	#local variables
	
	
	#create the circos/data directory to store the karyotye file
	if(!-d $storeLocation){
		mkpath($storeLocation, {mode => 0777, error => \my $error});
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
	
	
	$circosFiles{KARYOTYPE} = "$storeLocation/karyotype.txt";
	$circosFiles{ID_MAP} = "$storeLocation/idMap.txt";
	
	#read the genome file and prepare the karyotype
	#print "need to create karyotype from genome: $file\n";
	
	#create karyotype file and id_map file to write the karyotype info
	open(my $out, ">", $circosFiles{KARYOTYPE}) or die "Cant create file ",$circosFiles{KARYOTYPE},": $!";
	open(my $idMapOut, ">", $circosFiles{ID_MAP}) or die "Cant create file ",$circosFiles{ID_MAP},": $!";
	
	#read the genome file and extract the co-ordinates
	open(my $fh, "<", $file) or die "Cant open file $file: $!";
	
	#all patterns compiled so code is efficient
	my $fasta = qr/^>(\S+)/;
	my $emptyLine = qr/(^\s+?$)/;
	my $seq = qr/^\w+$/;
	
	
	my ($id, $lable, $start, $end, $color) = ("chr", 0, 0, 0, 0);
	my ($seqId, $oldId) = ("", "");
	while(<$fh>){
		#print $_;
		chomp;
		if($_=~ m/$fasta/){
			
			$seqId = $1;
			
			if($lable == 0){
				$lable++;
				$color++;
				$oldId = $seqId;
				next;
			}
			
			#reset the color to chr1 after chr24
			if($color == 25){
				$color = 1;
			}
			
			
			print $out $id." - ".$oldId." ".$oldId." ".$start." ".$end." ".$id.$color."\n";
			print $idMapOut "$oldId\t$id$color\t$end\n";
			$seqIdMap{$oldId} = [$id.$color, $end];
			
			$lable++;
			$color++;
			$oldId = $seqId;
			$end = 0;
						
		}
		elsif($_ =~ m/$emptyLine/){
			next;
		}
		elsif($_ =~ m/$seq/){
			$end += length($_);
		}
	
	
	}
	
	#print details for the last fasta record
	print $out $id." - ".$oldId." ".$oldId." ".$start." ".$end." ".$id.$color."\n";
	print $idMapOut "$oldId\t$id$color\t$end\n";
	$seqIdMap{$oldId} = [$id.$color, $end];
	#I=>[chr1, 230218]
	#II=>[chr2, 813184]
	
	close($fh);
	close($out);
	close($idMapOut);

	
	#foreach(sort keys %seqIdMap){
	#	print "$_\t@{$seqIdMap{$_}}\n";
	#}
	#out: the karyotype file will be created or copied at data/ and its complete path will be returned
	
	return \%seqIdMap;
}


#set the chromosomes_unit and multiplier parameters for tick lables
sub setChrUnit{
	my $chrUnit = 10;
	my $multiplierPower = -4;
	
	my $total = 0;
	my $lowerLim = 1000;
	my $upperLim = 5000;
	
	open(my $fh, $circosFiles{KARYOTYPE}) or die "Cannot open file $circosFiles{KARYOTYPE}: $!";
	
	while(<$fh>){
		my @temp = split(/\s+/, $_);
		$total += $temp[5];
	}
	
	close($fh);
	
	#decide the chromosomes_unit
	my @probableUnits;
	
	for(my $i=1; $i>0; $i*=10){
		#print "$i\n";
		
		if(($total/$i) > $upperLim){
			next;
		}
		elsif(($total/$i) < $lowerLim){
			last;
		}	
	
		my $tempLess = $i - $i/2;
		my $tempMore = $i*5;
		
		#print "$tempLess\t$i\t$tempMore\n";
		
		my @tempUnits = sort{$a->[1] <=> $b->[1]}map{[$_, abs($lowerLim + ($upperLim - $lowerLim)/2 - (int $total/$_))]}($tempLess, $i, $tempMore);
		push(@probableUnits, $tempUnits[0]);
	}
	
	@probableUnits = sort{$a->[1] <=> $b->[1]}@probableUnits;
	
	$chrUnit = $probableUnits[0]->[0];
	
	#decide the index for the multiplier
	for(my $i=1; $i>0; $i++){
		my $remainder = int($total % (10**$i));
		
		
		#if($remainder >= 10000 && $remainder < 100000){
		if($remainder >= $chrUnit){
			$multiplierPower = $i;
			last;
		}
		
		if($remainder > $total){
			last;
		}
	}
	
	#print "chromosome_units = $chrUnit\nmultiplier = 1e-$multiplierPower\n";
	
	return ($chrUnit, "1e-".$multiplierPower);
}



sub histogramAndLinks{
	
	my $blastOut = shift;
	my $storeLocation = shift;
	my $idRef = shift;
	my $histDataRef = shift;
	#IX => [[chr5	0	24999	0	0	0], [chr5	25000	49999	0	0	0], ....]
	
	#create the circos/data directory to store the histogram.txt and links.txt file
	if(!-d $storeLocation){
		mkpath($storeLocation, 1, 0777);
	}
	
	#create the links.txt file	
	$circosFiles{"LINKS"} = $storeLocation."/links.txt";
	
	print "Creating file to write LINKS data: ",$circosFiles{"LINKS"},"\n";
	open(my $lOut, ">", $circosFiles{"LINKS"}) or die "Cant open file ",$circosFiles{"LINKS"},": $!";
	
	#read BLAST outpt and prepare the histogram and links files
	open(my $in, "<", $blastOut) or die "Cant open file $blastOut: $!";
	
	#query_name | subject_name | percent_identities | aligned_length | #_mismatches | #_gaps | query_start | query_end | subject_start | subject_end | subject_start | e-value | bit_score
	
	my ($tmpQry, $newQry, $fileStart) = ("", 0, 1);
	my @tmpQryHits = ();
	
	my $commentLine = qr/^#/;
	
	my $i=0;
	while(my $line = <$in>){
		if($line =~m/$commentLine/){
			next;
		}
		
		my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore) = split(/\t/, $line);
		
		#calculate query coverage
		my $queryCov = 100 * $alnLength / $Query_Len{$queryId};
		
		#check if this hit to be considered based on the thresholds
		
		
		if($queryId ne $tmpQry){
			#process previous query hits for histogram and links
			if($#tmpQryHits != -1){
				#print "@{$tmpQryHits[0]}\n";
				&linksAndHistProcess(\@tmpQryHits, $histDataRef, $idRef, $lOut);
			}
			
			#reset the hits data for current query
			$tmpQry = $queryId;
			$newQry = 1;
			@tmpQryHits = ();
		}
		
		#
		push(@tmpQryHits, [$queryId, $subjectId, $percIdentity, $alnLength, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $queryCov]);
		
		
		#$i++;
		#if($i == 1200){
		#	last;
		#}
	}
	
	
	#process for the last BLAST output
	&linksAndHistProcess(\@tmpQryHits, $histDataRef, $idRef, $lOut);
	
	close($in);
	close($lOut);
	
	#create the Histogram files  for each of the %identity cutoffs and write the data to the files
	print "Writing histogram data to the histogram files\n";
	&writeHistogramData($storeLocation, $histDataRef);

	
	#out: file path for Histogram data files and links data file
	
}


#process the BLAST hits to get histogram and links data
sub linksAndHistProcess{
	my $hitsRef = shift;
	
	#IX => [[queryId0, subjectId0, percIdentity0, alnLength0, queryStart0, queryEnd0, subjectStart0, subjectEnd0, eVal0, queryCov0], [queryId1, subjectId1, percIdentity1, alnLength1, queryStart1, queryEnd1, subjectStart1, subjectEnd1, eVal1, queryCov1], ....]
	my $rawHistogram = shift;
	
	my $idsRef = shift;
	my $linkFh = shift;

	
	#print "$Hist_Coverage\n";
	
	for(my $i=1; $i<@{$hitsRef}; $i++){
		
		#print "$Query_Len{$hitsRef->[$i]->[0]}\t@{$hitsRef->[$i]}\n";
		#calculate the frequency for the intervals for plotting Histogram
		#this gives the index of the bin for a particular chromosome
		my $lowerLim = int($hitsRef->[$i]->[6] / $Hist_Interval_Len);
		my $upperLim = int($hitsRef->[$i]->[7] / $Hist_Interval_Len);
		

		my @incrementBins = ();
		
		if($lowerLim == $upperLim){
			#BLAST hit lies within a single interval
			
			push(@incrementBins, $lowerLim);
			#print "@{$rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]}\n";
			
		}
		elsif($lowerLim < $upperLim){
			#BLAST hit lies at the interval border so it is in two intervals
			#print join("\t", @{$hitsRef->[$i]}, $lowerLim, $upperLim, "\n");

			#($lowerLim < $upperLim) means a BLAST hit can span over two or more bins. The condition below will handle if a hit spans over more than 2 bins
			if(abs($lowerLim - $upperLim) > 1){
				#consider all the bins which are covered by the blast hit
				push(@incrementBins, $lowerLim+1..$upperLim-1);
				#print "($lowerLim - $upperLim):\t",join(",  ", $lowerLim+1..$upperLim-1),"\n";
			}
			
			#check if the hit can be placed in LOWER interval based on its coverage
			if(($rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]->[2] - $hitsRef->[$i]->[6]) > ($hitsRef->[$i]->[3] * $Hist_Coverage / 100)){
				
				unshift(@incrementBins, $lowerLim);
				#print "@{$rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]}\t",($rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]->[2] - $hitsRef->[$i]->[6]),"\n";
			}
			
			#check if the hit can be placed in UPPER interval based on its coverage
			if(($hitsRef->[$i]->[7] - $rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]->[1]) > ($hitsRef->[$i]->[3] * $Hist_Coverage / 100)){
				
				push(@incrementBins, $upperLim);
				#print "@{$rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]}\t",($hitsRef->[$i]->[7] - $rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]->[1]),"\n";
			}
			
			

			
		}
		elsif($lowerLim > $upperLim){
			#BLAST hit lies at the interval border so it is in two intervals
			#also this hit is on reverse complementary strand. if hit is on negative strand the BLAST output range is reversed
			
			#print join("\t", @{$hitsRef->[$i]}, $lowerLim, $upperLim, "\n");
			
			#($lowerLim < $upperLim) means a BLAST hit can span over just two or more bins. The condition below will handle if a hit spans over more than 2 bins
			if(abs($lowerLim - $upperLim) > 1){
				#consider all the bins which are covered by the blast hit
				push(@incrementBins, $upperLim+1..$lowerLim-1);
				#print "($lowerLim - $upperLim):\t",join(",  ", $upperLim+1..$lowerLim-1),"\n";
			}
			
			#check if the hit can be placed in LOWER interval based on its coverage
			if(($hitsRef->[$i]->[6] - $rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]->[1]) > ($hitsRef->[$i]->[3] * $Hist_Coverage / 100)){
				
				push(@incrementBins, $lowerLim);
				#print "** @{$rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]}\t",($hitsRef->[$i]->[6] - $rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]->[1]),"\n";
			}
			
			
			#check if the hit can be placed in UPPER interval based on its coverage
			if(($rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]->[2] - $hitsRef->[$i]->[7]) > ($hitsRef->[$i]->[3] * $Hist_Coverage / 100)){
				
				unshift(@incrementBins, $upperLim);
				#print "@{$rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]}\t",($rawHistogram->{$hitsRef->[$i]->[1]}->[$upperLim]->[2] - $hitsRef->[$i]->[7]),"\n" ;	
			}
			
		}

		
		#increment the histogram frequency for each of the Cutoffs
		foreach my $bin(@incrementBins){					
			
			for(my $k = 0; $k < @{$Hist_CutOffs}; $k++){
				if($hitsRef->[$i]->[2] >= $Hist_CutOffs->[$k]){
					#IX => [[queryId0, subjectId0, percIdentity0, alnLength0, queryStart0, queryEnd0, subjectStart0, subjectEnd0, eVal0, queryCov0], [queryId1, subjectId1, percIdentity1, alnLength1, queryStart1, queryEnd1, subjectStart1, subjectEnd1, eVal1, queryCov1], ....]
					$rawHistogram->{$hitsRef->[$i]->[1]}->[$bin]->[-@{$Hist_CutOffs} + $k]++;
				}
			}			
		}
		
		#print "@{$rawHistogram->{$hitsRef->[$i]->[1]}->[$lowerLim]}\n";
		
		#apply the query coverage filter
		if($hitsRef->[$i]->[-1] < $Links_Query_Coverage){
			next;
		}
		
		
		#find the linking relation of current hit with other hits 		
		for(my $j=0; $j<$i; $j++){
			
			#apply the query coverage filter
			if($hitsRef->[$j]->[-1] < $Links_Query_Coverage){
				next;
			}
			
			#establish links which are inter chromosomes only
			if($hitsRef->[$i][1] ne $hitsRef->[$j][1]){
				
				#check if the query sequences of the two hits has any overlap between them
				my $lowerMax = ([$hitsRef->[$i]->[4], $i], [$hitsRef->[$j]->[4], $j])[$hitsRef->[$i]->[4] < $hitsRef->[$j]->[4]];				
				my $upperMin = ([$hitsRef->[$i]->[5], $i], [$hitsRef->[$j]->[5], $j])[$hitsRef->[$i]->[5] > $hitsRef->[$j]->[5]];
				
				
				if($lowerMax->[0] < $upperMin->[0]){
					#overlap is present between two query sequences so the link can be formed between subjects
					#print "$i : @{$hitsRef->[$i]}","\t\t","$j : @{$hitsRef->[$j]}","\n";
					
					#apply the filter: minimum sequence to be matched between two chromosomes 
					if(abs($upperMin->[0] - $lowerMax->[0]) <= $Link_Threshold){
						next;
					}
					
					#print $linkFh "|",abs($upperMin->[0] - $lowerMax->[0]),"|  ";
					
					
					my @link1;
					my @link2;
					#find out the region to match
					
					if($lowerMax->[1] == $upperMin->[1]){
						#when one interval is subinterval of other: eg. (5, 30), (10, 25)
						my $secondIntercal = $lowerMax->[1];
						my $firstInterval = $i;
						if($secondIntercal == $firstInterval){
							$secondIntercal = $j;
						}
						
						#print "@{$hitsRef->[$secondIntercal]}\t\t@{$hitsRef->[$firstInterval]}\n";
						
						@link1 = ($hitsRef->[$secondIntercal]->[1], $hitsRef->[$secondIntercal]->[6], $hitsRef->[$secondIntercal]->[7]);
						@link2 = ($hitsRef->[$firstInterval]->[1], ($hitsRef->[$firstInterval]->[6] + $hitsRef->[$secondIntercal]->[4] - $hitsRef->[$firstInterval]->[4]), ($hitsRef->[$firstInterval]->[7] - ($hitsRef->[$firstInterval]->[5] - $hitsRef->[$secondIntercal]->[5])));
						
					
					}
					else{
						#when intervals have some common region as well as distinct region
						#print "$i : @{$hitsRef->[$i]}","\t\t","$j : @{$hitsRef->[$j]}","\n";
						
						my $firstInterval = $i;
						my $secondIntercal = $j;
						
						#possible conditions which are checked in the if-elsif-else conditions below
						#(5, 20), (5, 30)
						#(5, 30), (15, 30)
						#(5, 20), (5, 20)
						#(5, 20), (15, 30)
						
						if($hitsRef->[$i]->[4] == $hitsRef->[$j]->[4]){
							#check for (*5, 20), (*5, 30) i.e. same start
							
							if($hitsRef->[$i]->[5] == $hitsRef->[$j]->[5]){
								#check for (5, *30), (15, *30) i.e. same end
								#since same start is checked already, this will confirm (5, 20), (5, 20) i.e. same start and same end
								$firstInterval = $i;
								$secondIntercal = $j;
							}
							elsif($hitsRef->[$i]->[5] > $hitsRef->[$j]->[5]){
								#this is for same start but different ends i.e. (*5, 30), (*5, 20)
								$firstInterval = $j;
								$secondIntercal = $i;								
							}
							#no need to handle ($hitsRef->[$i]->[5] < $hitsRef->[$i]->[5]) i.e. (5, 20), (15, 30)
							#because this will endup setting $firstInterval = $i and $secondIntercal = $j which is default above
							
						}
						elsif($hitsRef->[$i]->[4] > $hitsRef->[$j]->[4]){
							#check for different start i.e. (15, 30), (5, 20)
							$firstInterval = $j;
							$secondIntercal = $i;
						}
						
						#no need to handle ($hitsRef->[$i]->[4] < $hitsRef->[$j]->[4]) i.e. (5, 20), (15, 30) and also 
						#no need to handle ($hitsRef->[$i]->[4] < $hitsRef->[$j]->[4]) AND ($hitsRef->[$i]->[5] == $hitsRef->[$j]->[5]) i.e. (5, 30), (15, 30)
						#because this will endup setting $firstInterval = $i and $secondIntercal = $j which is default above
					
						
						@link1 = ($hitsRef->[$firstInterval]->[1], ($hitsRef->[$firstInterval]->[6] + $hitsRef->[$secondIntercal]->[4] - $hitsRef->[$firstInterval]->[4]), $hitsRef->[$firstInterval]->[7]);
						@link2 = ($hitsRef->[$secondIntercal]->[1], $hitsRef->[$secondIntercal]->[6] , ($hitsRef->[$secondIntercal]->[7] - ($hitsRef->[$secondIntercal]->[5] - $hitsRef->[$firstInterval]->[5])));
						
						#print "@link1 @link2\n";
					}
					
					print $linkFh "@link1 @link2 color=$idsRef->{$link2[0]}->[0]\n";
					
				}
			}
			
		}
		
	}
	
}





#create the histogram files and write the data for histogram for each of the file
sub writeHistogramData{
	my $path = shift;
	my $histDataRef = shift;
	
	#print "***$path\n$histDataRef\n";
	#print "### ",join("\t", @{$Hist_CutOffs})," ###\n";
	#craete OUTPUT file handles for the data
	my @outFh;	
	foreach(@{$Hist_CutOffs}){
		push(@{$circosFiles{"HISTOGRAM"}}, "$path/histogram_".$_.".txt");
		
		print "CREATING HISTOGRAM file ",$circosFiles{"HISTOGRAM"}->[-1],"\n";
		open(my $out, ">", $circosFiles{"HISTOGRAM"}->[-1]) or die "Can not create file ",$circosFiles{"HISTOGRAM"}->[-1],": $!";

		push(@outFh, $out);
	}
	
	
	#print the histogram data
	foreach(keys %{$histDataRef}){
		
		foreach my $bin(@{$histDataRef->{$_}}){
			#for each of the %Identity Cutoff, write its frequency to the histogram
			for(my $i = 0; $i < @{$Hist_CutOffs}; $i++){
				print { $outFh[$i] } "@{$bin}[0..2] ",$bin->[3 + $i],"\n";
			}					
		}
	}
	
	
	
	#close the FH after printing
	foreach(@outFh){
		close($_);
	}
}





#read the query file and get the length for each query sequence
sub getQueryLength{
	
	my $file = shift;
	my $storeLocation = shift;
	
	my $outFile = "$storeLocation/queryLen.txt";
	
	#create the circos/data directory to store the karyotye file
	if(!-d $storeLocation){
		mkpath($storeLocation, 1, 0777);
	}
	
	
	#if the query length is already calculated then just read the data from that file directly instead of reading whole query sequences
	if(-e $outFile){
		open(my $fh, $outFile) or die "Cannot open file $outFile: $!";

		my $pattern = qr/(\S+)\s+(\d+)/;
		
		while(<$fh>){
			
			if($_ =~ m/$pattern/){
				$Query_Len{$1} = $2;
			}
		}
		
		close($fh);
		
		#########################		THIS RETURN IS VERY VERY IMPORTANT		############################
		return;
	}
	
	
	#the code below is for reading complete query file and calculate the query sequence length for each query
	#this code is only run if there is no query length file present from previous calculation
	
	#all patterns compiled so code is efficient
	my $fasta = qr/^>(\S+)/;
	my $emptyLine = qr/(^\s+?$)/;
	my $seq = qr/^\w+$/;
		
	open(my $fh, $file) or die "Can not open file $file: $!";
	open(my $out, ">", $outFile) or die "Cannot create file $outFile: $!";
	
	my $seqStart = 0;
	my $id = undef;
	my $len = 0;
	
	while(<$fh>){
		chomp;
		
		if($_ =~ m/$fasta/){
			my $tmpId = $1;
			if(!$id){
				$id = $tmpId;
				next;
			}
			
			
			$Query_Len{$id} = $len;
			print $out $id,"\t",$len,"\n";
			
			
			$len = 0;
			$id = $tmpId;
		}
		elsif($_ =~ m/$seq/){
			$len += length $_;
		}
		elsif($_ =~ m/$emptyLine/){
			next;
		}
	}
	
	#last fasta sequence
	$Query_Len{$id} = $len;
	print $out $id,"\t",$len,"\n";
	
	
	close($fh);
	close($out);
	

}






















1;
