use strict;
use warnings;

my $usage = '

';


my %barcodes = ();

foreach(<DATA>){
	chomp $_;
	my ($id, $seq) = split(/\t/, $_);
	$barcodes{$id} = $seq;
}

# # sample code for sequence mismatch count
# my $s1 = 'ATGCCATGC';
# my $s2 = 'ATGCAAT';

# my $x = $s1 ^ $s2;
# print "|",$x,"|\n";
# print $x | ("\x30" x length($s2)),"::\n";
# $x =~ tr/\0//d;
# print "|",$x,"|\n";
# print length($x),"**\n";
# my $y = $s1 =~ tr/AT//rd;
# print "||",$y,"||\n";

# print length($x),"\n";

# print length($s1) - ( ( $s1 ^ $s2 ) =~ tr/\0// ),"\n";
# my $z = ( ( $s1 ^ $s2 ) =~ tr/\0//dr );
# print "z = |$z|\n";


# print "####\n";
# my $n1 = '010001';
# my $n2 = '10101';
# my $len = 6;
# $n1 = substr((" " x $len).$n1, -$len);
# $n2 = substr((" " x $len).$n2, -$len);
# my $m = ($n1 ^ $n2);
# my $n = ($n1 ^ $n2) | ("\x30" x $len);
# print "\x30" x $len, ":\n";
# print "n1 = $n1\nn2 = $n2\n";
# print $m,"\n";
# print $n,"\n";
# print oct("\0"),"\n";

open(my $fh, $_[0]) or die "File not found: $_[0] ", $!;

while(my $line = <$fh>){
	chomp $line;
	my ($sp, $bc) = split(/\t/, $line);
	my $bcLen = length($sp =~s/^.*_($bc[ATGC]*)_.*/$1/r) + 1;

	# print "$sp\t$bc\t$bcLen\n";
	my $i=3;
	my %freq=();
	
	## frequency of first n base pairs
	open(my $fastq, "zcat $sp | ") or die;	
	while(<$fastq>){
		
		# print $_;
		# last;

		if($i % 4 == 0){
			$freq{substr($_, 0, $bcLen)}++;
		}
		
		$i++;		
	}
	
	my $readCount = ($i - 3) / 4;
	my $count = 0;
	foreach my $ad(sort{$freq{$b} <=> $freq{$a}}keys %freq){
		$count++;

		my $adMis = &index_match($ad, $bc);
		print "frag$i\t$bc\t$ad\t$freq{$ad}\t$readCount\t$adMis\t$ad\n";
		print join("\t", "frag$i", $bc, $ad, $freq{$ad}, $readCount, $adMis, $bc),"\n";
		
		foreach my $bar(sort keys %barcodes){
			my $barMis = &index_match($ad, $barcodes{$bar});
			if($barMis <= $adMis && $bar ne $bc){
				# my $matchPattern = ($ad ^ $barcodes{$bar}) | ("\x30" x length($bc));
				print join("\t", "frag$i", $bc,$ad, $freq{$ad}, $readCount, $barMis, $barcodes{$bar}),"\n";
			}
		}
		
		
		if($count >= 10){
			last;
		}
	}

	close($fastq);
	
	# print "\n";
}

close($fh)



sub index_match{
	my ($s1, $s2) = ('', '');
	
	## ensure that shorter sequence is $s1 and longer is $s2
	if(length($_[0]) < length($_[1])){
		$s1 = $_[0];
		$s2 = $_[1];
	}
	elsif(length($_[0]) > length($_[1])){
		$s1 = $_[1];
		$s2 = $_[0];
	}
	else{
		$s1 = $_[0];
		$s2 = $_[1];
	}
	
	
	my $mismatches = length($s1) - ( ( $s1 ^ $s2 ) =~ tr/\0// );
	
	return $mismatches;
}



## inhouse barcodes
__DATA__
bar01	ATCACGT
bar02	CGATGTT
bar03	TTAGGCT
bar04	TGACCAT
bar05	ACAGTGT
bar06	GCCAATT
bar07	CAGATCT
bar08	ACTTGAT
bar09	GATCAGT
bar10	TAGCTTT
bar11	GGCTACT
bar12	CTTGTAT
bar13	ATATAGGAT
bar14	AACCGTGTT
bar15	AGGTCAGTT
bar16	CTCTGTCTT
bar17	CCATACACT
bar18	CGCATTAAT
bar19	GTCTACATT
bar20	GAGTTAACT
bar21	GCAGCCTCT
bar22	TCGCGTACT
bar23	TATACCGTT
bar24	TGCGGTTAT
bar25	AACACCTACT
bar26	CCTTTACAGT
bar27	GGTCCTTGAT
bar28	TTGAGTGTT
bar29	ACTAACTGCT
bar30	CAGGAGGCGT
bar31	GTTGTCCCAT
bar32	TGACGCATT
bar33	ATCGCCAGCT
bar34	CATTCCAAGT
bar35	GCAAGTAGAT
bar36	TGATCCGAT
bar37	ACGTAGCTCT
bar38	CGAACTGTGT
bar39	TAGCTAGTAT
bar40	GTGGGATAT
bar41	ATCCTATTCT
bar42	CGGACGTGGT
bar43	GCGTTTCGAT
bar44	TATCTCCGT
bar45	ACAGTGCACT
bar46	CACAGTTGGT
bar47	GTGACTACAT
bar48	TGAGAGTGT
bar49	AATGCTGACT
bar50	CCGTCTGAGT
bar51	GGCAGACGAT
bar52	TTCTGATGT
bar53	AGTAGTGGCT
bar54	CTAGTCATGT
bar55	GACACTCTAT
bar56	TCATTAGGT
bar57	TCCAGCCTCT
bar58	CTAGATTCGT
bar59	GAACGCTGAT
bar60	AGAACACCT