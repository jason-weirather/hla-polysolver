#!/usr/bin/perl

use strict;



use Bio::DB::Sam;
use Dumpvalue;
use Data::Dumper;
use List::MoreUtils 'true';
use List::MoreUtils 'indexes';
use List::MoreUtils qw/ uniq /;
use List::Util qw/max/;

########

sub split_lik { # start sub split_lik

	my ($dir,$file,$ref1,$race,$scale) = @_;

	my ($line,@f,$lik,$likLog,%freqHash,%allelePairs,%allelePairsLog,@alleles,%alleleLik,
		%alleleLikLog,@alleles,$allele,$pair,%hashLik,%hashLikLog,$i,$freq,$freqLog);

	%freqHash = %{$ref1};

	open FILE, $dir."/".$file || die "split_lik: Cannot open $file\n";

	while($line = <FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		$allele = $f[0];	
		
		$pair = $f[1];

		$lik = $f[2];

		$likLog = $f[3];

		$allelePairs{$allele}->{$pair} = $lik;

		$allelePairsLog{$allele}->{$pair} = $likLog;

		$alleleLik{$allele} = $alleleLik{$allele} + $lik;

		$alleleLikLog{$allele} = $alleleLikLog{$allele} + $likLog;

	}
	
	close FILE;

	@alleles = keys %allelePairs;

	for($i = 0; $i <= $#alleles; $i++){

		$allele = $alleles[$i];

		#print "$allele\n";

		open FILE2, ">$dir/$allele.all.lik1";

		%hashLik = %{$allelePairs{$allele}};

		%hashLikLog = %{$allelePairsLog{$allele}};
					
		foreach $pair (keys %hashLik){

			print FILE2 $pair,"\t",$hashLik{$pair},"\t",$hashLikLog{$pair},"\n";

		}

		$freq = population_prior_likelihood_component(\%freqHash,$allele,$race,$scale);

		$freqLog = 0;

		if($freq > 0){

			$freqLog = log($freq);

		}
	
		print FILE2 "lik1\t",($alleleLik{$allele}+$freq),"\t",($alleleLikLog{$allele}+$freqLog),"\n";
		
		close FILE2;

		#die;
	}



} # end sub split_lik

########

sub get_likelihood_msf { # start sub get_likelihood_msf

	my ($dir,$pair1,$ref1,$ref2,$ref3,$ref4,$ref5,$s1,$e1,$s2,$e2,$ref6,$ref7,$insert,$ref8,$ref9,$ref10,$ref11) = @_;

	my (@seq1Msf,@qual1Msf,@seq2Msf,@qual2Msf,%msfArr,@allele1Msf,@allele2Msf,$i,%errHash,
		$indelVal,$lik,$likLog,$cSeq,$qSeq,$cRef,$iScore,$iScoreLog,%isizePval,$allele,
		%scoreHash,%scoreHashLog,$matchString,@matchArr,%msf,$len1,$len2,@seg1,@seg2,$seq,
		%scoreSeg1,%scoreSeg2,%scoreSeg1Log,%scoreSeg2Log,$j,$segment1,$segment2,@seg1Uniq,
		@seg2Uniq,@arr,@allelePair,@allelePairLog);

	@seq1Msf = @{$ref1};

	@qual1Msf = @{$ref2};

	@seq2Msf = @{$ref3};

	@qual2Msf = @{$ref4};

	%msfArr = %{$ref5};

	%errHash = %{$ref6};

	%isizePval = %{$ref7};

	%scoreHash = %{$ref8};

	%scoreHashLog = %{$ref9};
	
	%msf = %{$ref10};
	
	@arr = @{$ref11};

	# get the list of unique msf segments
	
	@seg1 = ();

	@seg2 = ();

	$len1 = $e1 - $s1 + 1;

	$len2 = $e2 - $s2 + 1;

	for($i = 0; $i <= $#arr; $i++){

		$seq = $msf{$arr[$i]};
		
		push @seg1, substr($seq,$s1,$len1);

		push @seg2, substr($seq,$s2,$len2);

	}

	@seg1Uniq = uniq(@seg1);

	@seg2Uniq = uniq(@seg2);


	print "unique segments:\t",($#seg1Uniq+1),"\t",($#seg2Uniq+1),"\n";

	#print "segs:\t@seg1\n@seg2\n";

	# score the unique segments

	for($j = 0; $j <= $#seg1Uniq; $j++){

		$segment1 = $seg1Uniq[$j];
	
		@allele1Msf = split(//,$segment1);

		$lik = 0;

		$likLog = 0;

		for($i = 0; $i <= $#qual1Msf; $i++){

			$lik = $lik + $scoreHash{$qual1Msf[$i]}->{$seq1Msf[$i].$allele1Msf[$i]};

			$likLog = $likLog + $scoreHashLog{$qual1Msf[$i]}->{$seq1Msf[$i].$allele1Msf[$i]};

		}

		$scoreSeg1{$segment1} = $lik;

		$scoreSeg1Log{$segment1} = $likLog;

		#print "$segment1\t$lik\t$likLog\n";

	}


	for($j = 0; $j <= $#seg2Uniq; $j++){

		$segment2 = $seg2Uniq[$j];

		@allele2Msf = split(//,$segment2);

		$lik = 0;

		$likLog = 0;

		for($i = 0; $i <= $#qual2Msf; $i++){

			$lik = $lik + $scoreHash{$qual2Msf[$i]}->{$seq2Msf[$i].$allele2Msf[$i]};

			$likLog = $likLog + $scoreHashLog{$qual2Msf[$i]}->{$seq2Msf[$i].$allele2Msf[$i]};

		}

		$scoreSeg2{$segment2} = $lik;

		$scoreSeg2Log{$segment2} = $likLog;

		#print "$segment2\t$lik\t$likLog\n";

	}

	#system("date");

	# get the score for each allele for the given read pair
	
	@allelePair = ();

	@allelePairLog = ();


	for($i = 0; $i <= $#arr; $i++) {

		$allele = $arr[$i];

		$segment1 = $seg1[$i];
		
		$segment2 = $seg2[$i];

		$iScore = $isizePval{abs($insert)};
	
		if($iScore > 0){
	
			$iScoreLog = log($iScore);
	
		}	
		
		$lik = $scoreSeg1{$segment1} + $scoreSeg2{$segment2};
	
		$likLog = $scoreSeg1Log{$segment1} + $scoreSeg2Log{$segment2};
		
		push @allelePair, $lik;

		push @allelePairLog, $likLog;
		
		#open FILE,">>$dir/$allele.all.lik1";

		#print FILE "$pair1\t$lik\t$likLog\n";

		#close FILE;

		#`echo "$pair1\t$lik\t$likLog" >> $dir/$allele.all.lik1`;

		#@allele1Msf = @{$msfArr{$allele}}[$s1..$e1];
	
		#@allele2Msf = @{$msfArr{$allele}}[$s2..$e2];
		
		
		#if($#allele1Msf != $#seq1Msf || $#allele2Msf != $#seq2Msf){
	
			#die "get_likelihood_msf: lengths are not equal $#allele1Msf != $#seq1Msf || $#allele2Msf != $#seq2Msf\n";
		
		#}
	
		#$iScore = $isizePval{abs($insert)};
	
		#if($iScore > 0){
	
		#	$iScoreLog = log($iScore);
	
		#}
	
		#$lik = 0;
	
		#$likLog = 0;
	
		#@matchArr = @{get_match_string(\@seq1Msf,\@allele1Msf)};

		#for($i = 0; $i <= $#qual1Msf; $i++){

		#	$lik = $lik + $scoreHash{$qual1Msf[$i]}->{$seq1Msf[$i].$allele1Msf[$i]};

		#	$likLog = $likLog + $scoreHashLog{$qual1Msf[$i]}->{$seq1Msf[$i].$allele1Msf[$i]};

		#}
	
		#@matchArr = @{get_match_string(\@seq2Msf,\@allele2Msf)};

		#for($i = 0; $i <= $#qual2Msf; $i++){

		#	$lik = $lik + $scoreHash{$qual2Msf[$i]}->{$seq2Msf[$i].$allele2Msf[$i]};

		#	$likLog = $likLog + $scoreHashLog{$qual2Msf[$i]}->{$seq2Msf[$i].$allele2Msf[$i]};

		#}

		#$lik = $lik + $iScore;
	
		#$likLog = $likLog + $iScoreLog;
		
		#`echo "$pair1\t$lik\t$likLog" >> $dir/$allele.all.lik1`;

	}

	return (\@allelePair,\@allelePairLog);

} # end sub get_likelihood_msf


########

sub plot_reads { # start sub plot_reads

	my ($file,$ref1,$ref2) = @_; # ex. file = NA18504/hla_b_reads.txt

	my (@lines,$i,$line1,$line2,$pair1,$allele1,$start1,$insert1,$seq1,$qual1,$len1,$pair2,$allele2,$start2,$insert2,$seq2,$qual2,$len2,$ref3,$ref4,@seq1Msf,@qual1Msf,@seq2Msf,@qual2Msf,@f1,@f2,%seq2MsfPos,%msf2SeqPos,$segment1,$segment2,$alleleStart1,$alleleEnd1,$alleleStartMsf1,$alleleEndMsf1,$alleleStart2,$alleleEnd2,$alleleStartMsf2,$alleleEndMsf2,$pad1,$pad2);

	%seq2MsfPos = %{$ref1};
	
	%msf2SeqPos = %{$ref2};	
	
	@lines =  `cat $file`;

	print "plot_reads: # lines = ",($#lines + 1),"\n";

	open FILE, $file;

	for($i=0;$i<=$#lines;$i++){
		#print "i=$i\n";
		#system("date");
		$line1 = $lines[$i];
		$line2 = $lines[$i+1];
		chomp($line1);
		chomp($line2);
		$i++;
		@f1 = split(/\t/,$line1);
		@f2 = split(/\t/,$line2);
		$pair1 = $f1[0];
		$allele1 = $f1[2];
		$start1 = $f1[3];
		$insert1 = $f1[8];
		$seq1 = $f1[9];
		$qual1 = $f1[10];
		$len1 = length($seq1);
		$pair2 = $f2[0];
		$allele2 = $f2[2];
		$start2 = $f2[3];
		$insert2 = $f2[8];
		$seq2 = $f2[9];
		$qual2 = $f2[10];
		$len2 = length($seq2);
		
		#print "allele1\t$allele1\t$start1\t$len1\t$start2\t$len2\n";

		# insert1,insert2 allele1,allele2 and pair1,pair2 should be the same; previously checked
	
		# get seq1Msf, qual1Msf, seq2Msf, qual2Msf
	
		($ref1,$ref2,$ref3,$ref4) = get_seq_msfs($allele1,$start1,$len1,$seq1,$qual1,
							$start2,$len2,$seq2,$qual2,\%seq2MsfPos,\%msf2SeqPos);
		@seq1Msf = @{$ref1};
		@qual1Msf = @{$ref2};
		@seq2Msf = @{$ref3};
		@qual2Msf = @{$ref4};
	
		#print "seq1Msf = ",@seq1Msf,"\n";
		#print "count = ", ($#seq1Msf+1),"\t",($#seq2Msf+1),"\n";
	
		$segment1 = join('',@seq1Msf);

		$segment2 = join('',@seq2Msf);	
	
		$segment1 =~ s/-1/\./g;

		$segment2 =~ s/-1/\./g;

		#print "$segment1\n$segment2\n";

		$alleleStart1 = $start1 - 1; 					# 0-base
		$alleleEnd1 = $start1 + $len1 - 2; 				# 0-base
		$alleleStartMsf1 = $seq2MsfPos{$allele1}->{$alleleStart1};	# 0-base 
		$alleleEndMsf1 = $seq2MsfPos{$allele1}->{$alleleEnd1};		# 0-base

		$alleleStart2 = $start2 - 1; 					# 0-base
		$alleleEnd2 = $start2 + $len2 - 2; 				# 0-base
		$alleleStartMsf2 = $seq2MsfPos{$allele2}->{$alleleStart2};	# 0-base 
		$alleleEndMsf2 = $seq2MsfPos{$allele2}->{$alleleEnd2};		# 0-base

		$pad1 = " " x $alleleStartMsf1;

		$pad2 = " " x $alleleStartMsf2;

		print "$pad1$segment1\n$pad2$segment2\n";

	}

	close FILE;



} # end sub plot_reads

########

sub get_match_string { # start sub get_match_string

	my ($ref1,$ref2) = @_;

	my (@seqMsf,@alleleMsf,$i,$matchString,@matchArr,$cSeq,$cRef);

	@seqMsf = @{$ref1};

	@alleleMsf = @{$ref2};

	@matchArr = ();

	for($i = 0; $i <= $#seqMsf; $i++){
		
		$cSeq = $seqMsf[$i];

		$cRef = $alleleMsf[$i];		
		
		if($cSeq ne "-1" & $cRef ne "."){

			if($cSeq eq $cRef){

				push @matchArr,"M";
	
			} else {

				push @matchArr,"m";

			}

		} elsif($cSeq eq "-1" & $cRef eq "."){

			push @matchArr,"b";	

		} else {

			push @matchArr,"b";

		}

	}

	return \@matchArr;

} # end sub get_match_string

########

sub get_seq_msfs { # start sub get_seq_msfs

	my ($allele1,$start1,$len1,$seq1,$qual1,$start2,$len2,$seq2,$qual2,$ref1,$ref2) = @_;

	my (%seq2MsfPos,%msf2SeqPos,$alleleStart,$alleleEnd,$alleleStartMsf,
			$alleleEndMsf,@seq1Arr,@qual1Arr,@seq2Arr,@qual2Arr,@seq1Msf,@qual1Msf,
			@seq2Msf,@qual2Msf,$k,$pos,$posMsf,@arr);

	%seq2MsfPos = %{$ref1};						# 0-base
	%msf2SeqPos = %{$ref2};						# 0-base

	$alleleStart = $start1 - 1; 					# 0-base
	$alleleEnd = $start1 + $len1 - 2; 				# 0-base
	$alleleStartMsf = $seq2MsfPos{$allele1}->{$alleleStart};	# 0-base 
	$alleleEndMsf = $seq2MsfPos{$allele1}->{$alleleEnd};		# 0-base

	@seq1Arr = split(//,$seq1);
	@qual1Arr = split(//,$qual1);
	@seq1Msf = ();
	@qual1Msf = ();
	$k = 0;

	#print "alleleStart = $alleleStart\talleleEnd = $alleleEnd\talleleStartMsf = $alleleStartMsf\talleleEndMsf = $alleleEndMsf\tposMsf = $posMsf\n";

	for($posMsf = $alleleStartMsf; $posMsf <= $alleleEndMsf; $posMsf++){
		$pos = $msf2SeqPos{$allele1}->{$posMsf};		# posMsf and pos are 0-base
		#print "$posMsf\t$pos\t$k\t$seq1Arr[$k]\n";
		if($pos != -1){
			push @seq1Msf, $seq1Arr[$k];
			push @qual1Msf, $qual1Arr[$k];
			$k++;
		} else {
			push @seq1Msf, -1;
			push @qual1Msf, -1;
			#print "$allele1\t$pair1\t$alleleStartMsf\t$alleleEndMsf\n";
		}
		#print "$seq1Msf[$#seq1Msf]\t$qual1Msf[$#qual1Msf]\n";

	}

	# get seq2Msf, qual2Msf

	$alleleStart = $start2 - 1; 
	$alleleEnd = $start2 + $len2 - 2;
	$alleleStartMsf = $seq2MsfPos{$allele1}->{$alleleStart};
	$alleleEndMsf = $seq2MsfPos{$allele1}->{$alleleEnd};

	@seq2Arr = split(//,$seq2);
	@qual2Arr = split(//,$qual2);
	@seq2Msf = ();
	@qual2Msf = ();
	$k = 0;

	for($posMsf = $alleleStartMsf; $posMsf <= $alleleEndMsf; $posMsf++){
		$pos = $msf2SeqPos{$allele1}->{$posMsf};
		if($pos != -1){
			push @seq2Msf, $seq2Arr[$k];
			push @qual2Msf, $qual2Arr[$k];
			$k++;
		} else {
			push @seq2Msf, -1;
			push @qual2Msf, -1;
			#print "$allele1\t$pair1\t$alleleStartMsf\t$alleleEndMsf\n";
		}

	}



	return (\@seq1Msf,\@qual1Msf,\@seq2Msf,\@qual2Msf);;

} # end sub get_seq_msfs


########

sub get_msf2Seq_pos{ # start sub get_msf2Seq_pos

	# function gets the mapping from msf to seq coordinates
		# msf2SeqPos{name}->{msfPos} = seqPos
	# also gets the mapping from seq to msf coordinate
		# seq2MsfPos{name}->{seqPos} = msfPos
	# function also returns the a hash of lengths of the different alleles
	# assumes that there is only sequence in each fasta file

	my ($ref1,$ref2,$ref3) = @_;

	my (@names,%msfHashArr,%msf2SeqPos,$i,$fastaFile,%hash,$seq,@seqArr,$key,@msfArr);

	my ($msfBase,$j,$k,$c,%seqLenHash,$length,%seq2MsfPos,%seqHash,$name);

	@names = @$ref1;

	%seqHash = %$ref2;

	%msfHashArr = %$ref3;

	%msf2SeqPos = ();

	for($i=0;$i<=$#names;$i++){

		$name = $names[$i];

		#print "name=$name\n";

		# get sequence

		$seq = $seqHash{$name};

		# get mapping from msf to seq coordinates

		@msfArr = @{$msfHashArr{$name}};
	
		@seqArr = split(//,$seq);

		$length = $#seqArr + 1;

		$seqLenHash{$name} = $length;

		#print "@msfArr\n@seqArr\n";
		#die;

		$k = 0;

		for($j=0;$j<=$#msfArr;$j++){

			$msfBase = $msfArr[$j];
			
			if($msfBase eq "."){

				$msf2SeqPos{$name}->{$j} = -1;

				next;

			}

			for(; $k<=$#seqArr;){

				if($msfBase eq $seqArr[$k]){
			
					$msf2SeqPos{$name}->{$j} = $k;

					$seq2MsfPos{$name}->{$k} = $j;

					$k++;

					last;

				}

			}

		}


	}


	return ( \%msf2SeqPos, \%seq2MsfPos, \%seqLenHash );

} # end sub get_msf2Seq_pos


########

sub get_read_alleles { # start sub get_read_alleles

	my ($dir,$bam) = @_;

	my (@lines,$line1,$line2,$i,@f1,@f2,%aReads,%bReads,%cReads,$read,$allele,$start1,$start2,
		$outFileA,$outFileB,$outFileC,$insert1,$insert2,$insertStatus,%alleleCounts,
		@aAlleles,@bAlleles,@cAlleles);

	$outFileA = $dir."/hla_a_reads.txt";
	$outFileB = $dir."/hla_b_reads.txt";
	$outFileC = $dir."/hla_c_reads.txt";

	open OUTFILEA, ">$outFileA" || die "Cannot open $outFileA\n";
	open OUTFILEB, ">$outFileB" || die "Cannot open $outFileB\n";
	open OUTFILEC, ">$outFileC" || die "Cannot open $outFileC\n";

	@lines = `samtools view $bam`;

	for($i=0;$i<=$#lines;$i++){
		$line1 = $lines[$i];
		$line2 = $lines[$i+1];
		chomp($line1);
		chomp($line2);
		@f1 = split(/\t/,$line1);
		@f2 = split(/\t/,$line2);
		if($f1[0] ne $f2[0]){
			next;
		}
		$i++;
		if($f1[2] ne $f2[2] || $f1[5]=~/I/ || $f1[5]=~/D/ || $f1[5]=~/\*/ || $f2[5]=~/I/ || $f2[5]=~/D/ || $f2[5]=~/\*/ || $f1[6] ne "=" || $f2[6] ne "="){
			next;
		}	
		$read = $f1[0];
		$allele = $f1[2];
		$start1 = $f1[3];
		$start2 = $f2[3];
		$insert1 = abs($f1[8]);
		$insert2 = abs($f2[8]);
		$insertStatus = 1;
		if($insert1 != $insert2){
			#print "get_read_alleles: ERROR: insertsizes unequal:\n$line1\n$line2\n";
			$insertStatus = 0;
		}

		$alleleCounts{$allele} = $alleleCounts{$allele} + 1;

		#if($allele eq "hla_a_01_07"){

			#print "$line1\t$line2\n";

		#}

		if($allele =~ /hla_a/ & $aReads{$read} != 1 & $insertStatus == 1){
			$aReads{$read} = 1;
			print OUTFILEA "$line1\n$line2\n";
		}
		if($allele =~ /hla_b/ & $bReads{$read} != 1 & $insertStatus == 1){
			$bReads{$read} = 1;
			print OUTFILEB "$line1\n$line2\n";
		}
		if($allele =~ /hla_c/ & $cReads{$read} != 1 & $insertStatus == 1){
			$cReads{$read} = 1;
			print OUTFILEC "$line1\n$line2\n";
		}

	}

	close OUTFILEA;
	close OUTFILEB;
	close OUTFILEC;

	$outFileA = $dir."/hla_a_alleles.txt";
	$outFileB = $dir."/hla_b_alleles.txt";
	$outFileC = $dir."/hla_c_alleles.txt";

	open OUTFILEA, ">$outFileA" || die "Cannot open $outFileA\n";
	open OUTFILEB, ">$outFileB" || die "Cannot open $outFileB\n";
	open OUTFILEC, ">$outFileC" || die "Cannot open $outFileC\n";


	foreach $allele (keys %alleleCounts){

		if($alleleCounts{$allele} > 0){

			if($allele =~ /hla_a/){

				push @aAlleles, $allele;
	
				print OUTFILEA $allele,"\t",$alleleCounts{$allele},"\n";

			}elsif($allele =~ /hla_b/){

				push @bAlleles, $allele;

				print OUTFILEB $allele,"\t",$alleleCounts{$allele},"\n";
	
			}elsif($allele =~ /hla_c/){

				push @cAlleles, $allele;

				print OUTFILEC $allele,"\t",$alleleCounts{$allele},"\n";
	
			}
		}		
		
	}

	close OUTFILEA;
	close OUTFILEB;
	close OUTFILEC;

	return 1;

} # end sub get_read_alleles

########

sub get_alleles { # start sub get_alleles

	my ($dir) = @_;

	my (@aAlleles,@bAlleles,@cAlleles,$i,$line,@f);

	open FILE, "$dir/hla_a_alleles.txt";

	while($line = <FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		push @aAlleles, $f[0];

	}

	close FILE;	

	open FILE, "$dir/hla_b_alleles.txt";

	while($line = <FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		push @bAlleles, $f[0];

	}

	close FILE;	


	open FILE, "$dir/hla_c_alleles.txt";

	while($line = <FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		push @cAlleles, $f[0];

	}

	close FILE;	


	return (\@aAlleles, \@bAlleles, \@cAlleles);

} # end sub get_alleles


########

sub get_msf { # start sub get_msf

	my ($file) = @_;

	my ($name,$i,%seq,%seqArr,%msfSeq,%msfArr,$s1,$s2,@f1,@f2);

	%seq  = ();

	%seqArr = ();

	%msfSeq = ();

	%msfArr = ();

	open FILE, $file || die "get_msf: Cannot open $file\n";

	while($name = <FILE>){

		chomp($name);
	
		$name = substr($name,1);

		$s1 = <FILE>;

		chomp($s1);

		$s2 = $s1;

		$s2 =~ s/\.+//g;	

		@f1 = split(//,$s1);

		@f2 = split(//,$s2);

		$seq{$name} = $s2;

		$msfSeq{$name} = $s1;

		$seqArr{$name} = [ @f2 ];
	
		$msfArr{$name} = [ @f1 ];

	}

	close FILE;

	return (\%seq,\%msfSeq,\%seqArr,\%msfArr);


} # end sub get_msf

########

sub get_seq_and_qual { # start sub get_seq_and_qual

	my ($dir) = @_;

	my ($seq1,$seq2,$qual1,$qual2,%seq1Hash,%seq2Hash,%qual1Hash,%qual2Hash,$i,$file1,$file2,$name1,$name2,@arr);

	$file1 = $dir."/chr6region.1.fastq";
	$file2 = $dir."/chr6region.2.fastq";

	open FILE1, $file1;
	open FILE2, $file2;

	%seq1Hash = ();
	%seq2Hash = ();
	%qual1Hash = ();
	%qual2Hash = ();
	while($name1 = <FILE1>){
		$name2 = <FILE2>;
		chomp($name1);
		chomp($name2);
		$name1 = substr($name1,1,length($name1)-3);
		$name2 = substr($name2,1,length($name2)-3);
		if($name1 ne $name2){
			die "get_seq_and_qual: name1 = $name1 ne name2 = $name2\n";
		}
		$seq1 = <FILE1>; 
		chomp($seq1);
		$seq2 = <FILE2>;
		chomp($seq2);
		<FILE1>; <FILE2>;
		$qual1 = <FILE1>;
		$qual2 = <FILE2>;
		chomp($qual1);
		chomp($qual2);
		$seq1Hash{$name1} = $seq1;
		$seq2Hash{$name1} = $seq2;
		$qual1Hash{$name1} = $qual1;
		$qual2Hash{$name1} = $qual2;
	}


	@arr = (\%seq1Hash,\%seq2Hash,\%qual1Hash,\%qual2Hash);


	return \@arr;

} # end sub get_seq_and_qual


########

sub get_likelihood_log { # start sub get_likelihood_log

	my ($ref1,$ref2,$q,$md) = @_;

	my (%scoreHash,%scoreHashLog,@qArr,@mdArr,$i,$lik,$likLog,@arr,$likMatch,$likMismatch,$likLogMatch,$likLogMismatch);

	%scoreHash = %{$ref1};

	%scoreHashLog = %{$ref2};
	
	@qArr = split(//,$q);
	
	@mdArr = split(//,$md);

	if($#qArr != $#mdArr){

		die "get_likelihood: $q is not the same length as $md\n";
	
	}

	$lik = 0;
	
	$likLog = 0;

        $likMatch = 0;

        $likLogMatch = 0;

        $likMismatch = 0;

        $likLogMismatch = 0;
	
	for($i = 0; $i <= $#qArr; $i++){

		#print $qArr[$i],"\t",$mdArr[$i],"\t",$scoreHash{$qArr[$i]}->{$mdArr[$i]},"\n";

		$lik = $lik + $scoreHash{$qArr[$i]}->{$mdArr[$i]};		
		
		$likLog = $likLog + $scoreHashLog{$qArr[$i]}->{$mdArr[$i]};

		if($mdArr[$i] eq "M"){

			$likMatch = $likMatch + $scoreHash{$qArr[$i]}->{$mdArr[$i]};

                      	$likLogMatch = $likLogMatch + $scoreHashLog{$qArr[$i]}->{$mdArr[$i]};

		} 

               if($mdArr[$i] eq "m"){

                       $likMismatch = $likMismatch + $scoreHash{$qArr[$i]}->{$mdArr[$i]};

                       $likLogMismatch = $likLogMismatch + $scoreHashLog{$qArr[$i]}->{$mdArr[$i]};

                }


		#print "$qArr[$i]\t$mdArr[$i]\t$lik\t$likLog\n";

	}
	
	@arr = ($lik,$likLog,$likMatch,$likLogMatch,$likMismatch,$likLogMismatch);

	return \@arr;

} # end sub get_likelihood_log


########

sub get_likelihood { # start sub get_likelihood

	my ($ref,$q,$md) = @_;

	my (%scoreHash,@qArr,@mdArr,$i,$lik,$likLog,@arr);

	%scoreHash = %{$ref};
	
	@qArr = split(//,$q);
	
	@mdArr = split(//,$md);

	if($#qArr != $#mdArr){

		die "get_likelihood: $q is not the same length as $md\n";
	
	}

	$lik = 0;
	
	$likLog = 0;
	
	for($i = 0; $i <= $#qArr; $i++){

		#print $qArr[$i],"\t",$mdArr[$i],"\t",$scoreHash{$qArr[$i]}->{$mdArr[$i]},"\n";

		$lik = $lik + $scoreHash{$qArr[$i]}->{$mdArr[$i]};		
		
		$likLog = $likLog + log($scoreHash{$qArr[$i]}->{$mdArr[$i]});

		print "$qArr[$i]\t$mdArr[$i]\t$lik\t$likLog\n";

	}
	
	@arr = ($lik,$likLog);

	return \@arr;

} # end sub get_likelihood

########

sub md_to_match_string { # start sub md_to_match_string

	my ($md) = @_;

	my ($i,$string,@f);	

	if($md =~ "MD:Z:0"){

                $md = substr($md,6);

        } else {

                $md = substr($md,5);

        }       

        while($md =~ /[A-Z]0[A-Z]/){

                $md =~ s/[A-Z]0[A-Z]/CC/g;

        }
        
        #print "$md\n";
        
        $md =~ s/[A-Z]/C/g;

        #print "$md\n";

        @f = split(/C/,$md);

        $string="";

        for($i=0;$i<=$#f;$i++){

                #print "$i\t$f[$i]\n";

                if($f[$i] eq ""){

                        $string = $string."m";

                }else{

                        $string = $string.("M" x $f[$i]);

                        $string = $string."m";

                }

                #print "$string\n";

        }       

        chop($string);

	return $string;	

} # end sub md_to_match_string

########

sub get_score_hash_both { # start sub get_score_hash_both

	my ($ref,$indelVal) = @_;

	my (%errHash,$key,$val,%scoreHash,%scoreHashLog);

	%errHash = %{$ref};

	foreach $key (keys %errHash){

		$val = $errHash{$key};

		$scoreHash{$key}->{'AA'} = 1 - $val;
		$scoreHash{$key}->{'CC'} = 1 - $val;
		$scoreHash{$key}->{'GG'} = 1 - $val;
		$scoreHash{$key}->{'TT'} = 1 - $val;

		$scoreHash{$key}->{'AC'} = $val/3;
		$scoreHash{$key}->{'AG'} = $val/3;
		$scoreHash{$key}->{'AT'} = $val/3;
		$scoreHash{$key}->{'A.'} = $indelVal;
		$scoreHash{$key}->{'CA'} = $val/3;
		$scoreHash{$key}->{'CG'} = $val/3;
		$scoreHash{$key}->{'CT'} = $val/3;
		$scoreHash{$key}->{'C.'} = $indelVal;
		$scoreHash{$key}->{'GA'} = $val/3;
		$scoreHash{$key}->{'GC'} = $val/3;
		$scoreHash{$key}->{'GT'} = $val/3;
		$scoreHash{$key}->{'G.'} = $indelVal;
		$scoreHash{$key}->{'TA'} = $val/3;
		$scoreHash{$key}->{'TC'} = $val/3;
		$scoreHash{$key}->{'TG'} = $val/3;
		$scoreHash{$key}->{'T.'} = $indelVal;
		$scoreHash{$key}->{'.A'} = $indelVal;
		$scoreHash{$key}->{'.C'} = $indelVal;
		$scoreHash{$key}->{'.G'} = $indelVal;
		$scoreHash{$key}->{'.T'} = $indelVal;
		$scoreHash{$key}->{'-1A'} =  $indelVal;
		$scoreHash{$key}->{'-1C'} =  $indelVal;
		$scoreHash{$key}->{'-1G'} = $indelVal;
		$scoreHash{$key}->{'-1T'} = $indelVal;
		$scoreHash{$key}->{'A-1'} =  $indelVal;
		$scoreHash{$key}->{'C-1'} =  $indelVal;
		$scoreHash{$key}->{'G-1'} = $indelVal;
		$scoreHash{$key}->{'T-1'} = $indelVal;

		$scoreHash{$key}->{'.-1'} = 0;
		$scoreHash{$key}->{'-1.'} = 0;


		if($val != 1) {

			$scoreHashLog{$key}->{'AA'} = log(1 - $val);
			$scoreHashLog{$key}->{'CC'} = log(1 - $val);
			$scoreHashLog{$key}->{'GG'} = log(1 - $val);
			$scoreHashLog{$key}->{'TT'} = log(1 - $val);

		} else {

			$scoreHashLog{$key}->{'AA'} = log(1-0.99);
			$scoreHashLog{$key}->{'CC'} = log(1-0.99);
			$scoreHashLog{$key}->{'GG'} = log(1-0.99);
			$scoreHashLog{$key}->{'TT'} = log(1-0.99);
				
		}

		$scoreHashLog{$key}->{'AC'} = log($val/3);
		$scoreHashLog{$key}->{'AG'} = log($val/3);
		$scoreHashLog{$key}->{'AT'} = log($val/3);
		$scoreHashLog{$key}->{'A.'} = log($indelVal);
		$scoreHashLog{$key}->{'CA'} = log($val/3);
		$scoreHashLog{$key}->{'CG'} = log($val/3);
		$scoreHashLog{$key}->{'CT'} = log($val/3);
		$scoreHashLog{$key}->{'C.'} = log($indelVal);
		$scoreHashLog{$key}->{'GA'} = log($val/3);
		$scoreHashLog{$key}->{'GC'} = log($val/3);
		$scoreHashLog{$key}->{'GT'} = log($val/3);
		$scoreHashLog{$key}->{'G.'} = log($indelVal);
		$scoreHashLog{$key}->{'TA'} = log($val/3);
		$scoreHashLog{$key}->{'TC'} = log($val/3);
		$scoreHashLog{$key}->{'TG'} = log($val/3);
		$scoreHashLog{$key}->{'T.'} = log($indelVal);
		$scoreHashLog{$key}->{'.A'} = log($indelVal);
		$scoreHashLog{$key}->{'.C'} = log($indelVal);
		$scoreHashLog{$key}->{'.G'} = log($indelVal);
		$scoreHashLog{$key}->{'.T'} = log($indelVal);
		$scoreHashLog{$key}->{'-1A'} =  log($indelVal);
		$scoreHashLog{$key}->{'-1C'} =  log($indelVal);
		$scoreHashLog{$key}->{'-1G'} = log($indelVal);
		$scoreHashLog{$key}->{'-1T'} = log($indelVal);
		$scoreHashLog{$key}->{'A-1'} =  log($indelVal);
		$scoreHashLog{$key}->{'C-1'} =  log($indelVal);
		$scoreHashLog{$key}->{'G-1'} = log($indelVal);
		$scoreHashLog{$key}->{'T-1'} = log($indelVal);

		$scoreHashLog{$key}->{'.-1'} = 0;
		$scoreHashLog{$key}->{'-1.'} = 0;

	}


	return (\%scoreHash, \%scoreHashLog);

} # end sub get_score_hash_both

########

sub get_score_hash_log { # start sub get_score_hash_log

	my ($ref,$scale) = @_;

	my (%errHash,$key,$val,%scoreHash,%scoreHashLog);

	%errHash = %{$ref};

	foreach $key (keys %errHash){

		$val = $errHash{$key};

		$scoreHash{$key}->{'M'} = 1 - $val;

		$scoreHash{$key}->{'m'} = $val/3;

		if($val != 1){

			$scoreHashLog{$key}->{'M'} = log($scale*(1 - $val));

			$scoreHashLog{$key}->{'m'} = log($scale*($val/3));

	                #$scoreHashLog{$key}->{'M'} = log(1 - $val) + log($scale);

                        #$scoreHashLog{$key}->{'m'} = log($val/3) + log($scale);

	
		} else {

			$scoreHashLog{$key}->{'M'} = 0;

			$scoreHashLog{$key}->{'m'} = log($scale*($val/3));

                        #$scoreHashLog{$key}->{'M'} = 0 + log($scale);

                        #$scoreHashLog{$key}->{'m'} = log($val/3) + log($scale);


		}	

	}

	return (\%scoreHash,\%scoreHashLog);

} # end sub get_score_hash_log


########

sub get_score_hash { # start sub get_score_hash

	my ($ref) = @_;

	my (%errHash,$key,$val,%scoreHash);

	%errHash = %{$ref};

	foreach $key (keys %errHash){

		$val = $errHash{$key};

		$scoreHash{$key}->{'M'} = 1 - $val;

		$scoreHash{$key}->{'m'} = $val/3;


	}

	return \%scoreHash;

} # end sub get_score_hash

########

sub get_err_hash_illumina { # start sub get_err_hash_illumina

	my (%errHash);

	$errHash{'!'} = 10**(-0/10);
	$errHash{'"'} = 10**(-1/10);
	$errHash{'#'} = 10**(-2/10);
	$errHash{'$'} = 10**(-3/10);
	$errHash{'%'} = 10**(-4/10);
	$errHash{'&'} = 10**(-5/10);
	$errHash{'\''} = 10**(-6/10);
	$errHash{'('} = 10**(-7/10);
	$errHash{')'} = 10**(-8/10);
	$errHash{'*'} = 10**(-9/10);
	$errHash{'+'} = 10**(-10/10);
	$errHash{','} = 10**(-11/10);
	$errHash{'-'} = 10**(-12/10);
	$errHash{'.'} = 10**(-13/10);
	$errHash{'/'} = 10**(-14/10);
	$errHash{'0'} = 10**(-15/10);
	$errHash{'1'} = 10**(-16/10);
	$errHash{'2'} = 10**(-17/10);
	$errHash{'3'} = 10**(-18/10);
	$errHash{'4'} = 10**(-19/10);
	$errHash{'5'} = 10**(-20/10);
	$errHash{'6'} = 10**(-21/10);
	$errHash{'7'} = 10**(-22/10);
	$errHash{'8'} = 10**(-23/10);
	$errHash{'9'} = 10**(-24/10);
	$errHash{':'} = 10**(-25/10);
	$errHash{';'} = 10**(-26/10);
	$errHash{'<'} = 10**(-27/10);
	$errHash{'='} = 10**(-28/10);
	$errHash{'>'} = 10**(-29/10);
	$errHash{'?'} = 10**(-30/10);
	$errHash{'@'} = 10**(-31/10);
	$errHash{'A'} = 10**(-32/10);
	$errHash{'B'} = 10**(-33/10);
	$errHash{'C'} = 10**(-34/10);
	$errHash{'D'} = 10**(-35/10);
	$errHash{'E'} = 10**(-36/10);
	$errHash{'F'} = 10**(-37/10);
	$errHash{'G'} = 10**(-38/10);
	$errHash{'H'} = 10**(-39/10);
	$errHash{'I'} = 10**(-40/10);
	$errHash{'J'} = 10**(-41/10);
	$errHash{'K'} = 10**(-42/10); # should this really be there?
	$errHash{'L'} = 10**(-43/10); # should this really be there?
	$errHash{'M'} = 10**(-44/10);
        $errHash{'N'} = 10**(-45/10);
	$errHash{'O'} = 10**(-46/10);
	$errHash{'P'} = 10**(-47/10);
	$errHash{'Q'} = 10**(-48/10);
	$errHash{'R'} = 10**(-49/10);
	$errHash{'S'} = 10**(-50/10);
	$errHash{'T'} = 10**(-51/10);
	$errHash{'U'} = 10**(-52/10);
	$errHash{'V'} = 10**(-53/10);
	$errHash{'W'} = 10**(-54/10);
	$errHash{'X'} = 10**(-55/10);
	$errHash{'Y'} = 10**(-56/10);
	$errHash{'Z'} = 10**(-57/10);
	$errHash{'['} = 10**(-58/10);
	$errHash{'\\'} = 10**(-59/10);
	$errHash{']'} = 10**(-60/10);
	$errHash{'^'} = 10**(-61/10);
	$errHash{'_'} = 10**(-62/10);
	$errHash{'`'} = 10**(-63/10);
	$errHash{'a'} = 10**(-64/10);
	$errHash{'b'} = 10**(-65/10);
	$errHash{'c'} = 10**(-66/10);
	$errHash{'d'} = 10**(-67/10);
	$errHash{'e'} = 10**(-68/10);
	$errHash{'f'} = 10**(-69/10);
	$errHash{'g'} = 10**(-70/10);
	$errHash{'h'} = 10**(-71/10);
	$errHash{'i'} = 10**(-72/10);
	$errHash{'j'} = 10**(-73/10);
	$errHash{'k'} = 10**(-74/10);
	$errHash{'l'} = 10**(-75/10);
	$errHash{'m'} = 10**(-76/10);
	$errHash{'n'} = 10**(-77/10);
	$errHash{'o'} = 10**(-78/10);
	$errHash{'p'} = 10**(-79/10);
	$errHash{'q'} = 10**(-80/10);
	$errHash{'r'} = 10**(-81/10);
	$errHash{'s'} = 10**(-82/10);
	$errHash{'t'} = 10**(-83/10);
	$errHash{'u'} = 10**(-84/10);
	$errHash{'v'} = 10**(-85/10);
	$errHash{'w'} = 10**(-86/10);
	$errHash{'x'} = 10**(-87/10);
	$errHash{'y'} = 10**(-88/10);
	$errHash{'z'} = 10**(-89/10);
	$errHash{'{'} = 10**(-90/10);
	$errHash{'|'} = 10**(-91/10);
	$errHash{'}'} = 10**(-92/10);
	$errHash{'~'} = 10**(-93/10);


	#$errHash{'!'} = 1.00000000;
	#$errHash{'"'} = 0.79432823;
	#$errHash{'#'} = 0.63095734;
	#$errHash{'$'} = 0.50118723;
	#$errHash{'%'} = 0.39810717;
	#$errHash{'&'} = 0.31622777;
	#$errHash{'\''} = 0.25118864;
	#$errHash{'('} = 0.19952623;
	#$errHash{')'} = 0.15848932;
	#$errHash{'*'} = 0.12589254;
	#$errHash{'+'} = 0.10000000;
	#$errHash{','} = 0.07943282;
	#$errHash{'-'} = 0.06309573;
	#$errHash{'.'} = 0.05011872;
	#$errHash{'/'} = 0.03981072;
	#$errHash{'0'} = 0.03162278;
	#$errHash{'1'} = 0.02511886;
	#$errHash{'2'} = 0.01995262;
	#$errHash{'3'} = 0.01584893;
	#$errHash{'4'} = 0.01258925;
	#$errHash{'5'} = 0.01000000;
	#$errHash{'6'} = 0.00794328;
	#$errHash{'7'} = 0.00630957;
	#$errHash{'8'} = 0.00501187;
	#$errHash{'9'} = 0.00398107;
	#$errHash{':'} = 0.00316228;
	#$errHash{';'} = 0.00251189;
	#$errHash{'<'} = 0.00199526;
	#$errHash{'='} = 0.00158489;
	#$errHash{'>'} = 0.00125893;
	#$errHash{'?'} = 0.00100000;
	#$errHash{'@'} = 0.00079433;
	#$errHash{'A'} = 0.00063096;
	#$errHash{'B'} = 0.00050119;
	#$errHash{'C'} = 0.00039811;
	#$errHash{'D'} = 0.00031623;
	#$errHash{'E'} = 0.00025119;
	#$errHash{'F'} = 0.00019953;
	#$errHash{'G'} = 0.00015849;
	#$errHash{'H'} = 0.00012589;
	#$errHash{'I'} = 0.00010000;
	#$errHash{'J'} = 0.00007943;


	return \%errHash;

} # end sub get_err_hash_illumina

########

sub convert_hla_name { # start sub convert_hla_name

	my ($name) = @_;

	my ($new,$t);

	$name =~ s/:/_/g;

	$name =~ s/\*/_/g;

	if($name=~/A/ || $name =~ /a/){

		$name =~ s/[Aa]/hla_a/g;

	}

        if($name=~/B/ || $name =~ /b/){

                $name =~ s/[Bb]/hla_b/g;

        }

        if($name=~/C/ || $name =~ /c/){

                $name =~ s/[Cc]/hla_c/g;

        }

        if($name=~/HLA-A/){

                $name =~ s/HLA-A/hla_a/g;

        }

        if($name=~/HLA-B/){

                $name =~ s/HLA-B/hla_b/g;

        }

        if($name=~/HLA-C/){

                $name =~ s/HLA-C/hla_c/g;

        }


	return ($name);

} # end sub convert_hla_name

########

sub get_4digit { # start sub get_4digit

	my ($target) = @_;

	my ($h1,@f2,$i);

	@f2 = split(/_/,$target);

	for($i = 0; $i < 4; $i++){

		if($i > $#f2){

		last;

		}

		if($i == 0){

			$h1 = $f2[$i];

			next;

		}

		$h1 = $h1."_".$f2[$i];

	}

	return $h1;

} # end sub get_4digit


########

sub population_prior_likelihood_component { # start sub population_prior_likelihood_component

	my ($ref2,$target,$race,$scale) =  @_;

	my (@arr,%freqHash,$likeScore,$likeScoreLog,$f1,$h1,$likeScoreFreq,$likeScoreLogFreq,
		$likeScoreFreqLog,$likeScoreLogFreqLog,$likeScoreFreqZero,$likeScoreLogFreqZero,
		$likeScoreFreqLogZero,$likeScoreLogFreqLogZero,$logf1,@f2,$i);

	%freqHash = %{$ref2};

	@f2 = split(/_/,$target);

	for($i = 0; $i < 4; $i++){

		if($i > $#f2){

			last;

		}		

		if($i == 0){

			$h1 = $f2[$i];

			next;

		}	

		$h1 = $h1."_".$f2[$i];

	}

	#print "h1=$h1\t$target\t$race\n";

	$f1 = $freqHash{$h1}->{$race};

	return $f1;

	} # end sub population_prior_likelihood_component

########


sub population_prior_likelihood { # start sub population_prior_likelihood

	my ($ref1,$ref2,$target,$race,$scale) =  @_;

	my (@arr,%freqHash,$likeScore,$likeScoreLog,$f1,$h1,$likeScoreFreq,$likeScoreLogFreq,
		$likeScoreFreqLog,$likeScoreLogFreqLog,$likeScoreFreqZero,$likeScoreLogFreqZero,
		$likeScoreFreqLogZero,$likeScoreLogFreqLogZero,$logf1,@f2,$i);

#	$likeScoreFreqAdd,$likeScoreFreqAdd2,$likeScoreFreqZero,$likeScoreFreqZero2,$f1,$h1);

	@arr = @{$ref1};

	%freqHash = %{$ref2};

	$likeScoreLog = $arr[0];

	$likeScore = $arr[1];

	#$h1 = substr($target,0,9);
	
	@f2 = split(/_/,$target);

	for($i = 0; $i < 4; $i++){

		if($i > $#f2){

			last;

		}		

		if($i == 0){

			$h1 = $f2[$i];

			next;

		}	

		$h1 = $h1."_".$f2[$i];

	}


	$f1 = $freqHash{$h1}->{$race};

	#print "h1=$h1\tf1=$f1\t$race\n";

	if($f1 > 0){

		$logf1 = log($scale*$f1);

	} else {

		$logf1 = 0;

	}	

	$likeScoreFreq = $likeScore + $f1;

	$likeScoreLogFreq = $likeScoreLog + $f1;

	$likeScoreFreqLog = $likeScore + $logf1;

	$likeScoreLogFreqLog = $likeScoreLog + $logf1;

	if($f1 > 0){
		$likeScoreFreqZero = $likeScore + $f1;
	} else {
		$likeScoreFreqZero = 0;
	}

	if($f1 > 0){
		$likeScoreLogFreqZero = $likeScoreLog + $f1;
	} else {
		$likeScoreLogFreqZero = 0;
	}
 
	if($f1 > 0){
		$likeScoreFreqLogZero = $likeScore + $logf1;
	} else {
		$likeScoreFreqLogZero = 0;
	}

	if($f1 > 0){
		$likeScoreLogFreqLogZero = $likeScoreLog + $logf1;
	} else{
		$likeScoreLogFreqLogZero = 0;
	}
	

	@arr = ($likeScore,$likeScoreLog,$likeScoreFreq,$likeScoreLogFreq,
		$likeScoreFreqLog,$likeScoreLogFreqLog,$likeScoreFreqZero,$likeScoreLogFreqZero,
		$likeScoreFreqLogZero,$likeScoreLogFreqLogZero,$f1,$logf1);

	return (\@arr);

#	if($f1 > 0){
#
#		# 2 = log
#
#		$likeScoreFreqAdd = $likeScore + $f1;
#
#		$likeScoreFreqAdd2 = $likeScore2 + log($f1);
#
#		$likeScoreFreqZero = $likeScore + $f1;
#
#		$likeScoreFreqZero2 = $likeScore2 + log($f1);
#
#	} else{
#
#		$likeScoreFreqAdd = $likeScore;
#
#		$likeScoreFreqAdd2 = $likeScore2;
#
#		$likeScoreFreqZero = 0;
#
#		$likeScoreFreqZero2 = 0;
#
#	}
#
#	print "score\t$target\t$f1\t$likeScoreFreqZero2\t$likeScoreFreqAdd2\t$likeScore2\t$likeScoreFreqZero\t$likeScoreFreqAdd\t$likeScore\n";
#	#die;
#	
#	@arr = ($likeScoreFreqZero2,$likeScoreFreqAdd2,$likeScoreFreqZero,$likeScoreFreqAdd);
#
#	return (\@arr);



} # end sub population_prior_likelihood

########

sub process_read_pairs_second { # start sub process_read_pairs_second

	# function takes a list of read pairs and calculates the quality-based likelihood scores

	my ($segment,$ref2,$ref3,$target,$ref4,$ref5,$allele1) = @_;

	my (@pairs,$pair,$pId,$name,@pairNames,@pairIds,$length,$first_mate,$second_mate,$cigar1,
		$cigar2,$ref1,$matches1,$query1,$ref2,$matches2,$query2,$start1,$start2,$end1,$end2,
		$strand1,$strand2,@scores1,@scores2,$lik,$lik2,$likProd,$iScore,%pairTarget,$likeScore,
		$likeScore2,$likeProdScore,$h1,$f1,$finalScore,$finalScore_2,%freqHash,
		$likeScoreFreqZero,$likeScoreFreqZero2,$likeScoreFreqAdd,
		$likeScoreFreqAdd2,@arr,@pairNames,@pairIds,@filterPairs,$c1,%isizePval,$iScoreLog,
		$useFactor,$target_name1,$target_name2,$likAllele1,$likTarget,$factor);

	@pairs = $segment->features(-type => 'read_pair');

	@pairNames = ();

	@pairIds = ();

	%pairTarget = ();

	%freqHash = %{$ref2};

	%isizePval = %{$ref3};

	@filterPairs = @{$ref4};

	%pairTarget = %{$ref5};

	print "filterPairs = ",$#filterPairs,"\n";
	
	#print "freqHash 2:\n";

	#print Dumper(%freqHash);

	#die;

	$likeScore = 0;

	$likeScore2 = 0;
	
	$likeProdScore = 1; 


 	for $pair (@pairs) {

		$useFactor  = 0;

		$pId = $pair->primary_id;

		$name = $pair->name;

		$c1 = true { /^$name$/ } @filterPairs;


		if($c1 > 0){

			$likAllele1 = $pairTarget{$name}->{$allele1};

			$likTarget = $pairTarget{$name}->{$target};

			print "factor1\t$name\t$target\t$allele1\t$likTarget\t$likAllele1\n";

			if($likTarget > 0 || $likAllele1 > 0){

				$factor = $likTarget / ($likTarget + $likAllele1);		

				print "factor2\t$name\t$target\t$allele1\t$factor\n";

				$useFactor = 1;

			} else {
					
				$useFactor = 0;
			}

		}
	
		push @pairNames,$name;

		push @pairIds,$pId;

     		$length = $pair->length;   # insert length

	
         	($first_mate,$second_mate) = $pair->get_SeqFeatures;

		if(!defined($first_mate) || !defined($second_mate)){

			next;

		}

		$target_name1 = $first_mate->seq_id;
		$target_name2 = $second_mate->seq_id;

		#print "target_name\t$target_name1\t$target_name2\n";

		#die;

		if($target_name1 ne $target_name2){

			next;

		}

		$cigar1 = $first_mate->cigar;
		$cigar2 = $second_mate->cigar;
		($ref1,$matches1,$query1) = $first_mate->padded_alignment;
		($ref2,$matches2,$query2) = $second_mate->padded_alignment;
		$start1 = $first_mate->start;
	        $start2 = $second_mate->start;
		$end1 = $first_mate->end;
		$end2 = $second_mate->end;
		$strand1 = $first_mate->strand;
		$strand2 = $second_mate->strand;
		@scores1 = $first_mate->qscore;
		@scores2 = $second_mate->qscore;

                ($lik,$lik2,$likProd) = @{likelihood_score($ref1,$matches1,$query1,\@scores1,$ref2,$matches2,$query2,\@scores2,76)};	

		$iScore = $isizePval{$length};
		if($iScore > 0){
			$iScoreLog = log($iScore);
		} else {
			$iScoreLog = 0;
		}
		#$iScore = 0;

		#$pairTarget{$name}->{$target} = $lik2;

		#print "pair\t$name\t$target\t$lik\t$lik2\t$likProd\t$iScore\n";

		$lik = $lik + $iScore;

		$lik2 = $lik2 + $iScoreLog;

		if($useFactor == 1){
			
			$lik = $lik * $factor;

			$lik2 = $lik2 * $factor;

			$likProd = $likProd * $factor;

		}


		$likeScore = $likeScore + $lik;

		$likeScore2 = $likeScore2 + $lik2;
	
		#print "#### $name\t$start1\t$end1\t$strand1\t$start2\t$end2\t$strand2\t$length\t$lik\t$lik2\t$likProd\t$iScore ####\n\n";
		#print "$ref1\n$matches1\n$query1\n\n";
		#print "@scores1\n\n";
		#print "$ref2\n$matches2\n$query2\n\n";
		#print "@scores2\n\n";

	}

	print "likeScore = $likeScore2\t$likeScore\n";

	@arr = ($likeScore2,$likeScore);

	return (\@arr,\@pairNames,\@pairIds);


} # end sub process_read_pairs_second


########

sub process_read_pairs { # start sub process_read_pairs

	# function takes a list of read pairs and calculates the quality-based likelihood scores

	my ($segment,$ref2,$ref3,$target) = @_;

	my (@pairs,$pair,$pId,$name,@pairNames,@pairIds,$length,$first_mate,$second_mate,$cigar1,
		$cigar2,$ref1,$matches1,$query1,$matches2,$query2,$start1,$start2,$end1,$end2,
		$strand1,$strand2,@scores1,@scores2,$lik,$lik2,$likProd,$iScore,%pairTarget,$likeScore,
		$likeScore2,$likeProdScore,$h1,$f1,$finalScore,$finalScore_2,%freqHash,
		$likeScoreFreqZero,$likeScoreFreqZero2,$likeScoreFreqAdd,
		$likeScoreFreqAdd2,@arr,@pairNames,@pairIds,%isizePval,$iScoreLog,$target_name1,
		$target_name2,$useFactor);

	@pairs = $segment->features(-type => 'read_pair');

	@pairNames = ();

	@pairIds = ();

	%pairTarget = ();

	%freqHash = %{$ref2};

	%isizePval = %{$ref3};

	#print "freqHash 2:\n";

	#print Dumper(%freqHash);

	#die;

	$likeScore = 0;

	$likeScore2 = 0;
	
	$likeProdScore = 1; 


 	for $pair (@pairs) {

		$useFactor  = 0;

		$pId = $pair->primary_id;

		$name = $pair->name;

		push @pairNames,$name;

		push @pairIds,$pId;

     		$length = $pair->length;   # insert length	

         	($first_mate,$second_mate) = $pair->get_SeqFeatures;

		if(!defined($first_mate) || !defined($second_mate)){

			next;

		}

		#print Dumper($first_mate);
		
		$target_name1 = $first_mate->seq_id;
		$target_name2 = $second_mate->seq_id;

		#print "target_name\t$target_name1\t$target_name2\n";

		#die;

		print "pair\t$name\t$target\n";

		if($target_name1 ne $target_name2){

			next;

		}

		$cigar1 = $first_mate->cigar;
		$cigar2 = $second_mate->cigar;
		($ref1,$matches1,$query1) = $first_mate->padded_alignment;
		($ref2,$matches2,$query2) = $second_mate->padded_alignment;
		$start1 = $first_mate->start;
	        $start2 = $second_mate->start;
		$end1 = $first_mate->end;
		$end2 = $second_mate->end;
		$strand1 = $first_mate->strand;
		$strand2 = $second_mate->strand;
		@scores1 = $first_mate->qscore;
		@scores2 = $second_mate->qscore;

                ($lik,$lik2,$likProd) = @{likelihood_score($ref1,$matches1,$query1,\@scores1,$ref2,$matches2,$query2,\@scores2,76)};	

		$iScore = $isizePval{$length};
		if($iScore > 0){
			$iScoreLog = log($iScore);
		} else {
			$iScoreLog = 0;
		}
		#$iScore = 0;

		#$pairTarget{$name}->{$target} = $lik2;

		if($lik > 0){

			print "pair\t$name\t$target\t$lik\t$lik2\t$likProd\t$iScore\t$length\n";

			#$lik = $lik + $iScore;

			$lik2 = $lik2 + $iScoreLog;

			$likeScore = $likeScore + $lik;

			$likeScore2 = $likeScore2 + $lik2;

			$pairTarget{$name}->{$target} = $lik;


		}
	
		#print "#### $name\t$start1\t$end1\t$strand1\t$start2\t$end2\t$strand2\t$length\t$lik\t$lik2\t$likProd\t$iScore ####\n\n";
		#print "$ref1\n$matches1\n$query1\n\n";
		#print "@scores1\n\n";
		#print "$ref2\n$matches2\n$query2\n\n";
		#print "@scores2\n\n";

	}

	#print "likeScore = $likeScore2\t$likeScore\n";

	@arr = ($likeScore2,$likeScore);

	return (\@arr,\@pairNames,\@pairIds,\%pairTarget);


} # end sub process_read_pairs

########

sub insert_score { #start sub insert_score

	my ($length,$median,$ref) = @_;

	my (@isizesMedianSub,$i,$num,$iScore,$diff);

	@isizesMedianSub = @{$ref};

	$num = 0;

	$diff = abs($length-$median);

	for($i = 0; $i <= $#isizesMedianSub; $i++){

		if($diff <= abs($isizesMedianSub[$i])){

			$num = $num + 1;

		}			
	
	}


	$iScore = $num / ($#isizesMedianSub + 1);


	return $iScore;


} # end sub insert_score

########

sub process_insert_size { # start sub process_insert_size

	# function takes the insert file generated by picard CollectInsertSizeMetrics and returns a hash
	# containing the probability of different insert lengths

	my ($iFile) = @_;

	my ($line,$i,@f,$isize,@isizes,%isizePval,$num,$p,$j,$totalCount,$count);

	@isizes = ();

	%isizePval = ();

	open FILE, $iFile || die  "Cannot open $iFile\n";

	while(($line=<FILE>) !~/^insert_size/){}
	
	$totalCount  = 0;

	while($line=<FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		$isize = $f[0];

		$count = $f[1];

		if($isize !~ /\d/){

			next;

		}	

		$totalCount = $totalCount + $count;

		$isizePval{$isize} = $count;		
		
		push @isizes, $isize;

	}

	close FILE;


	for($i = 0; $i <= $#isizes; $i++){

		$isize = $isizes[$i];

		$isizePval{$isize} = $isizePval{$isize}/$totalCount;

	}

	return \%isizePval;


} # end sub process_insert_size

########


sub likelihood_score{ # start sub likelihood_score

	# function takes a read pair alignment in terms of 2 sets of ref,matches,query and the 
	# corresponding quality scores as input. It then calculates the data-supported likelihood.
	# For now, any read pair which has any indel anywhere is discarded. In the future it would be
	# good to have an indel model

	my ($ref1,$matches1,$query1,$refS1,$ref2,$matches2,$query2,$refS2,$length) = @_;

	my (@scores1,@scores2,$len,@f1,@f2,$i,$s1,$s2,$e1,$e2,$m1,$m2,
		$lik,@refArr1,@refArr2,@qArr1,@qArr2,$c1,$c2,$likProd,@arr,$lik2,$scores2);

	@scores1 = @{$refS1};

	@scores2 = @{$refS2};

	@f1 = split(//,$matches1);

	@f2 = split(//,$matches2);

	@refArr1 = split(//,$ref1);

	@refArr2 = split(//,$ref2);

	@qArr1 = split(//,$query1);

	@qArr2 = split(//,$query2);

	$c1 = true { /\|/ } @f1;

	$c2 = true { /\|/ } @f2;

	#print "c1=$c1\tc2=$c2\n";

	if(($c1 != $length) || ($c2 != $length)){

		#return 0;

	}
	
	if($ref1=~/-/ || $query1=~/-/ || $ref2=~/-/ || $query2=~/-/){
	#if($ref1=~/\w-+\w/ || $query1=~/\w-+\w/ || $ref2=~/\w-+\w/ || $query2=~/\w-+\w/){

		@arr = (0,0,0);

		return \@arr;

	}

	$lik = 0;

	$lik2 = 0;
	
	$likProd = 1;

	for($i = 0; $i <= $#scores1; $i++){

		$s1 = $scores1[$i];

		$e1 = 10**(-$s1/10);
	
		if($refArr1[$i] ne "-"){
		#if($refArr1[$i] ne "-" && $qArr1[$i] ne "-"){

			if($qArr1[$i] eq $refArr1[$i]){

				$lik = $lik + 1 - $e1;
				$lik2 = $lik2 + log(1 - $e1);

				$likProd = $likProd*(1-$e1);

			} else {

				$lik = $lik + $e1/3;
				#$lik = $lik - (1 - $e1);
				$lik2 = $lik2 + log($e1/3);

				$likProd = $likProd*($e1/3);
	
			}

		} else {

			#$lik = $lik - (1 - $e1);

		}

		#print "s1=$s1\t$e1\t$s2\t$e2\t$qArr1[$i]\t$refArr1[$i]\t$lik\n";

	}

	#print "lik=$lik\n";

	for($i = 0; $i <= $#scores2; $i++){

		$s2 = $scores2[$i];

		$e2 = 10**(-$s2/10);	
	
		if($refArr2[$i] ne "-"){
		#if($refArr2[$i] ne "-" && $qArr2[$i] ne "-"){
	
			if($qArr2[$i] eq $refArr2[$i]){

				$lik = $lik + 1 - $e2;
				$lik2 = $lik2 + log(1 - $e2);

				$likProd = $likProd*(1-$e2);

			} else {
		
				$lik = $lik + $e2/3;
				#$lik = $lik - (1 - $e1);
				$lik2 = $lik2 + log($e2/3);

				$likProd = $likProd*($e2/3);

			}

		} else {

			#$lik = $lik - (1 - $e2);

		}
		
		#print "s1=$s1\t$e1\t$s2\t$e2\t$qArr1[$i]\t$refArr1[$i]\t$lik\n";


	}

	#print "lik=$lik\n";

	@arr = ($lik,$lik2,$likProd);
	
	#print "fn: likeProd=$likProd\n";
	
	return \@arr;

} # end sub likelihood_score


########

sub likelihood_score2{ # start sub likelihood_score2

	# function takes a read pair alignment in terms of 2 sets of ref,matches,query and the 
	# corresponding quality scores as input. It then calculates the data-supported likelihood.
	# For now, any read pair which has any indel anywhere is discarded. In the future it would be
	# good to have an indel model

	my ($ref1,$matches1,$query1,$refS1,$ref2,$matches2,$query2,$refS2,$length) = @_;

	my (@scores1,@scores2,$len,@f1,@f2,$i,$s1,$s2,$e1,$e2,$m1,$m2,
		$lik,@refArr1,@refArr2,@qArr1,@qArr2,$c1,$c2,$likProd,@arr,$lik2);

	@scores1 = @{$refS1};

	@scores2 = @{$refS2};

	@f1 = split(//,$matches1);

	@f2 = split(//,$matches2);

	@refArr1 = split(//,$ref1);

	@refArr2 = split(//,$ref2);

	@qArr1 = split(//,$query1);

	@qArr2 = split(//,$query2);

	$c1 = true { /\|/ } @f1;

	$c2 = true { /\|/ } @f2;

	#print "c1=$c1\tc2=$c2\n";

	if(($c1 != $length) || ($c2 != $length)){

		#return 0;

	}
	
	if($ref1=~/-/ || $query1=~/-/ || $ref2=~/-/ || $query2=~/-/){
	#if($ref1=~/\w-+\w/ || $query1=~/\w-+\w/ || $ref2=~/\w-+\w/ || $query2=~/\w-+\w/){

		@arr = (0,0,0);

		return \@arr;

	}

	$lik = 0;

	$lik2 = 0;
	
	$likProd = 1;

	for($i = 0; $i < $length; $i++){

		$s1 = $scores1[$i];

		$s2 = $scores2[$i];

		$e1 = 10**(-$s1/10);

		$e2 = 10**(-$s2/10);	

	
		#print "s1=$s1\t$e1\t$s2\t$e2\t$qArr1[$i]\t$refArr1[$i]\n";

		if($refArr1[$i] ne "-"){
		#if($refArr1[$i] ne "-" && $qArr1[$i] ne "-"){

			if($qArr1[$i] eq $refArr1[$i]){

				$lik = $lik + 1 - $e1;
				$lik2 = $lik2 + log(1 - $e1);

				$likProd = $likProd*(1-$e1);

			} else {

				$lik = $lik + $e1/3;
				#$lik = $lik - (1 - $e1);
				$lik2 = $lik2 + log($e1/3);

				$likProd = $likProd*($e1/3);
	
			}

		} else {

			#$lik = $lik - (1 - $e1);

		}

		if($refArr2[$i] ne "-"){
		#if($refArr2[$i] ne "-" && $qArr2[$i] ne "-"){
	
			if($qArr2[$i] eq $refArr2[$i]){

				$lik = $lik + 1 - $e2;
				$lik2 = $lik2 + log(1 - $e2);

				$likProd = $likProd*(1-$e2);

			} else {
		
				$lik = $lik + $e2/3;
				#$lik = $lik - (1 - $e1);
				$lik2 = $lik2 + log($e2/3);

				$likProd = $likProd*($e2/3);

			}

		} else {

			#$lik = $lik - (1 - $e2);

		}


	}

	@arr = ($lik,$lik2,$likProd);
	
	#print "fn: likeProd=$likProd\n";
	
	return \@arr;

} # end sub likelihood_score2

####

sub get_race { #start sub get_race

	my ($file) = @_;

	my (%raceHash, $line, @f);

	%raceHash  = ();

	open FILE, $file || die "get_race: Cannot open $file\n";

	while($line = <FILE>){

		chomp($line);
			
		@f = split(/\t/,$line);

		$raceHash{$f[0]} = $f[1];

	}

	close FILE;

	return \%raceHash;


} # end sub get_race

####

sub get_hla_frequencies { # start sub get_hla_frequencies

	my ($file,$FFILE) = @_;

	my ($line,@f,$line,$a,%hash);

	open FILE, $file || die "get_hla_frequencies: Cannot open $file\n";

	$line = <FILE>;

	while($line = <FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		$a = $f[0];

		if($FFILE eq "f2" || $FFILE eq "f1"){

			$a =~ s/\*//;

			$a =~ tr/A-Z/a-z/;

			$a = "hla_".$a;

		}elsif($FFILE eq "f3" || $FFILE eq "f4"){

			# do nothing

		} else {

			die "get_hla_frequencies: FFILE = $FFILE is invalid\n";
			
		}

		$hash{$a}->{'Caucasian'} = $f[1];

		$hash{$a}->{'Black'} = $f[2];

		$hash{$a}->{'Asian'} = $f[3];

	}


	close FILE;

	return \%hash;

} # end sub get_hla_frequencies


########

sub median {
 
    my @array = sort {$a <=> $b} @_;

    my $median;

    if ($#array % 2 == 0) {
        $median = $array[($#array / 2)];
    }else {
       $median = $array[int($#array / 2)] + (($array[int($#array / 2) + 1] - $array[int($#array / 2)]) / 2);
    }

    return $median;
}

########

sub get_winning_pair { # start sub get_winning_pair

	my ($ref1,$ref2,$ref3,$ref4) = @_;

	my (@targets,@pairs,%pairTarget,$i,$j,$t1,$t2,$k,$pair,%freqHash,$l1,$l2,$targetPairScore,$f1,$f2,
	$totalTargetPairScore,$h1,$h2,$t);

	@targets = @{$ref1};

	@pairs = @{$ref2};

	%pairTarget = %{$ref3};

	%freqHash = %{$ref4};

	for($i = 0; $i < $#targets; $i++){

		$t1 = $targets[$i];

		$h1 = substr($t1,0,9);

		$f1 = $freqHash{$h1}->{'white'};

		$targetPairScore = 0;

		for($j = ($i+1); $j <= $#targets; $j++){

			$t2 = $targets[$j];
		
			$h2 = substr($t2,0,9);

			$f2 = $freqHash{$h2}->{'white'};

			for($k = 0; $k <= $#pairs; $k++){

				$pair = $pairs[$k];

				$l1 = $pairTarget{$pair}->{$t1};	

				$l2 = $pairTarget{$pair}->{$t2};
		
				print "t1=$t1\t$t2\t$pair\tl1=$l1\tl2=$l2\tf1=$f1\tf2=$f2\n";

				$t = 0.5*$l1 + 0.5*$l2;
		
				if($t!=0){

					$targetPairScore = $targetPairScore + log($t);

				}
				

			}

			$totalTargetPairScore = $targetPairScore;

			if($f1!=0){

				$totalTargetPairScore = $totalTargetPairScore + log($f1);

			}

			if($f2!=0){

				$totalTargetPairScore = $totalTargetPairScore + log($f2);
			
			}

			print "totalTargetPairScore\t$t1\t$t2\t$targetPairScore\t$totalTargetPairScore\n";

		}

	}	
	


} #end sub get_winning_pair

########

sub process_insert_size_median { # start sub process_insert_size_median

	# function takes the insert file generated by picard CollectInsertSizeMetrics and calculates
	# the median insert size. It also returns a arrayof all median subtracted insert size lengths.

	my ($iFile) = @_;

	my ($line,$i,@f,@isizes,$median,$isize,$count,@isizesMedianSub,$t,$max,%isizePval,$found,
		$num,$p,$j);

	@isizes = ();

	open FILE, $iFile;

	while(($line=<FILE>) !~/^insert_size/){}
	
	while($line=<FILE>){

		chomp($line);

		@f = split(/\t/,$line);

		$isize = $f[0];

		$count = $f[1];

		for($i = 0; $i < $count; $i++){

			push @isizes,$isize;

		}

	}
	
	close FILE;
	

	#print "isizes=",$isizes[$#isizes],"\n";

	$median = median(@isizes);

	#print "median=$median\n";

	@isizesMedianSub = ();

	for($i = 0; $i <= $#isizes; $i++){

		$t = abs($isizes[$i] - $median);

		push @isizesMedianSub, $t;
	}


	$max = max @isizesMedianSub;

	@isizesMedianSub = sort { $a <=> $b } @isizesMedianSub;

	for($i = 0; $i <= $max; $i++){

		#@big_numbers = grep { $_ > 4 } @numbers;
	
		$num = 0;

		$found = 0;

		for($j = 0; $j <= $#isizesMedianSub; $j++){

			#if($


		}		

	}
		

	return ($median,\@isizesMedianSub);

} # end sub process_insert_size_median

########

###########################################

#function returns 1 if the name is found in the array;0 otherwise

#int found_arr(name,\array)

sub found_arr{ #START FOUND_ARR

	my ($fa_name,$fa_x1)=@_;

	my (@fa_name_arr,$fa_i,$fa_found);

	@fa_name_arr=@{$fa_x1};

	$fa_found=0;

	for($fa_i=0;$fa_i<=$#fa_name_arr;$fa_i++){ #START FA FOR 1

		if($fa_name eq $fa_name_arr[$fa_i]){

			$fa_found=1;

			return $fa_found;

		}

	} #END FA FOR 1

	return $fa_found;

} #END FOUND_ARR



1;
