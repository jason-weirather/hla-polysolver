#!/usr/bin/perl

# usage: $PSHOME/getWinnerScores.pl counts1.R0k6 $PSHOME


$file = $ARGV[0];
$PSHOME = $ARGV[1];

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;



$outFile = $file.".top2";
$sortedFile = $file.".sorted";

system("rm -f $sortedFile");
system("sort -k2,2rn $file > $sortedFile");

$freqFile = $PSHOME."/data/HLA_FREQ.txt";
$scale=exp(23);

print "sortedFile=$sortedFile\toutFile=$outFile\n";
`rm -f $outFile`;

# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};

open FILE, $sortedFile || die "Cannot open $sortedFile\n";
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";

# sort alleleScores by value

$aCount = 0;
$bCount = 0;
$cCount = 0;

@alleles = ();
while($line = <FILE>){
	chomp($line);
	@f = split(/\t/,$line);
	$allele = $f[0];
	$alleleScore = $f[1];
	#print "a=",$allele,"\t",$alleleScore,"\n";
	$allele4digit = get_4digit($allele);
	$freqWhite = population_prior_likelihood_component(\%freqHash,$allele,"Caucasian",$scale);
	$freqBlack = population_prior_likelihood_component(\%freqHash,$allele,"Black",$scale);
	$freqAsian = population_prior_likelihood_component(\%freqHash,$allele,"Asian",$scale);
	#if($freqWhite == 0 & $freqBlack == 0 & $freqAsian == 0){
	#	next;
	#}
	if($allele =~/hla_a/){
		if(found_arr($allele4digit,\@alleles)){
			next;
		}
		if($aCount == 0){
			$aFirst = $allele;
			$aFirstScore = $alleleScore;
			$aCount++;
			push @alleles,$allele4digit;
			print "aFirst=$aFirst\t$aFirstScore\n";
		}elsif($aCount == 1){
			$aSecond = $allele;
			$aSecondScore = $alleleScore;
			$aCount++;
                        push @alleles,$allele4digit;
			print "aSecond=$aSecond\t$aSecondScore\n";
		}
	}
        if($allele =~/hla_b/){
                if(found_arr($allele4digit,\@alleles)){
                        next;
                }
                if($bCount == 0){
                        $bFirst = $allele;
                        $bFirstScore = $alleleScore;
                        $bCount++;
                        push @alleles,$allele4digit;
                        print "bFirst=$bFirst\t$bFirstScore\n";
                }elsif($bCount == 1){
                        $bSecond = $allele;
                        $bSecondScore = $alleleScore;
                        $bCount++;
                        push @alleles,$allele4digit;
			print "bSecond=$bSecond\t$bSecondScore\n";
                }
        }
        if($allele =~/hla_c/){
                if(found_arr($allele4digit,\@alleles)){
                        next;
                }
                if($cCount == 0){
                        $cFirst = $allele;
                        $cFirstScore = $alleleScore;
                        $cCount++;
                        push @alleles,$allele4digit;
			print "cFirst=$cFirst\t$cFirstScore\n";
                }elsif($cCount == 1){
                        $cSecond = $allele;
                        $cSecondScore = $alleleScore;
                        $cCount++;
                        push @alleles,$allele4digit;
			print "cSecond=$cSecond\t$cSecondScore\n";
                }
        }



	#if($aCount == 2  & $bCount == 2 & $cCount == 2){
	#	last;
	#}

}

#$aDiff = $aFirstScore - $aSecondScore;
#$bDiff = $bFirstScore - $bSecondScore;
#$cDiff = $cFirstScore - $cSecondScore;
#$aNorm = 2*$aDiff/($aFirstScore+$aSecondScore);
#$bNorm = 2*$bDiff/($bFirstScore+$bSecondScore);
#$cNorm = 2*$cDiff/($cFirstScore+$cSecondScore);
#$aNorm2 = $aDiff/($aFirstScore);
#$bNorm2 = $bDiff/($bFirstScore);
#$cNorm2 = $cDiff/($cFirstScore);
($aDiff,$aNorm,$aNorm2) = @{get_stats($aFirstScore,$aSecondScore)};
($bDiff,$bNorm,$bNorm2) = @{get_stats($bFirstScore,$bSecondScore)};
($cDiff,$cNorm,$cNorm2) = @{get_stats($cFirstScore,$cSecondScore)};

print OUTFILE "$aFirst\t$aFirstScore\t$aSecond\t$aSecondScore\t$aDiff\t$aNorm\t$aNorm2\n";
print OUTFILE "$bFirst\t$bFirstScore\t$bSecond\t$bSecondScore\t$bDiff\t$bNorm\t$bNorm2\n";
print OUTFILE "$cFirst\t$cFirstScore\t$cSecond\t$cSecondScore\t$cDiff\t$cNorm\t$cNorm2\n";

close FILE;
close OUTFILE;

########

sub get_stats { # start sub get_stats

	my ($firstScore,$secondScore) = @_;
	my ($diff,$norm,$norm2,@arr);
	$diff = "NA";
	$norm = "NA";
	$norm2 = "NA";
	if($firstScore=~/^$/ & $secondScore=~/^$/){
		$diff="NA";
	} else{
		$diff = $firstScore - $secondScore; 
	}
	if(($firstScore+$secondScore) > 0){
		$norm = 2*$diff/($firstScore+$secondScore); 
	}
	if($firstScore > 0){
		$norm2 = $diff/$firstScore;
	}

	@arr = ($diff,$norm,$norm2);

	return \@arr;

} # end sub get_stats
