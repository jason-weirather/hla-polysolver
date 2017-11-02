#!/usr/bin/perl


use Bio::DB::Sam;
use Dumpvalue;
use Data::Dumper;
use List::MoreUtils 'true';
use List::MoreUtils 'indexes';
use List::MoreUtils qw/ uniq /;
use List::Util qw/max/;

$bam = $ARGV[0];
$PSHOME = $ARGV[1];
$SAMTOOLS = $ARGV[2];
$iFile = $ARGV[3];
$race = $ARGV[4];
$allele = $ARGV[5];
$includeFreq = $ARGV[6];
$outDir = $ARGV[7];


$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;


if(!($includeFreq == 0 || $includeFreq == 1)){
	die "includeFreq = $includeFreq is invalid\n";
}

#if(!($race eq "Caucasian" || $race eq "Black" || $race eq "Asian") && $includeFreq == 1){
#	die "With race = $race includeFreq must be 0\n";
#}

$scale = exp(23);
$freqFile = $PSHOME."/"."data/HLA_FREQ.txt"; 

if($includeFreq == 1){ # will always be 0 if race is other than Caucasian, Black, Asian
	$outFile = $outDir."/".$allele.".lik1";
} else {
	$outFile = $outDir."/".$allele.".nofreq.lik1";
}

open OUTFILE, ">$outFile";

# process insert size file


if($iFile ne "0"){
	%isizePval = %{process_insert_size($iFile)};
}

# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};

# get population frequency component

# if the race is Caucasian, Black or Asian:
	# includeFreq = 1 will give 97% accuracy
	# includeFreq = 0 will give 89%
# if race is Unknown
	# includeFreq = 1 will give 96%
	# includeFreq = 0 will give 89%

$freqComponent = 0;
$freqComponentLog = 0;

if($race eq "Caucasian" || $race eq "Black" || $race eq "Asian"){
	$freqComponent = population_prior_likelihood_component(\%freqHash,$allele,$race,$scale);
	print "race=$race\tallele=$allele\tfreq=$freqComponent\n";
	if($freqComponent > 0){
		$freqComponentLog = log($freqComponent);
	} else {
		$freqComponentLog = 0;
	}
	if($includeFreq == 1){
		if($freqComponent == 0){
			print OUTFILE "lik1\t0\t0\t0\n";
			close OUTFILE;
			exit 0;
			#die;
		}
	}
} else {
	$freqCaucasian = population_prior_likelihood_component(\%freqHash,$allele,"Caucasian",$scale);
        $freqBlack = population_prior_likelihood_component(\%freqHash,$allele,"Black",$scale);
        $freqAsian = population_prior_likelihood_component(\%freqHash,$allele,"Asian",$scale);
	if($includeFreq == 1) {
		if($freqCaucasian == 0 & $freqBlack == 0 & $freqAsian == 0){
        		print OUTFILE "lik1\t0\t0\t0\n";
			close OUTFILE;
			exit 0;
			#die;
		}
	}
}	 


# get error hash

%errHash = %{get_err_hash_illumina()};

# get score hash

($ref1,$ref2) = get_score_hash_log(\%errHash,$scale);
%scoreHash = %{$ref1};
%scoreHashLog = %{$ref2};


# calculate likelihood

@lines = `$SAMTOOLS view $bam`;

$likScore = 0;

$asScore = 0;

for($i=0;$i<=$#lines;$i++){
	$line1 = $lines[$i];
	$line2 = $lines[$i+1];
	@f1 = split(/\t/,$line1);
	@f2 = split(/\t/,$line2);
	if($f1[0] ne $f2[0]){
		next;
	}
	$i++;
	if($f1[2] ne $f2[2] || $f1[5]=~/I/ || $f1[5]=~/D/ || $f1[5]=~/\*/ || $f2[5]=~/I/ || $f2[5]=~/D/ || $f2[5]=~/\*/ || $f1[6] ne "=" || $f2[6] ne "="){
		next;
	}	
	$allele = $f1[2];
	$q1 = $f1[10];
	$q2 = $f2[10];
	$mapq1 = $f1[4];
	$mapq2 = $f2[4];
	$_=$line1;
	/(MD:Z:\w+)/;
	$md1 = $1;
	$_=$line2;
	/(MD:Z:\w+)/;
	$md2 = $1;
	$_ = $line1;
	/AS:i:(\w+)/;
	$as1 = $1;
	$_ = $line2;
	/AS:i:(\w+)/;
	$as2 = $1;
	$insert = $f1[8];
	$md1String = md_to_match_string($md1);
	$md2String = md_to_match_string($md2);
	$matchCount = ($md1String =~ tr/M//) + ($md2String =~ tr/M//);
	$mismatchCount = ($md1String =~ tr/m//) + ($md2String =~ tr/m//);
	$l1 = 0;
	$l2 = 0;
	$l1Log = 0;
	$l2Log = 0;
	$iScoreLog = 0;
	($l1,$l1Log,$l1Match,$l1LogMatch,$l1Mismatch,$l1LogMismatch) = @{get_likelihood_log(\%scoreHash,\%scoreHashLog,$q1,$md1String)};
	($l2,$l2Log,$l2Match,$l2LogMatch,$l2Mismatch,$l2LogMismatch) = @{get_likelihood_log(\%scoreHash,\%scoreHashLog,$q2,$md2String)};
	if($iFile ne "0"){
		$iScore = $isizePval{abs($insert)};
		if($iScore > 0){
			$iScoreLog = log($iScore); # not scaling the insert size p-val
		}
	}
	print "iScoreLog=$iScoreLog\n";
	$likLog = $l1Log + $l2Log;
	$likLogIscore = $likLog + $iScoreLog;
	print OUTFILE "$f1[0]\t$likLogIscore\t$likLog\t$iScoreLog\n";
	$likLogIscoreTotal = $likLogIscoreTotal + $likLogIscore;
}


if($includeFreq == 1){

	$likLogIscoreTotalFreq = $likLogIscoreTotal + $freqComponentLog;


} else {

	$likLogIscoreTotalFreq = $likLogIscoreTotal;

}

print OUTFILE "lik1\t$likLogIscoreTotalFreq\t$likLogIscoreTotal\t$freqComponentLog\n";


close OUTFILE;

