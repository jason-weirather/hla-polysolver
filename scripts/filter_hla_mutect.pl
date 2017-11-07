#!/usr/bin/env perl

#usage: ./filter_hla_mutect.pl HNSC-TCGA-BA-6873/HNSC-TCGA-BA-6873.mutect.unfiltered.annotated HNSC-TCGA-BA-6873 HNSC-TCGA-BA-6873 0 /cga/wu/sachet/hla/hla_caller/capture/polysolver_070914


use List::MoreUtils qw/ uniq /;

$annotFile = $ARGV[0];
$indiv = $ARGV[1];
$dir = $ARGV[2];
$N_ALT_THRESHOLD = $ARGV[3];
$PSHOME = $ARGV[4];

$outFile = $dir."/".$indiv.".mutect.filtered.annotated";
$chooseManuallyFile = $dir."/".$indiv.".mutect.ambiguous.annotated";

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;

$scale = exp(23);
$freqFile = $ENV{'DATA_DIR'}."/HLA_FREQ.txt";


# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};


@lines = ();
open AFILE, $annotFile;
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";


# get the column numbers of t_ref_count, t_alt_count, n_ref_count and n_alt_count

$header = <AFILE>;
chomp($header);
@f = split(/\t/,$header);
for($i = 0; $i <= $#f; $i++){
	$colName = $f[$i];
	if($colName eq "t_ref_count"){
		$trc = $i + 1;
	}
        if($colName eq "t_alt_count"){
                $tac = $i + 1;
        }
        if($colName eq "n_ref_count"){
                $nrc = $i + 1;
        }
        if($colName eq "n_alt_count"){
                $nac = $i + 1;
        }

}


print OUTFILE "$header\n";

while($line = <AFILE>){
	chomp($line);
	if($line=~/^#/ || $line !~ /KEEP/ || $line =~ /intron/){
		next;
	}
	@f = split(/\t/,$line);
	$sample = $f[0];
	$allele = $f[1];
	$pos = $f[2];
	$tRefCount = $f[$trc - 1];
        $tAltCount = $f[$tac - 1];
        $nRefCount = $f[$nrc - 1];
        $nAltCount = $f[$nac - 1];
	$peptideChange = $f[-1];
	$key = $sample.":".$allele;
	#print "$line\nsample=$sample\ttRefCount=$tRefCount\t$tAltCount\t$nRefCount\t$nAltCount\n";
	if($nAltCount > $N_ALT_THRESHOLD){
		next;
	}
	if($tAltCount < 2){
		next;
	}
        $keyCountHash{$key} =  $keyCountHash{$key} + 1;
	push @lines, $line;
}

# remove all lines where keyCount > 1

@lines2 = ();
%sampleChangeCountHash = (); # count of all sample:peptideChange occurrence for a particular sample and peptideChange
%sampleChangeLinesHash = (); # count of all lines corresponding to a particular sample:peptideChange
%sampleChangesHash = (); # list of all peptideChanges (including duplicates) for a sample
%sampleCountHash = (); # count of total peptideChanges (including duplicates) for a sample
@samples = ();
for($i = 0; $i <= $#lines; $i++){
	$line = $lines[$i];
	#print "LINE=$line\n";
	#next;
	@f = split(/\t/,$line);
        $sample = $f[0];
        $allele = $f[1];
        $pos = $f[2];
	$peptideChange = $f[-1];
	#print "PEPTIDECHANGE=$peptideChange\n";
	#next;
        $key = $sample.":".$allele;
	if($keyCountHash{$key} > 1){
		next;
	}
	push @samples, $sample;
	$sampleCountHash{$sample}++;
        $sampleChangeCountHash{$sample}->{$peptideChange}++;
	if($sampleCountHash{$sample} == 1){
		$sampleChangesHash{$sample} =  [ $peptideChange ];
	} else {
                @arr =  @{$sampleChangesHash{$sample}};
                push @arr, $peptideChange;
                $sampleChangesHash{$sample} =  [ @arr ];
	}
        if($sampleChangeCountHash{$sample}->{$peptideChange} == 1){
		$sampleChangeLinesHash{$sample}->{$peptideChange} =  [ $line ] ; 

        } else {
		@arr =  @{$sampleChangeLinesHash{$sample}->{$peptideChange}};
		push @arr, $line;
                $sampleChangeLinesHash{$sample}->{$peptideChange} = [ @arr ];

	}
	push @lines2, $line;
}


# choose only one allele 

open CMFILE,">$chooseManuallyFile" || die "Cannot open $chooseManuallyFile\n";

print CMFILE "$header\n"; 

@samples = uniq @samples;

for($i = 0; $i <= $#samples; $i++){
	$sample = $samples[$i];
	@changes = @{$sampleChangesHash{$sample}};
	@changes = uniq @changes;
	#next;
	for($j = 0; $j <= $#changes; $j++){
		$change = $changes[$j];
		#print "CHANGE=$change\t$sampleChangeCountHash{$sample}->{$change}\n";
		#next;
		if($sampleChangeCountHash{$sample}->{$change} == 1){
			@arr = @{$sampleChangeLinesHash{$sample}->{$change}};
			$line = $arr[0];
			print OUTFILE "$line\n";
		} else {
			@arr = @{$sampleChangeLinesHash{$sample}->{$change}};
			#print "#### choose 1 ####\n";
			for($k = 0; $k <= $#arr; $k++){
				$line = $arr[$k];
				@f = split(/\t/,$line);
				$tRefCount = $f[$trc - 1];
			        $tAltCount = $f[$tac - 1];
       				$nRefCount = $f[$nrc - 1];
        			$nAltCount = $f[$nac - 1];
				if($k == 0){
					$chosenLine = $line;
					$tAltCountMax = $tAltCount;	
				}
				if($tAltCount > $tAltCountMax){
					$chosenLine = $line;
					$tAltCountMax = $tAltCount;
				}

			}
			@lines2 = ();
			for($k = 0; $k <= $#arr; $k++){
				$line = $arr[$k];
                                @f = split(/\t/,$line);
                                $tRefCount = $f[$trc - 1];
                                $tAltCount = $f[$tac - 1];
                                $nRefCount = $f[$nrc - 1];
                                $nAltCount = $f[$nac - 1];
				if($tAltCount == $tAltCountMax){
					push @lines2, $line;
				}
			}
			if( ($#lines2+1) == 1){
				print OUTFILE "$chosenLine\n";
			} else{
				#print "can't choose\n";
				print OUTFILE "$chosenLine\n";
				for($k1 = 0; $k1 <= $#arr; $k1++){
					$line1 = $arr[$k1];
					print CMFILE "$line1\n";	
				}
			}
			#print "########\n";

		}
	}
	
}
close CMFILE;

close AFILE;
close OUTFILE;
