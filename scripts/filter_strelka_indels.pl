#!/usr/bin/env perl

#usage: ./filter_strelka_indels.pl out1.indels.071614.all.txt /cga/wu/sachet/hla/hla_caller/capture/polysolver_070914 choose_manually.strelka.txt > out1.indels.071614.all.filt.txt

use List::MoreUtils qw/ uniq /;

$annotFile = $ARGV[0];
$PSHOME = $ARGV[1];
$chooseManuallyFile = $ARGV[2];

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;

$scale = exp(23);
$freqFile = $ENV{'DATA_DIR'}."/HLA_FREQ.txt";

# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};


@lines = ();
open AFILE, $annotFile;
while($line = <AFILE>){
	chomp($line);
	if(!($line =~ /exon/ || $line =~ /splice/) || $line !~ /;QSI/){
		next;
	}
	@f = split(/\t/,$line);
	$sample = $f[1];
	$allele = $f[2];
	$pos = $f[3];
	$peptideChange = $f[16];
	$key = $sample.":".$allele;
        #print "$sample\t$key\n";
	$keyCountHash{$key} =  $keyCountHash{$key} + 1;
	push @lines, $line;
}

# remove all lines where keyCount > 1

@lines2 = ();
%sampleChangeCountHash = (); # count of all sample:peptideChange occurrenece for a particular sample and peptideChange
%sampleChangeLinesHash = (); # count of all lines corresponding to a particular sample:peptideChange
%sampleChangesHash = (); # list of all peptideChanges (including duplicates) for a sample
%sampleCountHash = (); # count of total peptideChanges (including duplicates) for a sample
@samples = ();
for($i = 0; $i <= $#lines; $i++){
	$line = $lines[$i];
	#print "LINE=$line\n";
	#next;
	@f = split(/\t/,$line);
        $sample = $f[1];
        $allele = $f[2];
        $pos = $f[3];
	$peptideChange = $f[16];
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
			print "$line\n";
		} else {
			@arr = @{$sampleChangeLinesHash{$sample}->{$change}};
			#print "#### choose 1 ####\n";
			for($k = 0; $k <= $#arr; $k++){
				$line = $arr[$k];
				@f = split(/\t/,$line);
				$tRefCount = $f[$trc];
			        $tAltCount = $f[$tac];
       				$nRefCount = $f[$nrc];
        			$nAltCount = $f[$nac];
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
                                $tRefCount = $f[$trc];
                                $tAltCount = $f[$tac];
                                $nRefCount = $f[$nrc];
                                $nAltCount = $f[$nac];
				if($tAltCount == $tAltCountMax){
					push @lines2, $line;
				}
			}
			if( ($#lines2+1) == 1){
				print "$chosenLine\n";
			} else{
				#print "can't choose\n";
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
