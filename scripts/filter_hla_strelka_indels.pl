#!/usr/bin/env perl

#usage: ./filter_hla_strelka_indels.pl LUAD-TCGA-95-7562/LUAD-TCGA-95-7562.strelka_indels.unfiltered.annotated LUAD-TCGA-95-7562 LUAD-TCGA-95-7562 /cga/wu/sachet/hla/hla_caller/capture/polysolver_070914


use List::MoreUtils qw/ uniq /;

$annotFile = $ARGV[0];
$indiv = $ARGV[1];
$dir = $ARGV[2];
$PSHOME = $ARGV[3];

$outFile = $dir."/".$indiv.".strelka_indels.filtered.annotated";
$chooseManuallyFile = $dir."/".$indiv.".strelka_indels.ambiguous.annotated";

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;

$scale = exp(23);
$freqFile = $PSHOME."/"."data/HLA_FREQ.txt";


# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};


@lines = ();
open AFILE, $annotFile;
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";


# get the header

$header = <AFILE>;
chomp($header);

print OUTFILE "$header\n";

while($line = <AFILE>){
        chomp($line);
        #if(!($line =~ /exon/ || $line =~ /splice/) || $line !~ /;QSI/){
        #        next;
        #}
        @f = split(/\t/,$line);
        $sample = $f[0];
        $allele = $f[1];
        $pos = $f[2];
        $peptideChange = $f[-1];
	if($peptideChange == -1){
		next;
	}
        $key = $sample.":".$allele;
        #print "$sample\t$key\n";
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
		print "CHANGE=$change\t$sampleChangeCountHash{$sample}->{$change}\n";
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
				if($k == 0){
					print OUTFILE "$line\n";
				}
				print CMFILE "$line\n";
			}	
			#print "########\n";
		}
	}
	
}
close CMFILE;

close AFILE;
close OUTFILE;
