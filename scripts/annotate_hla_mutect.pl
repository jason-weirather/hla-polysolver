#!/usr/bin/env perl

# usage: ./annotate_hla_mutect.pl HNSC-TCGA-BA-6873 HNSC-TCGA-BA-6873 /cga/wu/sachet/hla/hla_caller/capture/polysolver_based_muts_080814/data/a_complete.3100.new.eb.fasta /cga/wu/sachet/hla/hla_caller/capture/polysolver_based_muts_080814/data/b_complete.3100.new.eb.fasta /cga/wu/sachet/hla/hla_caller/capture/polysolver_based_muts_080814/data/c_complete.3100.new.eb.fasta $PSHOME


use POSIX;

$indiv = $ARGV[0];
$dir = $ARGV[1];
$aFile = $ARGV[2];
$bFile = $ARGV[3];
$cFile = $ARGV[4];
$PSHOME = $ARGV[5];

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;
$lib1 = $PSHOME."/scripts/common_functions_hla.pl";
require $lib1;

$outFile = $dir."/".$indiv.".mutect.unfiltered.annotated";

open INFILE, $inFile;
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";


($ref1,$ref2,$ref3) = @{get_seq_start_codon_splice_sites($aFile)};
%aH1 = %{$ref1};
%aH2 = %{$ref2};
%aH3 = %{$ref3};

($ref1,$ref2,$ref3) = @{get_seq_start_codon_splice_sites($bFile)};
%bH1 = %{$ref1};
%bH2 = %{$ref2};
%bH3 = %{$ref3};

($ref1,$ref2,$ref3) = @{get_seq_start_codon_splice_sites($cFile)};
%cH1 = %{$ref1};
%cH2 = %{$ref2};
%cH3 = %{$ref3};

%seqHash = (%aH1,%bH1,%cH1);
%startHash = (%aH2,%bH2,%cH2);
%splicePosHash = (%aH3,%bH3,%cH3);

# get the header file

@t = `ls $dir/call_stats*`;
$csFile = $t[0];
chomp($csFile);
open CSFILE, $csFile || die "Can not open $csFile\n";
while(($header = <CSFILE>) !~ /^contig/){}
chomp($header);
close CSFILE;

print OUTFILE "individual\t$header\texon_number\tintron_number\tprotein_change\n";

# get the list of call_stats calls

@t=`grep KEEP $dir/call_stats*`;

for($i1 = 0; $i1 <= $#t; $i1++){
	$line = $t[$i1];
	chomp($line);
	$line =~ s/.*://;
	@f = split(/\t/,$line);
	#$f[0] =~ s/.*://;
        $id = $f[0];
        $pos = $f[1];
        $context = $f[2];
        $ref = $f[3];
        $alt = $f[4];
        $seq = $seqHash{$id};
        @f2 = split(/\|/,$seq);
        $len = 0;
        $cdna = ();
        $status = ();
        $exonCount = 0;
        $intronCount = 0;
        $proteinChange = "-1";
        for($i = 0; $i <= $#f2; $i++){
                #print "i=$i\n";
                if(($i % 2) == 1){
                        $status = "exon";
                        $exonCount++;
                } else{
                        $status = "intron";
                        $intronCount++;
                }
                $len = $len + length($f2[$i]);
                        if($len >= $pos){
                                last;
                        }
                if(($i % 2) == 1){
                        #print "len: $len\t",length($f2[$i]),"\n$f2[$i]\n";
                        $cdna = $cdna.$f2[$i];
                }
        }
        #print "$indiv\t$id\tstatus=$status\tpos=$pos\texCount=$exonCount\t$intronCount\t";
        # check if it's a splice site mutation
        $isSplice = 0;
        @splicePos = @{$splicePosHash{$id}};
        #print "splicePos = @splicePos\n";
        if(found_arr($pos,\@splicePos)){
                $isSplice = 1;
        }

        if($status eq "exon"){
                $len = $len - length($f2[$i]);
                $diff = $pos - $len;
                $cdna2 = $cdna; #changed cdna
                $subseq = substr($f2[$i],0,$diff); # gets upto the changed base
                $subseq2 = substr($f2[$i],0,$diff - 1).$alt; # gets upto 1 base before the changed base
                $cdna = $cdna.$subseq;
                $cdna2 = $cdna2.$subseq2;
                $rem = length($cdna) % 3;
                if($rem == 1){
                        $aaNum = floor(length($cdna) / 3) + 1;
                        $cdna = $cdna.substr($f2[$i],$diff,2);
                        $cdna2 = $cdna2.substr($f2[$i],$diff,2);
                } elsif($rem == 2){
                        $aaNum = floor(length($cdna) / 3) + 1;
                        $cdna = $cdna.substr($f2[$i],$diff,1);
                        $cdna2 = $cdna2.substr($f2[$i],$diff,1);
                } elsif($rem == 0){
                        $aaNum = length($cdna) / 3;
                        $cdna = $cdna;
                        $cdna2 = $cdna2;
                }
                @arr1 = @{getPeptide($cdna,0)};
                $peptide1 = $arr1[0];
                $peptideLen1 = $arr1[1];
                $peptideStopPos1 = $arr1[2];
                $peptideTerminated1 = $arr1[3];
                @f3_1 = split(//,$peptide1);
                $aa1 = $f3_1[-1];
                @arr2 = @{getPeptide($cdna2,0)};
                $peptide2 = $arr2[0];
                $peptideLen2 = $arr2[1];
                $peptideStopPos2 = $arr2[2];
                $peptideTerminated2 = $arr2[3];
                @f3_2 = split(//,$peptide2);
                $aa2 = $f3_2[-1];
                $proteinChange = "p.".$aa1.$aaNum.$aa2;
                if($isSplice == 1){
                        $status = "splice_site";
                        $proteinChange = "p.".$aa1.$aaNum."_splice";

                }
                #print "i=$i\tdiff=$diff\tpos=$pos\n$cdna\n";
        } else {
                # check if it's a splice site mutation
                @splicePos = @{$splicePosHash{$id}};
                #print "splicePos = @splicePos\n";
                if(found_arr($pos,\@splicePos)){
                        $status = "splice_site";
                        #print "status = splice_site\n";
                        #print "cdna=",length($cdna),"\n";
                        $rem = length($cdna) % 3;
                        if($rem == 1){
                                $aaNum = floor(length($cdna) / 3) + 1;
                                $cdna = $cdna.substr($f2[$i+1],0,2);
                        }elsif($rem == 2){
                                $aaNum = floor(length($cdna) / 3) + 1;
                                $cdna = $cdna.substr($f2[$i+1],0,1);
                        } elsif($rem == 0) {
                                $aaNum = length($cdna) / 3;
                        }
                        @arr = @{getPeptide2($cdna,0)};
                        $peptide = $arr[0];
                        $peptideLen = $arr[1];
                        $peptideStopPos = $arr[2];
                        $peptideTerminated = $arr[3];
                        @f3 = split(//,$peptide);
                        $aa = $f3[-1];
                        $proteinChange = "p.$aa".$aaNum."_splice";
                        #print "proteinChange = $proteinChange\n";
                }
        }
        #print "$ref\t$alt\tproteinChange = $proteinChange\n";
        print OUTFILE "$indiv\t$line\t$status\t$exonCount\t$intronCount\t$proteinChange\n";

	#$line2 = join("\t",@f2);
	#print "$indiv\t$line2\n";
}


close INFILE;
close OUTFILE;



########

sub get_seq_start_codon_splice_sites { # start sub get_seq_start_codon_splice_sites

        # function takes a .eb.fasta file and returns a sequence hash, a start codon start pos hash,
        # and a splice site pos hash for all the alleles in the file
        # Note that the seq hash has |
        # the start codon start pos hash is base 1, in the original sequence with no |
        # the splice pos hash is base 1, in the original sequence with no |

        my ($file) = @_;

        my (%seqHash,%startHash,%splicePosHash,$i,@f,$line,$splicePos,$seq,$len,$subseq,$subseqLen,
                $start,@arr);

        open FILE, $file || die "get_seq_start_codon_splice_sites: Cannot open $file\n";
        %seqHash = ();
        %startHash = (); # 1-base, in original sequence without |
        %splicePosHash = ();
        while($line = <FILE>){
                chomp($line);
                $line = substr($line,1);
                $seq = <FILE>;
                chomp($seq);
                $seqHash{$line} = $seq;
                $start = index($seq,"|") + 1;
                @f = split(/\|/,$seq);
                @splicePos = ();
                $len = 0;
                for($i = 0; $i < $#f; $i++){
                        $subseq = $f[$i];
                        $subseqLen = length($subseq);
                        $len = $len + $subseqLen;
                        if($i==0){
                                next; # since the start codon is the start of the first exon
                        }
                        if(($i % 2)==1){
                                push @splicePos, $len;
                                push @splicePos, ($len + 1);
                                push @splicePos, ($len + 2);
                                push @splicePos, ($len + 3);
                                push @splicePos, ($len + 4);
                                push @splicePos, ($len + 5);
                                push @splicePos, ($len + 6);
                                push @splicePos, ($len + 7);
                                push @splicePos, ($len + 8);
                                #push @splicePos, ($len + 9);
                                #push @splicePos, ($len + 10);
                        } else{
                                #push @splicePos, ($len - 9);
                                #push @splicePos, ($len - 8);
                                #push @splicePos, ($len - 7);
                                #push @splicePos, ($len - 6);
                                #push @splicePos, ($len - 5);
                                #push @splicePos, ($len - 4);
                                #push @splicePos, ($len - 3);
                                #push @splicePos, ($len - 2);
                                push @splicePos, ($len - 1);
                                push @splicePos, $len;
                        }

                }
                $splicePosHash{$line} = [ @splicePos ];
                #print $line,"\t",$start,"\n";
        }

        close FILE;

        @arr = (\%seqHash,\%startHash,\%splicePosHash);

        return \@arr;

} # end sub get_seq_start_codon_splice_sites
                                                    

########

sub print_splice_pos {

        my ($ref) = @_;

        my (@arr,$i,$depth,$localCount,$lineState);

        $depth = 10;

        @arr = @{$ref};

        $localCount = 0;
        $lineState = 0;

        for($i = 0; $i <= $#arr; $i++){

                $localCount++;

                if($localCount < 10){

                print "$arr[$i] ";

                } else {

                        if($lineState %2  == 0){

                                print "$arr[$i] - ";

                                $localCount = 0;

                        } else{

                                print "$arr[$i]\n";

                                $localCount = 0;
                        }

                        $lineState++;

                }
        }
}

########

sub getPeptide2 { # start sub getPeptide

        # this function taken in a DNA sequence and cds start position and translates it
        # the cds start position has to be 0-base
        # it returns
                # peptide: the translated peptide
                # length: length in base pairs of the translated peptide
                # stopPos: the stop position of the codon in the original submitted sequence in bp
                        # 0-base
                        # inclusive of the stop codon
                        # is -1 is peptideTerminated = 0
                # peptideTerminated: 1 if an in-frame stop codon was found, 0 otherwise


        # it then takes the changePoint and pepWindow arguments to get a pepWindow # peptides
        # around the changePoint position
        # if pepWindow == 0, it returns a 9-mer and 10-mer amino acid chain around the
        # changePoint
        # for a 9-mer it takes 4 aa on either side of the changePoint + the changePoint
        # for a 10-mer it takes 5 aa upstream, changePoint and 4 aa downstream and also
        # 4 + 1 + 5
        # if a premature stop codon prevents the getting of the the peptide around the
        # changePoint, an error flag is also returned

        my ($seq,$cdsStart) = @_;

        my (%codonMap,$i,$seq2,$length,$found,$codon,$aa,$peptide,$seq3,@arr,$stopPos);

        my ($peptideTerminated);

        if($cdsStart=~/^$/ || $cdsStart < 0){

                print "ERROR: getPeptide: cdsStart = $cdsStart\n";

                $peptide = "NA";

                $length = 0;

                $stopPos = -1;

                $peptideTerminated = 0;

                @arr = ($peptide,$length,$stopPos,$peptideTerminated);

                return \@arr;
        }

        %codonMap = %{codonMap(1)};

        # get new orf

        $seq2 = substr($seq,$cdsStart);

        #print "getPeptide:cdsStart=$cdsStart\nseq=$seq\nseq2=$seq2\n";

        $seq3 = ();

        $peptide = ();

        $length = 0;

        $found = 0;

        $peptideTerminated = 0;

        #print "seq=$seq\nseq2=$seq2\n";

        while(!$found){

                $codon = substr($seq2,0,3);

                $codon =~ tr/a-z/A-Z/;

                if(length($codon) < 3){

                        last;

                }
                $aa = $codonMap{$codon};

                #print "codon=$codon\naa=$aa\n";

                if($aa=~/^$/){

                        die "ERROR: getPeptide:for cdsStart = $cdsStart codon=$codon aa=$aa in
                                seq2=$seq2 length=$length\n";

                }

                $length = $length + 3;

                $peptide = $peptide.$aa;

                #print "codon = $codon\n";

                if($codon eq "TGA" || $codon eq "TAG" || $codon eq "TAA"){

                        $found = 1;

                        last;

                }


                $seq2 = substr($seq2,3);

                if($seq2=~/^$/){

                        last;

                }

        }

        if($found==1){

                $peptideTerminated = 1;
        }

        # peptide now has the full amino acid sequence

        #print "getPeptide: peptide=$peptide\n";

        if($peptideTerminated){

                        $stopPos = $cdsStart + $length;

        } else{

                        $stopPos = -1;
        }

        @arr = ($peptide,$length,$stopPos,$peptideTerminated);

        return \@arr;

        #die;

        # find the hla9 and hla10 peptides


} # end sub getPeptide


