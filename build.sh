echo "Running the build.sh Script"
echo "Adding perl scripts"
cpanm -i Math::BaseCalc
mkdir -p $PREFIX/bin/
mkdir -p $PREFIX/jar/
mkdir -p $PREFIX/share/
mkdir -p $PREFIX/share/strelka
mkdir -p $PREFIX/scripts/
echo "Adding picard jar files"
cp $SRC_DIR/src/picard/* $PREFIX/jar/
echo "Adding polysolver data"
cp -r $SRC_DIR/data/ $PREFIX/share/polysolver_data
echo "Rebuilding novoalign index"
novoindex $PREFIX/share/polysolver_data/abc_complete.nix $PREFIX/share/polysolver_data/abc_complete.fasta
echo "Put scripts into the polysolver build"
cp $SRC_DIR/scripts/* $PREFIX/scripts/
echo "Build strelka"
echo "start building vcftools"
cd $SRC_DIR/include/ && \
  tar -xzf strelka-upstream-v1.0.11.tar.gz && \
  cd $SRC_DIR/include/strelka-upstream-v1.0.11/redist/vcftools-r837 && \
  make
echo "build main part"
cd $SRC_DIR/include/strelka-upstream-v1.0.11 && \
  ./configure --prefix=$PREFIX/share/strelka && \
  make
#echo "Trying to finish the vcftools install"
#cd $PREFIX/share/strelka/opt/vcftools && make
#mydisp=$(ls -l $PREFIX/share/strelka/opt/vcftools)
#echo "${mydisp}"
echo "Install hla-polysolver"
cp $SRC_DIR/scripts/shell_annotate_hla_mutations $PREFIX/bin/
cp $SRC_DIR/scripts/shell_call_hla_mutations_from_type $PREFIX/bin/
cp $SRC_DIR/scripts/shell_call_hla_type $PREFIX/bin/
