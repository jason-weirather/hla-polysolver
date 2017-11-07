echo "Running the build.sh Script"
echo "Adding perl scripts"
cpanm -i Math::BaseCalc
mkdir -p $PREFIX/bin/
mkdir -p $PREFIX/jar/
mkdir -p $PREFIX/share/
mkdir -p $PREFIX/scripts/
echo "Adding picard jar files"
cp $SRC_DIR/src/picard/* $PREFIX/jar/
echo "Adding polysolver data"
cp -r $SRC_DIR/data/ $PREFIX/share/polysolver_data
echo "Rebuilding novoalign index"
novoindex $PREFIX/share/polysolver_data/abc_complete.nix $PREFIX/share/polysolver_data/abc_complete.fasta
echo "Put scripts into the polysolver build"
cp $SRC_DIR/scripts/* $PREFIX/scripts/
