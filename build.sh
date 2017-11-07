echo "Running the build.sh Script"
cpanm -i Math::BaseCalc
mkdir -p $PREFIX/opt/
mkdir -p $PREFIX/shared/
cp $SRC_DIR/src/picard/* $PREFIX/opt/
cp -r $SRC_DIR/src/data/ $PREFIX/shared/polysolver_data
novoindex $PREFIX/shared/polysolver_data/abc_complete.nix $PREFIX/shared/polysolver_data/abc_complete.nix
