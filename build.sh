echo "Running the build.sh Script"
cpanm -i Math::BaseCalc
mkdir -p $PREFIX/opt/
cp $SRC_DIR/src/picard/* $PREFIX/opt/
