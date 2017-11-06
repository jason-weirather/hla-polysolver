echo "Running the build.sh Script"
export PERL_MM_USE_DEFAULT=1
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
perl -MCPAN -e 'install Math::BaseCalc'
