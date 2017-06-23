TOOL=sweepfinder2 ## MODIFY THIS

rm -r $TOOL
mkdir $TOOL
cd $TOOL
wget http://www.personal.psu.edu/mxd60/SF2.tar.gz
tar -zxvf SF2.tar.gz

cd SF2

#Sequential
make
cp SweepFinder2 ../ ## MODIFY THIS
