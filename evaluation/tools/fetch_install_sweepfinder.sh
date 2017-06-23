TOOL=sweepfinder ## MODIFY THIS

rm -r $TOOL
mkdir $TOOL
cd $TOOL
wget http://people.binf.ku.dk/rasmus/webpage/SF/SF.tar.gz
tar -zxvf SF.tar.gz

cd SF

#Sequential
make
cp SweepFinder ../ ## MODIFY THIS
