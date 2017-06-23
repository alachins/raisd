TOOL=raisd ## MODIFY THIS

rm -r $TOOL
mkdir $TOOL
cd $TOOL
wget https://github.com/alachins/$TOOL/archive/master.zip
unzip master.zip
cd $TOOL-master

#Sequential
make ## MODIFY THIS
cp bin/release/RAiSD ../ ## MODIFY THIS
