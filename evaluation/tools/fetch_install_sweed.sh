TOOL=sweed ## MODIFY THIS

rm -r $TOOL
mkdir $TOOL
cd $TOOL
wget https://github.com/alachins/$TOOL/archive/master.zip
unzip master.zip
cd $TOOL-master

#Sequential
make -f Makefile.gcc ## MODIFY THIS
cp SweeD ../ ## MODIFY THIS

make clean -f Makefile.gcc
make -f Makefile.PTHREADS.gcc
cp SweeD-P ../
