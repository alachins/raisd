wget http://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar -xvzf gsl-latest.tar.gz
curpath=$(pwd)
mkdir gsl
cd gsl-*
./configure --prefix=$curpath/gsl && make && make install
cd ..
rm Makefile
cp makefiles/Makefile.GSL Makefile
make clean
make
cd source && ln -s ../gsl/include/gsl/ .
