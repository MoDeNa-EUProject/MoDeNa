#!/bin/sh
mkdir temp
cd temp
wget http://math.lbl.gov/voro++/download/dir/voro++-0.4.6.tar.gz
tar -xf voro++-0.4.6.tar.gz
cd voro++-0.4.6
make
sudo make install
cd ../..
rm -rf temp
