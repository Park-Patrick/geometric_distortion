#!/bin/bash

if [ "$#" -lt 1 ];then
	echo "Usage: $0 <install folder (absolute path)> <checkout hash (optional)>"
	echo "For sudoer recommend: $0 /opt"
	echo "For normal user recommend: $0 $HOME/app"
	exit 0
fi

echo -n "installing NiftyReg git version..." #-n without newline

DEST=$1
mkdir -p $DEST

if [ "$#" -lt 2 ]
then
#June 14, 2016 release (perhaps most compatible with Lau NeuroImage Geo-distort paper)
#checkouthash=e175a1

#latest release as of Dec 4:   Nov 16, 2017
checkouthash=6bf84b

else
checkouthash=$2
fi
echo Using checkout hash $checkouthash

D_DIR=$DEST/niftyreg-git
if [ -d $D_DIR ]; then
	rm -rf $D_DIR
fi

NIFTY_SRC=$D_DIR/src
NIFTY_DIR=$D_DIR

mkdir -p $NIFTY_DIR
git clone https://git.code.sf.net/p/niftyreg/git $NIFTY_SRC
pushd $NIFTY_SRC
checkout $checkouthash
popd

pushd $NIFTY_DIR
cmake $NIFTY_SRC \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX=$NIFTY_DIR  && \
  make -j$(nproc) && \
  make install 
popd

#rename binaries to add _git
for file in `ls $NIFTY_DIR/bin/reg*`
do
  mv $file ${file}_git
done



PROFILE=~/.bashrc

if grep -xq "PATH=$NIFTY_DIR/bin:\$PATH" $PROFILE #return 0 if exist
then 
 	echo "PATH=$NIFTY_DIR/bin" in the PATH already.
else
 	echo "export PATH=$NIFTY_DIR/bin:\$PATH" >> $PROFILE
	echo "export LD_LIBRARY_PATH=$NIFTY_DIR/lib:\$LD_LIBRARY_PATH" >> $PROFILE
fi

source $PROFILE


