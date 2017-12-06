#!/bin/bash

if [ "$#" -lt 1 ];then
	echo "Usage: $0 <install folder (absolute path)>"
	echo "For sudoer recommend: $0 /usr/local or $0 /opt (recommend)"
	echo "For normal user recommend: $0 $HOME/app"
	exit 0
fi

echo -n "installing freesurfer minimal (mri_convert only)..." #-n without newline

DEST=$1
mkdir -p $DEST

S_DIR=$DEST
if [ -d $S_DIR/freesurfer_minimal ]; then
	rm -rf $S_DIR/freesurfer_minimal
else
    mkdir -p $S_DIR/freesurfer_minimal
fi

mkdir -p $S_DIR/freesurfer_minimal/bin

#get mri_convert
curl -L --retry 5 https://www.dropbox.com/s/rkd6q5ak4f0ugok/mri_convert.bin  -o $S_DIR/freesurfer_minimal/bin/mri_convert

#get license
curl -L --retry 5 https://www.dropbox.com/s/38gghuq2w7hl7pe/freesurfer_license -o $S_DIR/freesurfer_minimal/.license

chmod a+x $S_DIR/freesurfer_minimal/bin/mri_convert
