#!/bin/bash

cd ../
DIR=`pwd`

DIR=`basename $DIR`

COPYDIR=$DIR"_COPY"
mkdir ../"$COPYDIR"

find . ! -name \*.bin ! -name \*.out ! -name \*.err | cpio -p ../"$COPYDIR"

exit 0