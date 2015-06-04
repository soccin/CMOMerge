#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

python2.7 $SDIR/MergePortal $*
