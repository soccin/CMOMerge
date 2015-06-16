#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

/opt/common/CentOS_6/python/python-2.7.8/bin/python2.7 $SDIR/MergePortal $*
