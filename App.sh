#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"
PYTHON=/opt/common/CentOS_6/python/python-2.7.8/bin/python2.7
REPO=/ifs/res/socci/Work/CMO/CMOMerge/bic-mskcc

usage() {
    echo "CMOMerge/App.sh -p projectTag -b baseProject -r repository -s python_path"
    echo "    projectTag = string to grep for to find projects to merge"
    echo "    baseProject = full project ID for the base project of merge"
    echo
    $PYTHON $SDIR/MergePortal -h 2>1 | egrep -v Version

}

if [ "$#" == "0" ]; then
    usage
    exit
fi

LABARG=""
TUMORARG=""
MERGEARG=""
FORCEARG=""
rm -f merge_exclude
while getopts ":b:p:l:t:n:r:s:e:hf" opt; do
    case $opt in
        b)
            BASETAG=$OPTARG
            ;;
        p)
            PROJECTTAG=$OPTARG
            ;;
        l)
            LABARG="-l "$OPTARG
            ;;
        t)
            TUMORARG="-t "$OPTARG
            ;;
        n)
            MERGEARG="-n "${OPTARG/-/_}
            ;;
        r)
            REPO=$OPTARG
            ;;
        s)
            PYTHON=$OPTARG
            ;;
        e)
            echo "/"$OPTARG"$" >>merge_exclude
            ;;
        f)
            FORCEARG="--force"
            ;;
        h)
            usage
            exit
            ;;

    esac
done
shift $((OPTIND-1))

if [ "$PROJECTTAG" != "" ]; then
    echo PROJECTTAG=$PROJECTTAG
    PROJECTS=$(find -L $REPO -name '*meta_study.txt' \
        | fgrep $PROJECTTAG \
        | fgrep -v cbe \
        | sed 's/.meta.*//' \
        | sort -r)
fi

if [ -e merge_exclude ]; then
    echo "Need to exclude"
    cat merge_exclude
    PROJECTS=$(echo $PROJECTS | tr ' ' '\n' | egrep -vf merge_exclude)
fi

BASEPROJECT=""
OTHERPROJECTS=""
if [ "$BASETAG" != "" ]; then
    i=1
    for pi in $PROJECTS; do
        if [[ $pi =~ "${BASETAG}" ]]; then
            BASEPROJECT=$pi
        else
            OTHERPROJECTS="--project "$pi" "$OTHERPROJECTS
        fi
    done

    echo BASEPROJECT=$BASEPROJECT

else
    i=1
    echo "Need to set baseProject"
    echo
    for pi in $PROJECTS; do
        echo $pi | perl -pe  "s|${REPO}.||" | awk -v i=$i '{print i":\t"$1}'
        i=$((i+1))
    done
    exit
fi

echo "-----"
echo "Pulling and updating repository"
pushd .
cd $REPO
hg pull
hg update
popd
echo $PWD
echo "done"
echo "--"
echo

ARGS="$LABARG $TUMORARG $MERGEARG $FORCEARG"
echo ARGS=$ARGS

$PYTHON $SDIR/MergePortal --project $BASEPROJECT \
    --root $REPO \
    $OTHERPROJECTS \
    $ARGS
