#!/bin/bash
trap "exit 255" SIGINT SIGTERM
PROG=${HOME}/bin/sunspot

function testfunction {
    ${PROG} -l $1 --analyze $2
}

if [[ $# != 2 ]]; then
    echo "Usage: $0 <path> <tool>"
    echo "    path: location at which to search for checkpoints."
    echo "    tool: analysis tool to execute on checkpoints."
    exit -1
fi

for i in `find $1 -name "checkpoint-*.xml.gz"`; do
    RESULTS_DIR=`dirname $i`
    pushd ${RESULTS_DIR}
    testfunction `basename $i` $2
    popd > /dev/null
done
