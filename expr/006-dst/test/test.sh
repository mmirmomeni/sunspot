DIR=${HOME}/research/src/sunspot/expr/006-dst

function testfunction {
    ${HOME}/bin/sunspot -l $1 \
        --analyze sunspot_test_rmse
}

for i in `find ${DIR} -name "checkpoint-*.xml.gz"`; do
    IDIR=`dirname $i`
    TDIR=`basename ${IDIR}`
    mkdir -p ${TDIR}
    pushd ${TDIR}
    testfunction $i 
    popd
done
