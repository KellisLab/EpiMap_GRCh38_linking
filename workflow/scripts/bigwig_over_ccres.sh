#!/usr/bin/env bash
BSSID="${1}"
MARK="${2}"
ANNOT_TABLE="${3}"
BED="${4}"
OUT_GZ="${5}"
DNAME=${6:-.}

echo "bw: $*"
if [ ! -d ${DNAME} ]; then
    mkdir ${DNAME}
fi

cd ${DNAME}

BIGWIG=$(< ${ANNOT_TABLE} awk -v bssid=${BSSID} -v mark=${MARK} '$1 == bssid && $2 == mark { print $3 ".bigWig" }' )
echo "BIGWIG: ${BIGWIG}"
echo "OUT_GZ: ${OUT_GZ}"
if [ ! -f ${OUT_GZ} ]; then
    echo "${BIGWIG}" | grep -qE "^impute"
    if [ $? -eq 0 ]; then
        echo "Using imputed"
        URL="https://epigenome.wustl.edu/epimap/data/imputed/${BIGWIG}"
    else
        echo "Using observed"
        URL="https://epigenome.wustl.edu/epimap/data/observed/${BIGWIG}"
    fi
    wget ${URL} -O ${BIGWIG}
    mkdir -p `dirname ${OUT_GZ}`
    bigWigAverageOverBed ${BIGWIG} ${BED} ${OUT_GZ%.gz} -minMax
    rm ${BIGWIG} -f
    gzip -9 ${OUT_GZ%.gz}
else
	echo "File ${OUT_GZ} already exists"
	exit 1
fi
exit 0
