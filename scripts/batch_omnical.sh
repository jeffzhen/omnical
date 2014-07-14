#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -e grid_output
#$ -o grid_output
#$ -l paper
#$ -l h_vmem=7G
DIRS=`pull_args.py $*`
CALFILE=psa6240_v003
EXT=uvcRREcAC

for dir in ${DIRS}; do
    TAG1=`echo ${dir} | cut -d /  -f 6`
    TAG2=`echo ${dir} | cut -d /  -f 6 | cut -d a -f 2`
    echo ==============================================================================================================================
    echo TAG2 ${TAG2}
    echo ==============================================================================================================================
    for CHUNK in `seq 1 1 6`; do
        echo ./omnical_PSA64.py -C ${CALFILE} -p xx,yy --tag=${TAG1}_245${TAG2}.${CHUNK} $dir/zen.*.${CHUNK}*.${EXT} --skip --add --nadd 7
        echo ----------------------------------------------------------------------------------------------------------------------------
        ./omnical_PSA64.py -C ${CALFILE} -p xx,yy --tag=${TAG1}_245${TAG2}.${CHUNK} $dir/zen.*.${CHUNK}*.${EXT} --skip --add --nadd 7
    done;
done;

