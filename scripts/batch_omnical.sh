#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -e grid_output
#$ -o grid_output
#$ -l paper
#$ -l h_vmem=5G
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
        echo python scripts/omnical_PSA128.py $dir/zen.*.${CHUNK}*.${EXT} -C ${CALFILE} -p xx,yy -t aug2014  -d ${TAG1}_245${TAG2}.${CHUNK} -i doc/redundantinfo_PSA64_7ba_7bu_08-15-2014.bin --datapath /data2/home/hz2ug/omnical_old/results  -o /data4/paper/2012EoR/psa_live/PSA64_omnical_results_Aug_2014  --health 3,4 -k
        echo ----------------------------------------------------------------------------------------------------------------------------
        python scripts/omnical_PSA128.py $dir/zen.*.${CHUNK}*.${EXT} -C ${CALFILE} -p xx,yy -t aug2014  -d ${TAG1}_245${TAG2}.${CHUNK} -i doc/redundantinfo_PSA64_7ba_7bu_08-15-2014.bin --datapath /data2/home/hz2ug/omnical_old/results  -o /data4/paper/2012EoR/psa_live/PSA64_omnical_results_Aug_2014  --health 3,4 -k
    done;
done;

