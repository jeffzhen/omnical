#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o /data4/paper/hz2ug/2015PSA64/grid_output/
#$ -N JD6_v9
#$ -j y
#$ -l paper
#$ -l h_vmem=2.1G


DIRS=`/data4/paper/2013EoR/2456620_SAK/pull_args.py $*`
EXT=uvcRREcAC

CALFILE=psa6240_v003
INFO=/data2/home/hz2ug/omnical/doc/redundantinfo_first_cal_6285.57756.bin
P=/data2/home/hz2ug/omnical/doc/calpar_first_cal_6285.57756.p
BINDATADIR=/data4/paper/hz2ug/2015PSA64/binary_data/
OPDIR=/data4/paper/hz2ug/2015PSA64/

MODEL_TREASURE=/data4/paper/hz2ug/2015PSA64/2015PSA64_march2015v5_0.v0foreground.treasure
MODEL_NOISE=/data4/paper/hz2ug/2015PSA64/PSA64_chisq_model_feb2015_v2_xx.omnichisq,/data4/paper/hz2ug/2015PSA64/PSA64_chisq_model_feb2015_v2_yy.omnichisq

N_TREASURE=2

TAG=march2015v9

source /usr/global/paper/CanopyVirtualEnvs/JZ_Omnical/bin/activate

for dir in ${DIRS}; do
    TAG1=`echo ${dir} | cut -d /  -f 6`
    TAG2=`echo ${dir} | cut -d /  -f 6 | cut -d a -f 2`
    TREASURE=/data4/paper/hz2ug/2015PSA64/2015PSA64_${TAG}_`expr $TAG2 % $N_TREASURE`.treasure
        
        
        echo "==============================================================================================================="
        echo omnical_PSA128.py ${dir} -C ${CALFILE} -p xx,yy -t ${TAG}  -d `echo ${dir} | cut -d /  -f 7` -i ${INFO} -r ${P} -o ${OPDIR}  --health 3 --flag -f --treasure ${TREASURE} --flagsigma 5 --flagf 3 -s --mem 1e8 --model_noise ${MODEL_NOISE} --chemo 2 --chemot 2 --flagt 50 --skip_sun --chisq_leniency 1.5 --model_treasure ${MODEL_TREASURE}
        echo ----------------------------------------------------------------------------------------------------------------------------
        omnical_PSA128.py ${dir} -C ${CALFILE} -p xx,yy -t ${TAG}  -d `echo ${dir} | cut -d /  -f 7` -i ${INFO} -r ${P} -o ${OPDIR}  --health 3 --flag -f --treasure ${TREASURE} --flagsigma 5 --flagf 3 -s --mem 1e8 --model_noise ${MODEL_NOISE} --chemo 2 --chemot 2 --flagt 50 --skip_sun --chisq_leniency 1.5 --model_treasure ${MODEL_TREASURE}

        echo "==============================================================================================================="
done;



