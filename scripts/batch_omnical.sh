#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o /data4/paper/hz2ug/2015PSA64/grid_output/
#$ -N JD624_625
#$ -j y
#$ -l paper
#$ -l h_vmem=4G


DIRS=`/data4/paper/2013EoR/2456620_SAK/pull_args.py $*`
EXT=uvcRREcAC

CALFILE=psa6240_v003
INFO=/data2/home/hz2ug/omnical/doc/redundantinfo_first_cal_6285.57756.bin
P=/data2/home/hz2ug/omnical/doc/calpar_first_cal_6285.57756.p
BINDATADIR=/data4/paper/hz2ug/2015PSA64/binary_data/
OPDIR=/data4/paper/hz2ug/2015PSA64/
TREASURE=/data4/paper/hz2ug/2015PSA64/2015PSA64_624_625.treasure

source /usr/global/paper/CanopyVirtualEnvs/JZ_Omnical/bin/activate
for dir in ${DIRS}; do
    TAG1=`echo ${dir} | cut -d /  -f 6`
    TAG2=`echo ${dir} | cut -d /  -f 6 | cut -d a -f 2`
        echo omnical_PSA128.py ${dir} -C ${CALFILE} -p xx,yy -t feb2015  -d `echo ${dir} | cut -d /  -f 7` -i ${INFO} -r ${P} -o ${OPDIR}  --health 3 --flag -f --treasure ${TREASURE} --flagsigma 0.5 -s
        echo ----------------------------------------------------------------------------------------------------------------------------
        omnical_PSA128.py ${dir} -C ${CALFILE} -p xx,yy -t feb2015  -d `echo ${dir} | cut -d /  -f 7` -i ${INFO} -r ${P} -o ${OPDIR}  --health 3 --flag -f --treasure ${TREASURE} --flagsigma 0.5 -s

done;

