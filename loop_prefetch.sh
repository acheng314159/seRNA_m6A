#!bin/bash
echo "usage: bash script.sh <SRP_ID> <THREAD_NUM> [loop_time]"
echo "P.S. The name of SRR_ID accession list must be SRP_ID.txt"
# "bash loop_prefetch.sh SRP039396 10 20" which means SRP039396 in the dir.
ID=$1
THREAD_NUM=$2
loop_time=$3
i=0
mkdir  $ID
while [ "$i" != "${loop_time:=100}" ]
do
 echo "loop_time: "$(($i+1))
 ls $ID| xargs -iFILE basename FILE .sra |sort| uniq | cat - $ID*.txt |sort |uniq -u |awk 'NF{print}' |xargs -iSRR_ID -P${THREAD_NUM} /Sum/lirui/soft/sratoolkit.2.9.6-1-ubuntu64/bin/prefetch  SRR_ID -O ./$ID --max-size 30G
 i=$(($i+1))
done
