#! /bin/sh
##########################################
# get RF list file from data as 
# input file for stack/pierc_new_n.inp 
##########################################

rm -rf datalist.txt; touch datalist.txt

ls -d ???? > sta.tmp
#ls -d ??.???? >> sta.tmp

number=0
while read sta
do
echo ${sta}
num=`ls ${sta}/*nsl.sac.cut100rv.mout | awk -F/ '{print $2}' | wc -l`
echo "${sta}" >> datalist.txt
echo "${num}" >> datalist.txt
ls ${sta}/*nsl.sac.cut100rv.mout | awk -F/ '{print $2}' >> datalist.txt
echo " " >> datalist.txt
number=`echo "${number}+${num}" | bc`
done < sta.tmp
rm -rf sta.tmp
echo " " >> datalist.txt
echo " " >> datalist.txt
echo " " >> datalist.txt
echo " " >> datalist.txt
echo "Total RF number: ${number}" >> datalist.txt
echo " " >> datalist.txt
