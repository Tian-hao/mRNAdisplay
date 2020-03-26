#!/bin/bash
#call mutations
echo $SGE_TASK_ID
for i in `seq 0 9`; do
   my_task_id=$((SGE_TASK_ID + i))
   echo $my_task_id
   python /u/home/t/tianhao/display/script/Mapper1_stringent.py `sed -n ${my_task_id}p /u/home/t/tianhao/display/ref/filename.txt`
done


