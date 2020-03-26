#!/bin/bash
#call mutations
echo $SGE_TASK_ID
python /u/home/t/tianhao/display/script/Mapper1.py `sed -n ${SGE_TASK_ID}p /u/home/t/tianhao/display/ref/filename.txt`


