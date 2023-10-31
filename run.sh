#!/bin/bash

source /home/dbarila/miniconda3/bin/activate devo-grn  && python script.py 
wait
echo "script done" | mail -s "script done" danabr93@gmail.com


# run from qlogin: qsub -cwd -m ea -l m_mem_free=128G ./run.sh
