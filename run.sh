#!/bin/bash

source /home/ybahar/miniconda3/bin/activate devo-grn  && python script.py 
wait
echo "script done" | mail -s "script done" danabr93@gmail.com
