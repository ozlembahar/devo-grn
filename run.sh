#!/bin/bash

source /home/ybahar/miniconda3/bin/activate devo-grn  && ./script.py 
wait
echo "script done" | mail -s "script done" danabr93@gmail.com
