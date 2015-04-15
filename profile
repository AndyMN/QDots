#!/bin/bash

python -m cProfile -o $2 -s time $1
python profileRead.py $2
