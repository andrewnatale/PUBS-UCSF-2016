#!/bin/bash

cat list.txt | while read a 
do
    python nd2.py $a
    done 
