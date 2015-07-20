#!/bin/bash

curdate="$(date +'%m-%d-%Y')"
nohup doit --verbosity 2 -n 10 > log_doit_$curdate.txt &

