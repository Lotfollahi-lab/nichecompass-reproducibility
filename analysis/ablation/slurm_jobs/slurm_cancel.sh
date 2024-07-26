#!/bin/bash

for j in `seq 13093202 13093329` ; do
	scancel $j
	echo $j
done
