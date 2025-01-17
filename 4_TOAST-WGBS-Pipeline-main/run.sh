#!/bin/bash

publishdir=./results/
echo "INFO: Running nextflow. workdir is $PWD. publishdir is $publishdir"

export NXF_OPTS="-Xmx8G"
export _JAVA_OPTIONS="-Xmx8G"
ulimit -u 2048

nohup nextflow run main.nf -c nextflow.config -profile nscc -resume &
