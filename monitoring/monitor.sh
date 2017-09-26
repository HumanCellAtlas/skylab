#!/bin/bash
echo --- General Information ---
echo \#CPU: $(nproc)
echo Total Memory \(MB\): $(free -m | grep Mem | awk '{ print $2 }')  # prints in mb
echo Total Disk Space \(KB\): $(df -k | grep cromwell_root | awk '{ print $2}')  # prints in kb
echo
echo --- Runtime Information ---

function runtimeInfo() {
        # echo [$(date)]  # we don't really care about the date that much.

        # print memory usage as a percentage with up to 4 significant digits
        echo \* Memory usage \(%\): $(free -m | grep Mem | awk '{ OFMT="%.2f"; print ($3/$2)*100; }')%
        # print memory usage in megabytes
        echo \* Memory usage \(MB\): $(free -m | grep Mem | awk '{ print $3 }')
        # print disk usage as a percentage with up to 4 significant digits
        echo \* Disk usage \(%\): $(df -k | grep cromwell_root | awk '{ OFMT="%.2f"; print ($3/$2)*100; }')%
        # print disk usage in 1024 k blocks
        echo \* Disk usage \(KB\): $(df -k | grep cromwell_root | awk '{ print $3; }')

}

while true; do runtimeInfo; sleep 5; done