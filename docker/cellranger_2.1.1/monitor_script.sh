#!/bin/bash
echo ==================================
echo =========== MONITORING ===========
echo ==================================
echo --- General Information ---
echo \#CPU: $(nproc)
echo Total Memory: $(free -h | grep Mem | awk '{ print $2 }')
echo Total Disk space: $(df -h | grep cromwell_root | awk '{ print $2}')
echo 
echo --- Runtime Information ---

function runtimeInfo() {
        echo [$(date)]
        echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
}

while true; do runtimeInfo; sleep 200; done