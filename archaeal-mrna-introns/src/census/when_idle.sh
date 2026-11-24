#!/bin/bash
# run "$@" once the 1-minute load is under 8 (checks every 60 s)
while :; do
  l=$(sysctl -n vm.loadavg | awk '{print $2}')
  awk -v l="$l" 'BEGIN{exit !(l<8)}' && break
  sleep 60
done
echo "$(date '+%H:%M:%S') start (load $l)"
exec "$@"
