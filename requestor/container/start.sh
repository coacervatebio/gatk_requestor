#!/bin/bash
  
# turn on bash's job control
set -m
  
# Maybe screen -d -m
# Start the yagna daemon and put it in the background
yagna service run &

# Wait for daemon init and then init requestor
sleep 5
yagna payment init --sender

/bin/bash