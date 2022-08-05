#!/bin/bash
# maybe take some optional args to run snakefile
  
# Start the yagna daemon and put it in the background
screen -d -m yagna service run
sleep 5

# Wait for daemon init and then init requestor
# TODO: logic to check if it has already been initialized
# if not fund it with testnet money, see one fo the suggested projects in disco
yagna payment init --sender # datadir?

