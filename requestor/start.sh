#!/bin/bash
while getopts m: flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
    esac
done

# Start the yagna daemon and put it in the background
echo "Starting yagna daemon in screen session.."
screen -d -m -S yagna_daemon yagna --datadir /yagna service run
sleep 2

get_appkey () {
    echo $(yagna --datadir /yagna app-key list --json | jq -r .values[0][1])
}

if [ $(get_appkey) = "null" ]
then
    echo "No appkey found with --datadir /yagna, creating requestor.."
    yagna --datadir /yagna app-key create requestor
    sleep 5
    echo "Funding.."
    yagna --datadir /yagna payment fund
    sleep 10
else
    echo "Found existing appkey with --datadir /yagna"
fi

echo "Exporting app key and initializing yagna as sender.."
export YAGNA_APPKEY=$(get_appkey)
sleep 3
yagna --datadir /yagna payment init --sender # init sender

echo "Running script in ${mode} mode.."
case "${mode}" in
    interactive) /bin/bash;;
    automatic) snakemake -c1 -s=/mnt/workflow/rules/Snakefile;;
    *) echo "Please select automatic or interactive mode";;
esac
