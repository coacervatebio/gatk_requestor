#!/bin/bash
while getopts m: flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
    esac
done

# Start the yagna daemon and put it in the background
# Using `-v /host/yagna/datadir:/home/requestor/.local/share/yagna`
# will persist yagna datadir between containers
echo "Starting yagna daemon in screen session.."
screen -d -m -S yagna_daemon yagna service run
sleep 15

get_appkey () {
    echo $(yagna app-key list --json | jq -r .values[0][1])
}

if [ $(get_appkey) = "null" ]
then
    echo "No appkey found, creating requestor.."
    yagna app-key create requestor
    sleep 5
    echo "Funding.."
    yagna payment fund
    sleep 10
else
    echo "Found existing appkey"
fi

echo "Exporting app key and initializing yagna as sender.."
export YAGNA_APPKEY=$(get_appkey)
sleep 3
yagna payment init --sender # init sender

echo "Running script in ${mode} mode.."
case "${mode}" in
    interactive) /bin/bash;;
    local) snakemake -c2 -s=/mnt/workflow/rules/local.smk;;
    golem) snakemake -c2 -s=/mnt/workflow/rules/hc_golem.smk;;
    req_only) python /mnt/workflow/scripts/requestor.py;;
    *) echo "Please select interactive, local, or golem mode";;
esac
