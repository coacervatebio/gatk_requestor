#!/bin/bash
while getopts m:y:o:a: flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        y) yagna=${OPTARG};;
        o) outputs=${OPTARG};;
        a) arb_sm=${OPTARG};;
    esac
done

if [[ $yagna = "on" ]]
    then
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
fi

echo "Running script in ${mode} mode.."
case "${mode}" in
    interactive) /bin/bash;;
    local) snakemake -c2 -s=/data/workflow/rules/local.smk;;
    golem) snakemake -c2 -s=/data/workflow/Snakefile;;
    specific) snakemake ${outputs} -c2 -s=/data/workflow/Snakefile ${arb_sm};;
    req_only) python /data/workflow/scripts/requestor.py;;
    *) echo "Please select interactive, local, or golem mode";;
esac
