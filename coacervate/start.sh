#!/bin/bash
while getopts m:d:o:a: flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        d) yagna_daemon=${OPTARG};;
        o) outputs=${OPTARG};;
        a) arb_sm=${OPTARG};;
    esac
done

if [[ $yagna_daemon = "on" ]]
    then
    # Start the yagna daemon and put it in the background
    # Using `-v yagna_datadir:/home/coacervate/.local/share/yagna`
    # will persist yagna datadir between containers
    echo "Starting yagna daemon in screen session.."
    screen -dmS yagna_daemon yagna service run
    sleep 15

    get_appkey () {
        echo $(yagna app-key list --json | jq -r .values[0][1])
    }

    if [ -z $(get_appkey) ]
    then
        echo "Problem initializing requestor, check permissions.."
        exit 1
    elif [ $(get_appkey) = "null" ]
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
