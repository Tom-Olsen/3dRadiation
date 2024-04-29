#!/bin/bash
# Get the current user's username
USER_NAME=$(whoami)

# Get a list of all jobs belonging to the current user
JOB_IDS=$(squeue -u $USER_NAME -h -o %A)

# Cancel each job
for JOB_ID in $JOB_IDS
do
    scancel $JOB_ID
done