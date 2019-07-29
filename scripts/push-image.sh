#!/bin/sh
# Pushes docker image for service to docker repo. Takes 1 arg:
# arg1: docker image tag
# Example usage: scripts/push-image.sh v1.0.0
set -e

# cd to project root directory
cd "$(dirname "$(dirname "$0")")"

# import common stuff
source scripts/lib/common.sh

# Exit if two args not given
if [[ -z "$1" ]]; then
    echo "One argument required.";
    exit;
fi

docker push ${DOCKER_ORG}/${IMAGE_NAME}:${1}

docker tag ${DOCKER_ORG}/${IMAGE_NAME}:${1} ${DOCKER_ORG}/${IMAGE_NAME}:latest
docker push ${DOCKER_ORG}/${IMAGE_NAME}:latest