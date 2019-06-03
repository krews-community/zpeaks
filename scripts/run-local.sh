#!/bin/bash
set -e

# cd to project root directory
cd "$(dirname "$(dirname "$0")")"

java -jar build/*.jar $@