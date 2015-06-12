#!/usr/bin/env bash

export LD_LIBRARY_PATH=/usr/local/openmm/lib:$LD_LIBRARY_PATH

python "$@"
