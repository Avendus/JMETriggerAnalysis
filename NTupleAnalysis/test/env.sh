#!/bin/bash

if [ ! -z ${JMEANA_BASE} ]; then
  echo "environment already set: JMEANA_BASE=${JMEANA_BASE}"
  return
fi

export JMEANA_BASE=${PWD}

export PATH=${JMEANA_BASE}:${PATH}

#export PYTHON27PATH=${PWD}:${PYTHON27PATH}
#export PYTHON27PATH=${PWD}/common:${PYTHON27PATH}
#export PYTHONDONTWRITEBYTECODE=1

export PYTHON3PATH=${PWD}:${PYTHON3PATH}
export PYTHON3PATH=${PWD}/common:${PYTHON3PATH}
export PYTHONDONTWRITEBYTECODE=1
