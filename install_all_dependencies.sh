#!/usr/bin/env bash
# Check if python is available
if (! hash python3 2>/dev/null); then
  echo "You do not have python3 installed, please install it before continuing"
  exit 1
fi

if (! hash Rscript 2>/dev/null); then
  echo "You do not have R installed, please install it before continuing"
  exit 1
fi

if (! hash mafft 2>/dev/null); then
  echo "You do not have mafft installed, please install it before continuing"
  exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Install R dependencies
Rscript "${DIR}"/install_R_dependencies.R
# Install python dependencies
python3 -m pip install -r "${DIR}"/requirements.txt
