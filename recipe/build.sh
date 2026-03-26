#!/bin/bash
set -ex

$PYTHON -m pip install . --no-deps -vv
$PYTHON -m pip install thirdparty/simple_slurm-0.3.6-py3-none-any.whl
