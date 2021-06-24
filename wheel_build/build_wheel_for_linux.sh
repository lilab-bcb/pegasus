#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$whl" --plat "$PLAT" -w /wheelhouse/
        rm "$whl"
    fi
}

declare -a PythonVersions=("cp37-cp37m" "cp38-cp38" "cp39-cp39")

for val in ${PythonVersions[@]}; do
    /opt/python/$val/bin/pip install -r /src/requirements.txt
    /opt/python/$val/bin/pip wheel /src/ --no-deps -w /wheelhouse/
done

suffix=*_x86_64.whl

for whl in /wheelhouse/$suffix; do
    repair_wheel "$whl"
done
