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

if [[ $PLAT =~ .*_x86_64 ]]; then
    declare -a PythonVersions=("cp36-cp36m" "cp37-cp37m" "cp38-cp38")
else
    declare -a PythonVersions=("cp36-cp36m" "cp37-cp37m")
fi

for val in ${PythonVersions[@]}; do
    /opt/python/$val/bin/pip install -r /src/requirements.txt
    /opt/python/$val/bin/pip wheel /src/ --no-deps -w /wheelhouse/
done


if [[ $PLAT =~ .*_x86_64 ]]; then
    suffix=*_x86_64.whl
else
    suffix=*_i686.whl
fi

for whl in /wheelhouse/$suffix; do
    repair_wheel "$whl"
done
