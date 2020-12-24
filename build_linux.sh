docker pull quay.io/pypa/manylinux2014_x86_64
docker run -it --rm -e PLAT=manylinux2014_x86_64 -v `pwd`:/src -v `pwd`/dist:/wheelhouse quay.io/pypa/manylinux2014_x86_64 /bin/sh -c "/src/wheel_build/build_wheel_for_linux.sh"
