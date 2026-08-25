ARG BASE_IMAGE=noaaepic/ufschem-spack-base-ubuntu-gcc-13-dev:latest
FROM ${BASE_IMAGE}

SHELL ["/bin/bash", "-c"]

USER root
RUN apt-get update --yes && apt-get install --yes libkokkos-dev libyaml-cpp-dev

COPY . /opt/catchem/src

RUN source /opt/ufschem/spack-stack/setup.sh && \
    spack env activate ufschem && \
    cd /opt/catchem/src && \
    mkdir build && \
    cd build && \
    cmake .. \
      -DCATCHEM_BUILD_TESTING=ON \
      -DCATCHEM_BUILD_NUOPC=ON \
      -DCATCHEM_BUILD_MUSICA=ON \
      -DCATCHEM_ENABLE_KOKKOS=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/catchem/install \
      -DCMAKE_BUILD_TYPE=Debug && \
    make -j$(nproc) && \
    make install && \
    ctest --output-on-failure

RUN cat <<EOF >> ~/.bashrc
# enable spack env by default
source /opt/ufschem/spack-stack/setup.sh
spack env activate ufschem
EOF

WORKDIR /opt

CMD ["/bin/bash"]
