FROM ubuntu:22.04

RUN echo "Updating Ubuntu"
RUN apt-get update && apt-get upgrade -y

RUN echo "Installing dependencies..."
RUN apt install -y \
            ccache \
            clang \
            clang-format \
            clang-tidy \
            cppcheck \
            curl \
            doxygen \
            gcc \
            git \
            graphviz \
            make \
            ninja-build \
            python3 \
            python3-pip \
            tar \
            unzip \
            vim

RUN echo "Installing dependencies not found in the package repos..."

RUN apt install -y wget tar build-essential libssl-dev && \
            wget https://github.com/Kitware/CMake/releases/download/v3.25.1/cmake-3.25.1.tar.gz && \
            tar -zxvf cmake-3.25.1.tar.gz && \
            cd cmake-3.25.1 && \
            ./bootstrap && \
            make && \
            make install

# Disabled pthread support for GTest due to linking errors
RUN git clone https://github.com/google/googletest.git --branch release-1.12.1 && \
        cd googletest && \
        cmake -Bbuild -Dgtest_disable_pthreads=1 && \
        cmake --build build --config Release && \
        cmake --build build --target install --config Release
