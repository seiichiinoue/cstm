# image
FROM ubuntu:14.04

RUN apt-get update && apt-get install -y libboost-all-dev python3-dev git cmake g++ gdb python3-dbg vim

# language environment settings
RUN apt-get -y install language-pack-ja-base language-pack-ja

# env var
ENV LANG ja_JP.UTF-8

WORKDIR /workspace/