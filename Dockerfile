# Dockerfile for cnaplotr
# Image: ghcr.io/cchmc-research-mgps/cnaplotr:<tag>

FROM python:3.12-slim-bullseye

LABEL maintainer="Somak Roy<somak.roy@cchmc.org>" \
    function="Docker image with cnaplotr"

ENV NON_ROOT_USER="bioseq"

RUN apt update && \
    adduser --gecos '' --disabled-password $NON_ROOT_USER && \
    apt clean

USER ${NON_ROOT_USER}
WORKDIR /home/${NON_ROOT_USER}

COPY . .

RUN pip3 install -r requirements.txt && \
    pip3 cache purge && \
    rm -rf .vscode .github plots .gitignore Dockerfile LICENSE README.md requirements.txt