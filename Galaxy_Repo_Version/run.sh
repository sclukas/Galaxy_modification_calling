#!/bin/sh

ROOT_DIR="$(dirname $(realpath "$0"))"

docker run \
	--rm \
	-p 8080:80 \
	-v ${ROOT_DIR}/shared:/export \
	galaxy-helm
