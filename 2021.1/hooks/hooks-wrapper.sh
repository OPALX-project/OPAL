#!/bin/bash

if [ -x $0.local ]; then
    $0.local "$@" || exit $?
fi

if [ -x hooks/$(basename $0) ];then
    hooks/$(basename $0) "$@" || exit $?
fi