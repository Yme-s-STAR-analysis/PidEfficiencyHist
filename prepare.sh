#!/bin/bash
if [ ! $# -eq 1 ];then
    echo "Need target path."
else
    target=$1
    echo "target path is $1."
    rm Csubmit.xml
    sed -e "s|TARGET|$1|" .tpl.xml > Csubmit.xml

    mkdir -p $1/report
    mkdir -p $1/csh
    mkdir -p $1/err
    mkdir -p $1/list
    mkdir -p $1/log
    mkdir -p $1/out

    rm -rf  $1/report/* 
    rm -rf  $1/csh/* 
    rm -rf  $1/err/* 
    rm -rf  $1/list/* 
    rm -rf  $1/log/* 
    rm -rf  $1/out/*
    
    rm -rf  libZip*.zip
    rm -rf  libZip*.package
    rm -rf  *dataset
    rm -rf  *session.xml
fi
