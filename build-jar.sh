#!/bin/bash
thisDir=$(dirname $0) || false

pushd ${thisDir} > /dev/null
    gradle clean build shadowJar
    cp build/libs/StatAlign-3.2-all.jar ./StatAlign.jar
popd > /dev/null
