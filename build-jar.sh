#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

version=$( cat release_version.txt | sed 's/v//' )
pushd $DIR > /dev/null
    gradle clean build shadowJar -DmainClass=statalign.StatAlign -DbaseName=StatAlign
    cp build/libs/StatAlign-$version-all.jar ./StatAlign.jar
popd > /dev/null
