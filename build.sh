#!/bin/bash

if [ -f ".env" ]; then
    . ".env"
else
    echo "WORNING: .env file does not exist."
fi

if [ -d "_build" ]; then
    rm -rf _build
fi

echo -en "Q. Which BUILD_TYPE do you select? (release/debug) : "
read _build_type
case $_build_type in
r*)
    build_type=Release
;;
d*)
    build_type=Debug
;;
*)
    echo "ERROR: '$_build_type' is invalid. expected 'r*' or 'd*'"
    exit 1
;;
esac

echo -e "> '$build_type' is selected"
echo -e ""


cmake -S . -B _build \
    -D CMAKE_BUILD_TYPE=${build_type} \
&& cmake --build _build