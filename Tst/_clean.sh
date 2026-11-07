#!/bin/sh
#delete processing results
find . -type f -name "*_??.raw" -delete
find . -type f -name "*_??_????.raw" -delete
find . -type f -name "*_1??_????.raw" -delete
find . -type f -name "*_2??_????.raw" -delete
find . -mindepth 1 -maxdepth 1 -type d -name '$*' -exec rm -rf {} +