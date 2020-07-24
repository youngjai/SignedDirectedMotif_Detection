#!/bin/sh

git add . && wait
git commit -m "$1" && wait
git push -u motif master && wait
