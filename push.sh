#! /bin/bash

git add --all .
git commit -am "auto"
git push
git push github

cd ~/git_project/netsnap
git add --all .
git commit -am "auto"
git push
git push github

