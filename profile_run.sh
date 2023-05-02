#!bin/bash
# profile_run.sh

python3 -m cProfile -o run_profile.prof test_run.py

snakeviz run_profile.prof