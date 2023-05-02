#!bin/bash
# profile_post.sh

python3 -m cProfile -o post_profile.prof test_post.py

snakeviz post_profile.prof