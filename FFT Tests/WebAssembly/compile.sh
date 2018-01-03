#!/bin/bash
emcc ../C/time_ffts.c -s WASM=1 -O3 -o time_ffts.html
