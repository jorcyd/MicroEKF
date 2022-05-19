#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License

SRC = ../../src

all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.c $(SRC)/tiny_ekf.c $(SRC)/tiny_ekf.h
	gcc -Wall -I. -I$(SRC) -o gps_ekf gps_ekf.c $(SRC)/tiny_ekf.c -lm

edit:
	vim gps_ekf.c

clean:
	rm -f gps_ekf *.o *~ ekf.csv

#EMCC for testing.
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 --post-js gps_ekf_download.js -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv -Wno-unused-function
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 --post-js gps_ekf_download.js -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv -Wno-unused-function
#emcc -Wall -O0 -g -Wdouble-promotion -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#clean_emcc:
#	rm -f gps_ekf.html gps_ekf.js gps_ekf.wasm