#
# Makefile for TinyEKF GPS example
#
# Copyright (C) 2015 Simon D. Levy
#
# MIT License
# For generating profiling code per line use -g -pg

CC=gcc
CFLAGS=-Wall -Wno-unused-function -Wdouble-promotion 
PROF=-g -pg
SRC=../../src

all: gps_ekf

run: gps_ekf
	./gps_ekf

gps_ekf: gps_ekf.c $(SRC)/tiny_ekf.c $(SRC)/tiny_ekf.h
	$(CC) $(CFLAGS) -I. -I$(SRC) -o gps_ekf gps_ekf.c $(SRC)/tiny_ekf.c -lm
	
prof:	gps_ekf.c $(SRC)/tiny_ekf.c $(SRC)/tiny_ekf.h
	$(CC) $(CFLAGS) $(PROF) -I. -I$(SRC) -o gps_ekf gps_ekf.c $(SRC)/tiny_ekf.c -lm

edit:
	vim gps_ekf.c

clean:
	rm -f gps_ekf gps_ekf.html gps_ekf.js gps_ekf.wasm *.o *.exe *~ ekf.csv gmon.out

#EMCC for testing.
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 --post-js gps_ekf_download.js -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv -Wno-unused-function
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 --post-js gps_ekf_download.js -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv -Wno-unused-function
#emcc -Wall -O0 -g -Wdouble-promotion -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#emcc -Wall -O0 -g -I. -I ../../src -o gps_ekf.html gps_ekf.c ../../src/tiny_ekf.c -lm -s WASM=1 -s EXPORTED_RUNTIME_METHODS=FS --preload-file gps.csv
#clean_emcc:
#	rm -f gps_ekf.html gps_ekf.js gps_ekf.wasm
