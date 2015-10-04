#!/bin/bash

echo Timing: ON
printf "\n"

echo Cliquet HIGH
echo ============
for m in 10 20 50 100 200 500 1000
do
	echo m = $m
	python -m timeit -s "from runCliquet import run_cliquet_high" "run_cliquet_high([$m])"
	echo ----------------------------------------
done
printf "\n"

echo Cliquet LOW
echo ===========
for m in 10 20 50 100 200 500 1000
do
	echo m = $m
	python -m timeit -s "from runCliquet import run_cliquet_low" "run_cliquet_low([$m])"
	echo ----------------------------------------
done
printf "\n"

echo Cliquet VAR
echo ===========
for m in 10 20 50 100 200 500 1000
do
	echo m = $m
	python -m timeit -s "from runCliquet import run_cliquet_var" "run_cliquet_var([$m])"
	echo ----------------------------------------
done
printf "\n"

echo Timing: OFF
