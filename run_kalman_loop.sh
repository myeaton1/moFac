#!/bin/bash

for y in {1..48}
do
  for m in {1..12}
  do
    dmatlab12b -r "pca_kalman_oos_loop_it($y, $m, 5)"

    sleep 0
    done

  sleep 0
done