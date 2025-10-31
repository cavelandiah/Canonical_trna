#!/usr/bin/env bash

./generate_ss_shape_V2.py g "" NAI 75
./generate_ss_shape_V2.py g "" DMS 75
./generate_ss_shape_V2.py g "compatible" NAI 75
./generate_ss_shape_V2.py g "compatible" DMS 75
./generate_ss_shape_V2.py g "noncompatible" DMS 75
./generate_ss_shape_V2.py g "noncompatible" NAI 75
