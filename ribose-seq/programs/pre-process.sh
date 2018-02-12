#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

trim_galore --gzip--length $minimum --clip_R1 3 -a $adapter - -o $output
