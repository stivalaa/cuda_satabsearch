#!/bin/sh
#
# mkranktab.sh - make table of ranks of structure that has nonlinear
#                match to query structure from Abyzov & Ilyin (2007)
#                examples.
#
# No formatting at the moment, just outputs the data for summary purposes.
#
# $Id: mkranktab.sh 1790 2008-08-04 02:17:59Z astivala $

for result in d1di6a*.out
do
    echo -n $result '\t'
    getrank.sh d1vi6a_ $result
done

echo

for result in d1ki9a*.out
do
    echo -n $result '\t'
    getrank.sh d1gcaa_ $result
done

