#!/bin/sh
#
# tabcolmax.sh - indicate the maximum value in a LaTeX table column
#                by enclosing it in \colmax{}
#
# Alex Stivala, originally 2006
# modified for real numbers not integers October 2008,
# works on room0219pos09.cs.mu.oz.au (maybe not charikar.cs.mu.oz.au,
# note strange problem with differnece in handling \\ and \\\\ etc.) even
# though both Linux
#
# Usage: tablcolmax <inputfilename> <col1> <col2> ...
#
# The <inputfilename> is the file the table is read from (must be able
# to be read twice, so we can't use stdin).
#
# <col1> <col2> etc are the column (field) numbers (from 1) to process. 
# The field separator is '&'
#
# Define something like:
#
# \newcommand{\colmax}[1]{\textbf{1}} % indicate maximum value in a table column
#
# in the LaTeX file to make the maximum value in boldface.
#
# $Id: tabcolmax.sh 1961 2008-10-07 06:08:06Z astivala $
#



if [ $# -lt 2 ]; then
    echo "Usage: $0 <intputfilename> <col1> <col2> ..." >&2
    exit 1
fi

infile=$1
shift
col_list=$@

delim='&'

tempfile=/var/tmp/tabcolmax.$$.tmp
cat /dev/null > $tempfile

# we are assuming every row has the same numbef of fields
numfields=`head -1 $infile | sed "s/[^${delim}]//g" | wc -c`

for col_num in $col_list ; do
    # we assume that the field we are processing is numeric
    maxval=`cut -d $delim -f $col_num < $infile | sort -n -r | head -1`
    if [ $col_num -eq $numfields ]; then
        # last field is messy: has whitesapce then \\\ so fix it
        maxval=`echo "$maxval" | cut -f1 -d' '`
    fi
    /bin/echo -e "$col_num\t$maxval" >> $tempfile
done


# tricky: sed needed to replace \ with \\ since read removes single \
#sed 's/\\/\\\\/g' $infile | while read line ; do
cat $infile | while read line ; do
    i=1
    while [ $i -le $numfields ]; do
        done_field=0
        field=`echo "$line" | cut -d $delim -f $i`
#        if [ $i -eq $numfields ]; then
#            # last field is messy: has whitesapce then \\\ so fix it
#            field=`echo "$field" | cut -f1 -d' '`
#        fi
        for col_num in $col_list ; do
            if [ $col_num -eq $i ]; then
                maxline=`grep "^${col_num}[ ]*" $tempfile`
                if [ $? -eq 0 ]; then
                    maxval=`echo "$maxline" | cut -f2`
                    iseq=`echo "$field == $maxval" | bc`
                    if [ $iseq -ne 0 ]; then
                        /bin/echo -En "\colmax{$field}"
                    else
                        echo -n "$field"
                    fi
                    done_field=1
                fi
            fi
        done
        if [ $done_field -eq 0 ]; then
            echo -n "$field"
        fi
        if [ $i -eq $numfields ]; then
            echo \\\\
        else
            echo -n $delim
        fi
        i=`expr $i + 1`
    done
done 

rm $tempfile

