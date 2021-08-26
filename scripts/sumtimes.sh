#!/bin/sh
#
# File:    sumtimes.sh
# Author:  Alex Stivala
# Created: August 2008
#
# sumtimes.sh - add up all the elapsed times from /usr/bin/time in input files
#
# Usage: sumtimes.sh [-m|-h] file_list
#
# -m : output only minutes and seconds, no hours
# -h : only output hours and minutes, no seconds
#
# file_list if list of files containing in them /usr/bin/time output.
# assume times are in this sort of format, as generated by /usr/bin/time
# on Linux:
#49576.89user 2907.19system 14:45:02elapsed 98%CPU (0avgtext+0avgdata 0maxresident)k
#0inputs+0outputs (0major+504890689minor)pagefaults 0swaps
# or (for less than an hour):
#2603.51user 217.97system 47:10.00elapsed 99%CPU (0avgtext+0avgdata 0maxresident)k
#0inputs+0outputs (0major+26943522minor)pagefaults 0swaps
#
# Also handles the bash builtin time format e.g.
#
#real    48m8.320s
#user    34m50.430s
#sys     7m20.210s
#
# Output, sum of all the elpased times, is to stdout.
#
# $Id: sumtimes.sh 1878 2008-09-08 07:59:32Z astivala $
#

if [ $# -lt 1 ]; then
    echo "Usage: $0 [-m|-h] file list" >&2
    exit 1
fi

use_hours=1
use_seconds=1

if [ $# -ge 1 ] ; then
    if [ `expr substr $1 1 1` = "-" ]; then
        if [ "$1" = "-m" ]; then
            use_hours=0
        elif [ "$1" = "-h" ]; then
            use_seconds=0
        else
            echo "Usage: $0 [-m|h] file list" >&2
            exit 1
        fi
        shift 1
    fi
fi

total_seconds=0
builtinformat=0
for errfile in $*
do
  grep --text elapsed ${errfile} >/dev/null 2>&1
  if [ $? -eq 0 ]; then
      elapsed=`grep --text elapsed ${errfile} | awk '{print $3}' | sed 's/elapsed//' | tail -1`
      dotindex=`expr index ${elapsed} '.'`
      if [ ${dotindex} -ne 0 ]; then
              # less than an hour
              colonindex=`expr index ${elapsed} ':'`
              colonindex=`expr $colonindex - 1`
              hours=0
              mins=`expr substr ${elapsed} 1 ${colonindex}`
              secindex=`expr $colonindex + 2`
              secs=`expr substr ${elapsed} ${secindex} 2`
      else
              colonindex=`expr index ${elapsed} ':'`
              colonindex=`expr $colonindex - 1`
              hours=`expr substr ${elapsed} 1 ${colonindex}`
              next=`expr $colonindex + 2`
              rest=`expr substr ${elapsed} $next 999`
              colonindex=`expr index ${rest} ':'`
              colonindex=`expr $colonindex - 1`
              mins=`expr substr ${rest} 1 $colonindex`
              next=`expr $colonindex + 2`
              rest=`expr substr ${rest} $next 999`
              secs=$rest
      fi
      total_seconds=`expr $total_seconds + $hours \* 3600 + $mins \* 60 + $secs`
  else
      builtinformat=1
      elapsed=`grep --text real ${errfile} | awk '{print $2}'`
      mindex=`expr index ${elapsed} 'm'`
      mindex=`expr $mindex - 1`
      mins=`expr substr ${elapsed} 1 ${mindex}`
      secindex=`expr $mindex + 2`
      sindex=`expr index ${elapsed} 's'`
      seclen=`expr ${sindex} - ${secindex}`
      secs=`expr substr ${elapsed} ${secindex} ${seclen}`
      
      total_seconds=`echo "$total_seconds + $mins * 60 + $secs" | bc -l`
  fi
done

if [ ${builtinformat} -eq 1 ]; then
      total_seconds=`printf "%.0f" $total_seconds`
fi

if [ $use_hours -eq 1 ]; then
    hours=`expr $total_seconds / 3600`
    mins=`expr $total_seconds - $hours \* 3600`
    mins=`expr $mins / 60`
    rsecs=`expr $total_seconds - $hours \* 3600`
    rsecs=`expr $rsecs - $mins \* 60`
    if [ $use_seconds -eq 1 ]; then
        printf '%d h %02d m %02d s' ${hours} ${mins} ${rsecs}
    else
        if [ $rsecs -ge 30 ]; then
            mins=`expr $mins + 1`
        fi
        printf '%d h %02d m' ${hours} ${mins}
    fi
else
    mins=`expr $total_seconds / 60`
    rsecs=`expr $total_seconds - $mins \* 60`
    printf '%02d m %02d s' ${mins} ${rsecs}
fi
