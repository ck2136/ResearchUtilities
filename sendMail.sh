#!/bin/bash - 
#===============================================================================
#
#          FILE: sendMail.sh
# 
#         USAGE: ./sendMail.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: C.K. (), chong.kim@ucdenver.edu
#  ORGANIZATION: 
#       CREATED: 04/11/2019 01:01:16 AM
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

#===============================================================================
# 1. Send Email of Pubmed Abstract
#===============================================================================
mutt -s "Asthma Abstracts" -a /home/ck/Downloads/abstracts/abstractdf.csv -- chong.kim@ucdenver.edu < testmail.txt
