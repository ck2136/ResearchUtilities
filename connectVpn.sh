#!/bin/bash - 
#===============================================================================
#
#          FILE: connectVpn.sh
# 
#         USAGE: ./connectVpn.sh 
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
# 1. VPN Connect to University Network
#===============================================================================
globalprotect connect -p amc-vpn.ucdenver.edu -u kimchon
