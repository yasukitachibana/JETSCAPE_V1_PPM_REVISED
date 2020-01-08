#!/usr/bin/env bash

############################################################################

# Report issues via email to yasuki.tachibana@gmail.com

##############################################################################
 
# Download the Tables for Lattice Equation of State
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "http://yasuki.sakura.ne.jp/PPM_JETSCAPE_V1_EOS/EOS.tar.gz"

tar xvzf EOS.tar.gz

