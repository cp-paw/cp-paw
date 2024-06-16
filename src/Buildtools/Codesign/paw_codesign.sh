#!/bin/bash
export THISDIR=$(pwd)
################################################################################
#  paw_codesign
#
#  Usage: paw_codesign executable
#
#         where executable is the name of the executable to be codesigned
#
#  on macos, the parallel executable of the code needs to be signed
#  to avoid pop-up windows asking whether the code may communicate
#
#  This shell script assumes that there is a signature paw_codesign
#  in the keychain and the file ~/bin/paw-entitlement.xml exists.
#  (see Administration/Howtos/codesign-howto on how to generate them)
#
#  written by Peter Bloechl 13.09.2021
################################################################################

##############################################################
# codesign -vv $(which ppaw_dbg.x) 
# success
#..ppaw_dbg.x: valid on disk
#..ppaw_dbg.x: satisfies its Designated Requirement
# failure 
#.. ppaw_dbg.x: code object is not signed at all
#.. In architecture: x86_64
##############################################################
X=$1
export TMPFILE=$(mktemp /tmp/paw_codesign.XXXXXX)
#-------------------------------------------------------------------------------
#--  verify code signature                                                    --
#-------------------------------------------------------------------------------
codesign -vv $X 2> $TMPFILE 
Y1=$(grep 'valid on disk' $TMPFILE)
Y2=$(grep 'not signed' $TMPFILE)
if [[ -z $Y1 ]]; then
    #---------------------------------------------------------------------------
    # do the code signing
    #---------------------------------------------------------------------------
    echo code-signing $X
#
#   ----------------------------------------------------------------------------
#   -- write a entitlement.xml file to TMPENTITLEMENT
#   ----------------------------------------------------------------------------
    export TMPENTITLEMENT=$(mktemp).xml
    echo '<?xml version="1.0" encoding="UTF-8"?>'            > ${TMPENTITLEMENT}
    echo '<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" '\
      ' "http://www.apple.com/DTDs/PropertyList-1.0.dtd">'  >> ${TMPENTITLEMENT}
    echo '<plist version="1.0">'                            >> ${TMPENTITLEMENT}
    echo '<dict>'                                           >> ${TMPENTITLEMENT}
    echo '  <key>com.apple.security.cs.debugger</key>'      >> ${TMPENTITLEMENT}
    echo '  <true/>'                                        >> ${TMPENTITLEMENT}
    echo '</dict>'                                          >> ${TMPENTITLEMENT}
    echo '</plist>'                                         >> ${TMPENTITLEMENT}
#   ----------------------------------------------------------------------------
#   --  paw_codesign refers to a certificate in the Keychain with that name   --
#   ----------------------------------------------------------------------------
    codesign --entitlements ${TMPENTITLEMENT} -fs paw_codesign $X
    if [[ -f ${TMPENTITLEMENT} ]] ; then rm ${TMPENTITLEMENT} ; fi

    #---------------------------------------------------------------------------
    #-- verify code signature to check whether code-signing succeeded
    #---------------------------------------------------------------------------
    codesign -vv $X 2> $TMPFILE 
    Y1=$(grep 'valid on disk' $TMPFILE)
    if [[ ! -z $Y1 ]]; then
      echo code-signing $X succeeded
    else
      echo code-signing $X failed
    fi
else
  echo no codesigning required for X=$X
fi
rm $TMPFILE
