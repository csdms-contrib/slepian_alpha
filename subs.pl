#!/bin/csh -f
# Substitute (old) $1 with (new) $2 in $3-$n
#
# Don't need quotes or anything - give path names
# as 
# subs.pl /home/fjsimons /home/backup file
#
# -p loops over input files
# -e evaluates the command line
# /i ignores the case
# /g performs global substitution
# -i does it in-place (can provide back-up extension)
# Using # instead of / allows for use of unmodified / in pathnames
# Use backslash in front of special characters.
# WATCH OUT WHEN USING DOTS AND @ SIGNS!! 
#
# LINKS become files!
#
# EXAMPLE: Use "" to use wildcards!
#
# subs.pl http://geoweb.princeton.edu/people/simons/index.html http://www.frederik.net *-ci.html
#
# subs.pl "@" "-at-" "*.m"
#
# GET RID OF FUNNY STUFF BEFORE COPYRIGHT SIGN
# subs.pl "\\302\\251" "\251" "*.eps"
# GET RID OF FUNNY STUFF BEFORE DEGREE SIGN, SEE ALSO DEGS
# subs.pl "\\302" "" "*.eps"
# PUT SPACE AFTER DEGREE SIGN...
#
# OTHER EXAMPLES:
# subs.pl /brolga/data21/ftp/pub/fjsimons/THESIS/FIGURES/ /home/fjsimons/GifPix/EPS/ "*.tex"
# subs.pl /rosella/data14/newton_tmp/ /home/fjsimons/IFILES/EARTHMODELS/EIGENFUNCTIONS/ "*/results/oldpar.*"
#
# Last modified by fjsimons-at-alum.mit.edu, November 11, 2003

perl -i -pe 's#'$1'#'$2'#ig' $3


