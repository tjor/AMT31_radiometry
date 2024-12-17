#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Copyright 2021 Florida State University

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

###########################################################################
# Calculate true winds from vessel speed, course and relative wind
#
# These routines will compute meteorological true winds (direction from
# which wind is blowing, relative to true north; and speed relative to
# the fixed earth).
#
# Created: 12/17/96
# Comments updated: 10/01/97
# Last updated:  4/12/2921
# Developed by: Shawn R. Smith and Mark A. Bourassa
# Programmed by: Mylene Remigio
# Converted to Matlab from C true wind computation code truewind.c .
# Converted to python from Matlab with SMOP 0.1 
# Direct questions to:  samos@coaps.fsu.edu
#
# 9/30/2014 : If the true wind has speed and its coming from the north
#	      then its direction should be 360deg. The problem was fixed
#	      in which the program's output showed 0deg instead of 360deg.
# 8/08/2017 : Converted from Matlab to python by David Pablo Cohn
#             (david.cohn@gmail.com); verified on python 2.7-3.5
# 4/12/2021 : Added a check to ensure the input lists to truewinds() have the
#             same length. Change made by Homer McMillan (hmcmillan@coaps.fsu.edu)

from math import pi, cos, sin, atan2, sqrt

DEFAULT_ZRL = 0.0  # clockwise angle between bow and anemometer reference line
DEFAULT_MISSING_VALUES = [-1111.0, # missing val for course_over_ground
                          -9999.0, # missing val for speed_over_ground
                          1111.0,  # missing val for wind_dir
                          9999.0,  # missing val for wind_speed
                          5555.0,  # missing val for heading
                         ]

###########################################################################
# FUNCTION truew() - calculates true winds from vessel speed, course and
# relative wind
#
# INPUTS
#
# crse	real	Course TOWARD WHICH the vessel is moving over
#			the ground. Referenced to true north and the
#                       fixed earth.
# cspd	real	Speed of vessel over the ground. Referenced
#			to the fixed earth.
# hd	real	Heading toward which bow of vessel is pointing.
#			Referenced to true north.
# zlr	real		Zero line reference -- angle between bow and
#			zero line on anemometer.  Direction is clockwise
#			from the bow.  (Use bow=0 degrees as default
#			when reference not known.)
# wdir 	real	Wind direction measured by anemometer,
#			referenced to the ship.
# wspd	real	Wind speed measured by anemometer,referenced to
#			the vessel's frame of reference.
# wmis	real	Five element array containing missing values for
#			crse, cspd, wdir, wspd, and hd. In the output,
#                       the missing value for tdir is identical to the
#                       missing value specified in wmis for wdir.
#                       Similarly, tspd uses the missing value assigned
#			to wmis for wspd.

# *** WDIR MUST BE METEOROLOGICAL (DIRECTION FROM)! CRSE AND CSPD MUST
#     BE RELATIVE TO A FIXED EARTH! ***

# OUTPUT VALUES:

# tdir	real	True wind direction - referenced to true north
#                       and the fixed earth with a direction from which
#			the wind is blowing (meteorological).
# tspd	real	True wind speed - referenced to the fixed earth.
# adir	real	Apparent wind direction (direction measured by
#			wind vane, relative to true north). IS
#                       REFERENCED TO TRUE NORTH & IS DIRECTION FROM
#                       WHICH THE WIND IS BLOWING. Apparent wind
#			direction is the sum of the ship relative wind
#                       direction (measured by wind vane relative to the
#                       bow), the ship's heading, and the zero-line
#			reference angle.  NOTE:  The apparent wind speed
#			has a magnitude equal to the wind speed measured
#		  	by the anemometer.

# DIAGNOSTIC OUTPUT:

# nw	integer		Number of observation times for which tdir and
#                       tspd were calculated (without missing values)
# nwpm	integer		Number of observation times with some values
#  			(crse, cspd, wdir, wspd, hd) missing.  tdir,
#			tspd set to missing value.
# nwam	integer		Number of observation times with all values
#  			(crse, cspd, wdir, wspd, hd) missing. tdir,
#			tspd set to missing value.
# nwf	integer		Number of observation times where the program
#			fails -- at least one of the values (crse, cspd,
#			wdir, wspd, hd) is invalid

def truew(crse=None,
          cspd=None,
          wdir=None,
          zlr=DEFAULT_ZRL,
          hd=None,
          wspd=None,
          wmis=DEFAULT_MISSING_VALUES,
          ):
    # INITIALIZE VARIABLES
    adir = 0
    nw = 0
    nwam = 0
    nwpm = 0
    nwf = 0
    dtor = pi / 180

    # Check course, ship speed, heading, wind direction, and
    # wind speed for valid values (i.e. neither missing nor
    # outside physically acceptable ranges).
    if ((((crse < 0) or (crse > 360)) and (crse != wmis[0])) or
        ((cspd < 0) and (cspd != wmis[1])) or
        (((wdir < 0) or (wdir > 360)) and (wdir != wmis[2])) or
        ((wspd < 0) and (wspd != wmis[3])) or
        (((hd < 0) or (hd > 360)) and (hd != wmis[4]))):
        # When some or all of input data fails range check, true
        # winds are set to missing. Step index for input
        # value(s) being out of range
        nwf += 1
        tdir = wmis[2]
        tspd = wmis[3]
        if ((crse != wmis[0]) and (cspd != wmis[1]) and
            (wdir != wmis[2]) and (wspd != wmis[3]) and (hd != wmis[4])):
            # Step index for all input values being non-missing
            nw += 1
        else:
            if ((crse != wmis[0]) or (cspd != wmis[1]) or
                (wdir != wmis[2]) or (wspd != wmis[3]) or (hd != wmis[4])):
                # Step index for part of input values being missing
                nwpm += 1
            else:
                # Step index for all input values being missing
                nwam += 1
    # When course, ship speed, heading, wind direction, and wind speed
    # are all in range and non-missing, then compute true winds.
    else:
        if ((crse != wmis[0]) and (cspd != wmis[1]) and
            (wdir != wmis[2]) and (wspd != wmis[3]) and (hd != wmis[4])):
            nw += 1
            # Convert from navigational coordinates to
            # angles commonly used in mathematics
            mcrse = 90 - crse
            # Keep the value between 0 and 360 degrees
            if (mcrse <= 0.0):
                mcrse = mcrse + 360.0
            # Check zlr for valid value.  If not valid, set equal to 0
            if ((zlr < 0.0) or (zlr > 360.0)):
                zlr = 0.0
            # Calculate apparent wind direction
            adir = hd + wdir + zlr
            # Keep adir between 0 and 360 degrees
            while adir >= 360.0:
                adir = adir - 360.0

            # Convert from meteorological coordinates to angles
            # commonly used in mathematics
            mwdir = 270.0 - adir
            # Keep mdir between 0 and 360 degrees
            if (mwdir <= 0.0):
                mwdir = mwdir + 360.0
            if (mwdir > 360.0):
                mwdir = mwdir - 360.0
            # Determine the east-west vector component and the
            # north-south vector component of the true wind
            x = wspd*cos(mwdir*dtor) + cspd*cos(mcrse*dtor)
            y = wspd*sin(mwdir*dtor) + cspd*sin(mcrse*dtor)
            # Use the two vector components to calculate the true wind
            # speed
            tspd = sqrt(x*x + y*y)
            calm_flag = 1
            # Determine the angle for the true wind
            if (abs(x) > 1e-05):
                mtdir = (atan2(y,x)) / dtor
            else:
                if (abs(y) > 1e-05):
                    mtdir = 180.0 - (90.0*y) / abs(y)
                else:
                    # The true wind speed is essentially zero: winds
                    # are calm and direction is not well defined
                    mtdir = 270.0
                    calm_flag = 0
            # Convert from the common mathematical angle coordinate to
            # the meteorological wind direction
            tdir = 270.0 - mtdir
            # Make sure that the true wind angle is between
            # 0 and 360 degrees
            while tdir < 0.0:
                tdir = (tdir + 360.0)*calm_flag

            while tdir > 360.0:
                tdir = (tdir - 360.0)*calm_flag

            # Ensure wmo convention for tdir = 360 for win
            # from north and tspd > 0
            if (calm_flag == 1 and (tdir < 0.0001)):
                tdir = 360.0
            x = 0.0
            y = 0.0
        else:
            if ((crse != wmis[0]) or (cspd != wmis[1]) or
                (wdir != wmis[2]) or (wspd != wmis[3]) or (hd != wmis[4])):
                nwpm = nwpm + 1
                tdir = wmis[2]
                tspd = wmis[3]
                # When course, ship speed, apparent direction, and
                # wind speed are all in range but all of these input
                # values are missing, then set true wind direction and
                # speed to missing.
            else:
                nwam += 1
                tdir = wmis[2]
                tspd = wmis[3]

    return (tdir, tspd, adir, nw, nwam, nwpm, nwf)

###########################################################################
###########################################################################
# FUNCTION truewinds() - calculates true winds for a list of inputs
#
# INPUT VALUES:
#
# sel	integer		Sets option for diagnostic output.  There are
#			four settings:
#
#			Option 4:  Calculates true winds from input
#				   arrays with no diagnostic output or
#				   warnings. NOT RECOMMENDED.
#			Option 3:  [DEFAULT] Diagnostic output lists the
#                                  array index and corresponding
#				   variables that either violate the
#    				   range checks or are equal to the
#				   missing value. An additional table
#                                  lists the number of observation times
#                                  with no missing values, some (but not
#				   all)  missing values, and all missing
#				   values; as well as similar totals for
#				   the observation times that fail the
#				   range checks. Range checks identify
#				   negative input values and verify
#				   directions to be between 0 and 360
#				   degrees.
#			Option 2:  In addition to the default
#				   diagnostics (option 3), a table of
#				   all input and output values for
#				   observation times with missing data
#                                  is provided.
#			Option 1:  Full diagnostics -- In addition to
#				   the diagnostics provided by option 2
#				   and 3, a complete data chart is
#                                  output. The table contains input and
#                                  output values for all observation
#                                  times passed to truewind.

# crse	real list	Course TOWARD WHICH the vessel is moving over
#			the ground. Referenced to true north and the
#                       fixed earth.
# cspd	real list	Speed of vessel over the ground. Referenced
#			to the fixed earth.
# hd	real list	Heading toward which bow of vessel is pointing.
#			Referenced to true north.
# zlr	real		Zero line reference -- angle between bow and
#			zero line on anemometer.  Direction is clockwise
#			from the bow.  (Use bow=0 degrees as default
#			when reference not known.)
# wdir 	real list	Wind direction measured by anemometer,
#			referenced to the ship.
# wspd	real list	Wind speed measured by anemometer,referenced to
#			the vessel's frame of reference.
# wmis	real list	Five element array containing missing values for
#			crse, cspd, wdir, wspd, and hd. In the output,
#                       the missing value for tdir is identical to the
#                       missing value specified in wmis for wdir.
#                       Similarly, tspd uses the missing value assigned
#			to wmis for wspd.

# OUTPUT VALUES:

# tdir	real list	True wind direction - referenced to true north
#                       and the fixed earth with a direction from which
#			the wind is blowing (meteorological).
# tspd	real list	True wind speed - referenced to the fixed earth.
# adir	real list	Apparent wind direction (direction measured by
#			wind vane, relative to true north). IS
#                       REFERENCED TO TRUE NORTH & IS DIRECTION FROM
#                       WHICH THE WIND IS BLOWING. Apparent wind
#			direction is the sum of the ship relative wind
#                       direction (measured by wind vane relative to the
#                       bow), the ship's heading, and the zero-line
#			reference angle.  NOTE:  The apparent wind speed
#			has a magnitude equal to the wind speed measured
#		  	by the anemometer.

# DIAGNOSTIC OUTPUT:

# nw	integer		Number of observation times for which tdir and
#                       tspd were calculated (without missing values)
# nwpm	integer		Number of observation times with some values
#  			(crse, cspd, wdir, wspd, hd) missing.  tdir,
#			tspd set to missing value.
# nwam	integer		Number of observation times with all values
#  			(crse, cspd, wdir, wspd, hd) missing. tdir,
#			tspd set to missing value.
# nwf	integer		Number of observation times where the program
#			fails -- at least one of the values (crse, cspd,
#			wdir, wspd, hd) is invalid

###########################################################################
###########################################################################
def truewinds(sel=None,crse=None,cspd=None,
              wdir=None,zlr=None,hd=None,wspd=None,
              wmis=None):
    # INITIALIZE VARIABLES
    tdir = []
    tspd = []
    adir = []
    nw = 0
    nwam = 0
    nwpm = 0
    nwf = 0

    if not all(len(l) == len(crse) for l in (cspd, wdir, hd, wspd)):
        raise Exception('All input lists must have the same length!')

    for i in range(len(crse)):
        (tdir_tmp,
         tspd_tmp,
         adir_tmp,
         nw_count,
         nwam_count,
         nwpm_count,
         nwf_count) = truew(crse=crse[i],
                            cspd=cspd[i],
                            wdir=wdir[i],
                            zlr=zlr,
                            hd=hd[i],
                            wspd=wspd[i],
                            wmis=wmis)
        tdir.append(tdir_tmp)
        tspd.append(tspd_tmp)
        adir.append(adir_tmp)
        nw += nw_count
        nwam += nwam_count
        nwpm += nwpm_count
        nwf += nwf_count

    #   OUTPUT SELCTION PROCESS
    if 1 == (sel):
        full(crse,cspd,wdir,zlr,hd,adir,wspd,tdir,tspd)
        missing_values(crse,cspd,wdir,hd,wspd,tdir,tspd,wmis)
        truerr(crse,cspd,hd,wdir,wspd,wmis,nw,nwpm,nwam,nwf)
    elif 2 == (sel):
        missing_values(crse,cspd,wdir,hd,wspd,tdir,tspd,wmis)
        truerr(crse,cspd,hd,wdir,wspd,wmis,nw,nwpm,nwam,nwf)
    elif 3 == (sel):
        truerr(crse,cspd,hd,wdir,wspd,wmis,nw,nwpm,nwam,nwf)
    else:
        print('Selection not valid. Using selection #3 by default. ')
        truerr(crse,cspd,hd,wdir,wspd,wmis,nw,nwpm,nwam,nwf)

    return (tdir, tspd, adir, nw, nwam, nwpm, nwf)

if __name__ == '__main__':
    pass

# NOTE: definitions below are for diagnostic purposes only, and have not
# been prettified or adapted from the mechanical Matlab-to-Python conversion.

###########################################################################
###########################################################################
# **********************************************************************
#                      OUTPUT SUBROUTINES
# **********************************************************************

# Function:  FULL
#  Purpose:  Produces a complete data table with all values.
#            Accessed only when selection #1 is chosen.

def full(crse=None,cspd=None,wdir=None,zlr=None,hd=None,adir=None,wspd=None,tdir=None,tspd=None,*args,**kwargs):
    print('\n------------------------------------------------------------------------------------\n')
    print('                                   FULL TABLE')
    print('                                  ************')
    print('  index  course  sspeed  windir  zeroln  shiphd |  appspd |  appdir  trudir  truspd')
    for j in range(len(tdir)):
        print('%7d %7.1f %7.1f %7.1f %7.1f %7.1f | %7.1f | %7.1f %7.1f %7.1f' % ((j+1),crse[j],cspd[j],wdir[j],zlr,hd[j],wspd[j],adir[j],tdir[j],tspd[j]))

    print('\n                   NOTE:  Wind speed measured by anemometer is identical')
    print('                          to apparent wind speed (appspd).')
    print('\n------------------------------------------------------------------------------------\n')
    return

if __name__ == '__main__':
    pass

    # **********************************************************************

#    Function:  MISSING_VALUES
# Purpose:  Produces a data table of the data with missing values.
#           Accessed when selection #1 or #2 is chosen.
def missing_values(crse=None,cspd=None,wdir=None,hd=None,wspd=None,tdir=None,tspd=None,wmis=None,*args,**kwargs):
    print('                               MISSING DATA TABLE')
    print('                              ********************')
    print('          index  course  sspeed  windir  shiphd  appspd  trudir  truspd')
    for j in range(len(tdir)):
        if ((crse[j] != wmis[0]) and (cspd[j] != wmis[1]) and (wdir[j] != wmis[2]) and (wspd[j] != wmis[3]) and (hd[j] != wmis[4])):
            continue
        else:
            print('        %7d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f' % ((j+1),crse[j],cspd[j],wdir[j],hd[j],wspd[j],tdir[j],tspd[j]))

    print('\n------------------------------------------------------------------------------------\n')
    return

if __name__ == '__main__':
    pass

# **********************************************************************
#    Function:  TRUERR
#    Purpose:  List of where range tests fail and where values are
#               invalid.  Also prints out number of records which are
#               complete, incomplete partially, incomplete entirely, and
#               where range tests fail.  Accessed when selection #1, #2,
#               #3, or the default is chosen.
def truerr(crse=None,cspd=None,hd=None,wdir=None,wspd=None,wmis=None,nw=None,nwpm=None,nwam=None,nwf=None,*args,**kwargs):
    print('                               TRUEWINDS ERRORS')
    print('                              ******************')
    for i in range(len(crse)):
        if (((crse[i] < 0) or (crse[i] > 360)) and (crse[i] != wmis[0])):
            print('        Truewinds range test failed.  Course value #%d invalid.' % (i+1))
        if ((cspd[i] < 0) and (cspd[i] != wmis[1])):
            print('        Truewinds range test failed.  Vessel speed value #%d invalid.' % (i+1))
        if (((wdir[i] < 0) or (wdir[i] > 360)) and (wdir[i] != wmis[2])):
            print('        Truewinds range test failed.  Wind direction value #%d invalid.' %(i+1))
        if ((wspd[i] < 0) and (wspd[i] != wmis[3])):
            print('        Truewinds range test failed.  Wind speed value #%d invalid.' % (i+1))
        if (((hd[i] < 0) or (hd[i] > 360)) and (hd[i] != wmis[4])):
            print('        Truewinds range test failed.  Ship heading value #%d invalid.' % (i+1))

    print('')
    for i in range(len(crse)):
        if (crse[i] == wmis[0]):
            print('        Truewinds data test:  Course value #%d missing.' %(i+1))
        if (cspd[i] == wmis[1]):
            print('        Truewinds data test:  Vessel speed value #%d missing.' % (i+1))
        if (wdir[i] == wmis[2]):
            print('        Truewinds data test:  Wind direction value #%d missing.' % (i+1))
        if (wspd[i] == wmis[3]):
            print('        Truewinds data test:  Wind speed value #%d missing.' % (i+1))
        if (hd[i] == wmis[4]):
            print('        Truewinds data test:  Ship heading value #%d missing.' % (i+1))

    print('\n------------------------------------------------------------------------------------\n')
    print('                                 DATA REVIEW')
    print('                                *************')
    print('                            no data missing = %4d' % nw)
    print('                       part of data missing = %4d' % nwpm)
    print('                           all data missing = %4d' % nwam)
    print('                         failed range tests = %4d' % nwf)
    return

if __name__ == '__main__':
    pass