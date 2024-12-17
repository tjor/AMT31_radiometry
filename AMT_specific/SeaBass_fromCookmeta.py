#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 17:36:57 2024

@author: tjor
"""

 ## The relative azimuth convention adopted in HyperCP (and expected in ancillary inputs)
 ## defines relAz as the angle between the sensor viewing angle (from the sensor to the target)
##  and the solar azimuth angle (from the sensor to the sun). 

# lat, lon, ship heading, ship speed, relative sensor azimuth,

# general imports
import pandas as pd
import numpy as np
import xarray as xr
import datetime
import scipy
from scipy import signal as sg
import glob

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as cl

import ephem
import math 
import pyproj
from truewind import truew
from truewind import truewinds


def calc_true_wind(heading, gps_speed, app_wind_dir, app_wind_speed):
     
    print('calculating true wind')
    # calculates true_wind between headings and wind directions 
    # approximate course of ship by heading 
    
    true_wind_angle = []
    true_wind_speed = []
    zlr = 0 # zero reference line
    for i in range(len(heading)):
            
            crse_i = heading[i] # approximate course of ship via heading
            cspd_i = gps_speed[i]
            wdir_i = app_wind_dir[i]
            wspd_i = app_wind_speed[i]
            hd_i = heading[i]
           
            true_wind_vec = truew(crse_i, cspd_i, wdir_i, zlr, hd_i , wspd_i)
            true_wind_angle.append(true_wind_vec[0])
            true_wind_speed.append(true_wind_vec[1])

    return true_wind_angle, true_wind_speed


def calc_solar_angles(time, lat, lon):
    'solar elevatuon and azimuth angle function from ephem library. Acts on timestamp, and lat-lon vectors'''
    #
    print('calculating solar angles')
    
    solar_elevation = np.nan*np.ones(len(time))
    solar_azimuth = np.nan*np.ones(len(time))
    for i in range(len(time)):    
       
        obs = ephem.Observer()
        sun = ephem.Sun()
        obs.date = time[i]
        obs.lat, obs.lon = str(lat[i]), str(lon[i])
        sun.compute(obs)
        solar_elevation[i] = (sun.alt * 180. / np.pi)
        solar_azimuth[i] =  (sun.az* 180. / np.pi)
   
    return solar_elevation, solar_azimuth



def export_2_seabass(header, amt2csv, fnout):
    print('writing SeaBASS file...')

    with open(fnout, 'w') as ict:
        # Write the header lines, including the index variable for
        # the last one if you're letting Pandas produce that for you.
        # (see above).
        for key in header.keys():
            ict.write(key + header[key] + "\n")

        # Just write the data frame to the file object instead of
        # to a filename. Pandas will do the right thing and realize
        # it's already been opened.
        amt2csv.to_csv(ict, header = False,
                            index = False,
                            na_rep = header['/missing='])

    print('...done')

    return fnout


def plot_geometry():
  
        plt.figure(figsize=(14,20))
        plt.rcParams.update({'font.size': 16})
        plt.subplot(3,1,1)
        plt.plot_date(gps.time, gps.heading, ms=2, alpha= 0.5, label ='Ship heading')  
        plt.plot_date(gps.time, gps.heading - 60, ms=2, label='HSAS azimuth')
        plt.plot_date(gps.time, solar_azimuth, ms=2,label='Solar azimuth')
       # plt.ylim(0,360)
        plt.ylabel('Bearing [deg clockwise]')
        plt.legend(markerscale=3)
        plt.xticks(rotation=45)
        
        plt.subplot(3,1,2)
        plt.plot_date(gps.time, abs(solar_azimuth - (gps.heading - 60)), ms=2,alpha= 0.5,label='|HSAS azimuth - Solar aziumth|')
        plt.plot_date(gps.time, np.ones(len(gps_time))*150, ms=1,alpha=0.4, color='green', label ='150 deg.')
        plt.plot_date(gps.time, np.ones(len(gps_time))*135, ms=1,alpha=0.4, color='pink', label ='135 deg.')
        plt.plot_date(gps.time, np.ones(len(gps.time))*90, ms=1, alpha=0.4, color='gray', label = '90 deg.')
        plt.ylabel('Angle [deg]')
        plt.legend(markerscale=3 )
        plt.xticks(rotation=45)
         
        plt.subplot(3,1,3)
        plt.plot_date(gps.time, solar_elevation, ms=2,label='Solar elevation')
        plt.plot_date(gps.time, np.ones(len(gps.time))*70, ms=1, alpha=0.4, color='green', label = '70 deg.')
        plt.plot_date(gps.time, np.ones(len(gps_time))*30, ms=1, alpha=0.4, color='pink', label ='30 deg.')
        plt.plot_date(gps.time, np.ones(len(gps.time))*20, ms=1, alpha=0.4, color='gray', label = '20 deg.')
        plt.ylabel('Angle [deg]')
        plt.xlabel('Time [MM-DD HH]')
        plt.legend(markerscale=3 )
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        
        filename  =  fig_dir + '/'   + 'JC_geometry_' + yyyymmdd + '.png'
        plt.savefig(filename,dpi=300)
        
        
        return
    
def plot_geometry_east():
    
            heading = 45
      
            plt.figure(figsize=(14,20))
            plt.rcParams.update({'font.size': 16})
            plt.subplot(3,1,1)
            plt.plot_date(gps.time, np.ones(len(gps.time))*heading, ms=2, alpha= 0.5, label ='Ship heading')  
            plt.plot_date(gps.time, np.ones(len(gps.time))*(heading-60), ms=2, label='HSAS azimuth')
            plt.plot_date(gps.time, solar_azimuth, ms=2,label='Solar azimuth')
           # plt.ylim(0,360)
            plt.ylabel('Bearing [deg clockwise]')
            plt.legend(markerscale=3 )
            
            plt.subplot(3,1,2)
            plt.plot_date(gps.time, abs(solar_azimuth - np.ones(len(gps.time))*(heading - 60)), ms=2,alpha= 0.5,label='|HSAS azimuth - Solar aziumth|')
            plt.plot_date(gps.time, np.ones(len(gps_time))*270, ms=1,alpha=0.4, color='orange', label ='270 (-90) deg.')
            plt.plot_date(gps.time, np.ones(len(gps_time))*225, ms=1,alpha=0.4, color='pink', label ='225 (-235) deg.')
            plt.plot_date(gps.time, np.ones(len(gps.time))*210, ms=1, alpha=0.4, color='gray', label = '210 (-150) deg.')          
            plt.plot_date(gps.time, np.ones(len(gps_time))*150, ms=1,alpha=0.4, color='orange', label ='150 deg.')
            plt.plot_date(gps.time, np.ones(len(gps_time))*135, ms=1,alpha=0.4, color='pink', label ='135 deg.')
            plt.plot_date(gps.time, np.ones(len(gps.time))*90, ms=1, alpha=0.4, color='gray', label = '90 deg.')
            plt.ylabel('Angle [deg]')
            plt.legend(markerscale=3 )
             
            plt.subplot(3,1,3)
            plt.plot_date(gps.time, solar_elevation, ms=2,label='Solar elevation')
            plt.plot_date(gps.time, np.ones(len(gps.time))*70, ms=1, alpha=0.4, color='green', label = '70 deg.')
            plt.plot_date(gps.time, np.ones(len(gps_time))*30, ms=1, alpha=0.4, color='pink', label ='30 deg.')
            plt.plot_date(gps.time, np.ones(len(gps.time))*20, ms=1, alpha=0.4, color='gray', label = '20 deg.')
            plt.ylabel('Angle [deg]')
            plt.xlabel('Time [MM-DD HH]')
            plt.legend(markerscale=3 )
            plt.xticks(rotation=45)
            
            filename  =  fig_dir + '/'   + 'JC_geometry_NorthEast_hdg' + yyyymmdd + '.png'
            plt.savefig(filename,dpi=300)
            
            return
    
def plot_windsonic():
    
        plt.figure(figsize=(18,20))
        plt.rcParams.update({'font.size': 16})
        plt.subplot(2,1,1)
        plt.plot_date(gps.time, gps.heading, ms=2, alpha= 0.5, label ='Ship heading')  
        plt.plot_date(gps.time, app_wind_dir, ms=2, alpha= 0.5, label='Apparent wind direction (w-sonic)')
        plt.plot_date(gps.time, true_wind_dir, ms=2,alpha= 0.5,label='True wind direction (w-sonic)')
        plt.ylim(0,360)
        plt.ylabel('Bearing [deg. clockwise]')
        plt.legend(markerscale=3 )
        
        plt.subplot(2,1,2)
        plt.plot_date(gps.time, gps_speed, ms=2,alpha= 0.5,label ='Ship speed') 
        plt.plot_date(gps.time, app_wind_speed, ms=2,alpha= 0.5,label='Apparent wind speed (w-sonic)')
        plt.plot_date(gps.time, true_wind_speed, ms=2,alpha= 0.5,label='True wind speed (w-sonic)')
        plt.ylabel('Speed [m/s]')
        plt.xlabel('Time [MM-DD HH]')
        plt.legend(markerscale=3 )
        plt.xticks(rotation=45)
        
        filename  =  fig_dir + '/'   + 'JC_windsonic_' + yyyymmdd + '.png'
        plt.savefig(filename,dpi=300)
     
        return
    
def plot_windsurf():
    
    plt.figure(figsize=(18,20))
    plt.rcParams.update({'font.size': 16})
    plt.subplot(3,1,1)
    plt.plot_date(gps.time, gps.heading, ms=2, alpha= 0.5, label ='Ship heading')  
    plt.plot_date(gps.time, app_wind_dir_surf, ms=2, alpha= 0.5, label='Apparent wind direction (surf)')
    plt.plot_date(gps.time, true_wind_dir_surf, ms=2,alpha= 0.5,label='True wind direction (surf)')
    plt.ylim(0,360)
    plt.ylabel('Bearing [deg. clockwise]')
    plt.legend(markerscale=3 )
    
    plt.subplot(3,1,2)
    plt.plot_date(gps.time, gps_speed, ms=2,alpha= 0.5,label ='Ship speed') 
    plt.plot_date(gps.time, app_wind_speed_surf, ms=2,alpha= 0.5,label='Apparent wind speed (surf)')
    plt.plot_date(gps.time, true_wind_speed_surf, ms=2,alpha= 0.5,label='True wind speed (surf)')
    plt.ylabel('Speed [m/s]')
    plt.xlabel('Time [MM-DD HH]')
    plt.legend(markerscale=3 )
    plt.xticks(rotation=45)
        
    filename  =  fig_dir + '/'  + 'JC_windsurf_' + yyyymmdd + '.png'
    plt.savefig(filename,dpi=300)
    
    return

def plot_att():

    plt.figure()
    plt.rcParams.update({'font.size': 16})
    plt.subplot(2,1,1)
    plt.plot_date(gps.time, roll, ms=2,alpha= 0.5,label ='Roll') # assumes, 
    plt.plot_date(gps.time, pitch, ms=2,alpha= 0.5,label='Pitch')
    plt.ylabel('[deg.]')
    plt.xlabel('Time [MM-DD HH]')
    plt.legend(markerscale=3)
 
    plt.subplot(2,1,2)
    plt.plot_date(gps.time, tilt, ms=2,alpha= 0.5,label ='Tilt') # assumes, 
    plt.ylabel('[deg.]')
    plt.xlabel('Time [MM-DD HH]')
    plt.legend(markerscale=3)
    plt.xticks(rotation=45)
    
    filename  =  fig_dir + '/'  + 'JC_Att_' + yyyymmdd + '.png'
    plt.savefig(filename,dpi=300)

    return



if __name__ == '__main__':
    
   
    
    # directry to ship net CDF files
    nc_dir = '/mnt/d/AMT31/Optics_all/ShipsData/NetCDF/'
    fig_dir = '/mnt/d/AMT31/Optics_all/Data/HSAS/HSAS_anc/Figs/'
    
    # select data hours to process
    first_hour = 9
    last_hour= 19
    
    index = -11 #processes previous day of data with -1
    
    ########################################
    # 1.  gps data
    ############################################
    
    gps_dir = nc_dir + 'GPS/'
    gps_files = glob.glob(gps_dir + '*position-POSMV_GPS*')
    gps = xr.open_dataset(gps_files[index])
    gps = gps.load()
    gps = gps.interp({'time':gps.time[first_hour*3600:last_hour*3600]}) # reduce to daylight hours
    
    lat = gps['lat'].values  
    lon = gps['long'].values
    #heading = gps['heading'].values
    heading = 90*np.ones(len(gps['heading'].values))
    gps_speed =  gps['gndspeed']*0.514444 # speed m/s (converted from knots)  
  

    # convert gps time to seabass time fields
    gps_time = [] 
    year = [] 
    month = [] 
    day = []
    hour = []
    minute = []
    second =[]
    print('converting to seabass time fields')
    for i in range(len(gps.time.values)):
  
        timestamp_i = pd.to_datetime(gps.time.values[i])
        gps_time.append(timestamp_i)
        year.append(timestamp_i.year)
        month.append(timestamp_i.month)
        day.append(timestamp_i.day)
        hour.append(timestamp_i.hour)
        minute.append(timestamp_i.minute)
        second.append(timestamp_i.second)
    year = np.array(year)
    month = np.array(month)
    day = np.array(day)
    hour = np.array(hour)
    minute = np.array(minute)
    second = np.array(second)

    # calculate solar geometry
    solar_elevation, solar_azimuth = calc_solar_angles(gps_time, lat, lon)
    HSAS_offset = 60 # mounting on hexagonal scaffold
    RelAz = solar_azimuth - (gps.heading - HSAS_offset) 
    


    ########################################
    # 2. Wind data
    #######################################

    # wind data - sonic (digital)
    wind_dir = nc_dir + 'WINDSONIC/'
    wind_files = glob.glob(wind_dir + '*WINDSONIC*')
    wind = xr.open_dataset(wind_files[index])
    wind = wind.load()
    wind = wind.interp({'time':gps_time}) # interpolate wind to gps time

    app_wind_speed = wind.speed.values*0.514444  # apparent windspeed in m/s (converted from knots)
    app_wind_dir = wind.direct.values
    true_wind_dir, true_wind_speed = calc_true_wind(heading, gps_speed, app_wind_dir, app_wind_speed) # calc
    

    # wind data surf (analog)
    wind_dir = nc_dir + 'SURFMETV3/'
    wind_surf_files = glob.glob(wind_dir + '*MET-SURFMET*')
    wind_surf = xr.open_dataset(wind_surf_files[index])
    wind_surf = wind_surf.load()
    wind_surf = wind_surf.interp({'time':gps_time}) # interpolate wind to gps time
    
    app_wind_speed_surf = wind_surf.speed.values*0.514444  # apparent windspeed in m/s (converted from knots)
    app_wind_dir_surf = wind_surf.direct.values
    true_wind_dir_surf, true_wind_speed_surf = calc_true_wind(heading, gps_speed, app_wind_dir_surf, app_wind_speed_surf) # calculate true wind speed


    print('Windspeed sonic minus windspeed surf = ')
    print(str(np.nanmean(np.array(true_wind_speed_surf)-np.array(true_wind_speed))))
    print('m/s')


    ###################
    # 3. Att - pitch roll tilt
    ###################
     
    att_dir = nc_dir + 'ATT/'
    att_files = glob.glob(att_dir + '*shipattitude-POSMV_ATT*')
    # att_files = glob.glob(att_dir + '*psxn-Seapath330*');
    att = xr.open_dataset(att_files[index]);
    att = att.load()    
    att = att.interp({'time':gps_time}) # interpolate app to gps time
    roll = att['roll'].values
    pitch = att['pitch'].values
    tilt = (180/np.pi)*np.arctan(np.sqrt((roll*np.pi/180)**2 + (pitch*np.pi/180)**2)) # small angle approx

    
    ################
    # 4.  Structure SeaBASS file
    ###############

    # seabass header
    yyyymmdd = str(year[0]).zfill(2) + str(month[0]).zfill(2) + str(day[0]).zfill(2)
    start_time = str(gps_time[0].hour).zfill(2) + ':' + str(minute[0]).zfill(2) + ':' +  str(minute[0]).zfill(2) + '[GMT]' 
    end_time = str(gps_time[-1].hour).zfill(2) + ':' + str(minute[-1]).zfill(2) + ':' +  str(minute[-1]).zfill(2) + '[GMT]' 
    min_lat = str(min(lat))[0:7]
    max_lat = str(max(lat))[0:7]
    min_lon = str(min(lon))[0:7] 
    max_lon = str(max(lon))[0:7]
    
    station = -9999*np.ones(len(gps.time)).astype(int)
    
    sb_header = {
    "/begin_header": "",
    "/investigators=": "Tom_Jordan,Gavin_Tilstone",
    "/affiliations=": "Plymouth_Marine_Laboratory",
    "/contact=": "tjor@pml.ac.uk,ghti@pml.ac.uk",
    "/experiment=": "AMT_HSAS",
    "/cruise=": "AMT_31",
    "/station=": "NA",
    "/data_file_name=": "JamesCookAnc_" + yyyymmdd,
    "/documents=": "NA", 
    "/calibration_files=": "NA", #
    "/data_type=": "above_water",
    "/data_status=": "final",
    "/start_date=":  yyyymmdd,
    "/end_date=":  yyyymmdd,
    "/start_time=":  start_time,
    "/end_time=": end_time,
    "/north_latitude=":  max_lat,
    "/south_latitude=":  min_lat,
    "/east_longitude=":  max_lon,
    "/west_longitude=":  min_lon,
    "/missing=": "-9999",
    "/delimiter=": "comma",
    "/fields=":"station,year,month,day,hour,minute,second,lat,lon,wind,wdir,RelAz,heading,tilt,roll,pitch",
    "/units=":"none,yyyy,mo,dd,hh,mm,ss,degrees,degrees,m/s,degrees,degrees,degrees,degrees,degrees,degrees",
    "/end_header": "",
          }
    
    # seabass data
    sb_df = pd.DataFrame()
    sb_df['station'] = station
    sb_df['year'] = year
    sb_df['month'] = month
    sb_df['day'] = day
    sb_df['hour'] = hour
    sb_df['minute'] = minute
    sb_df['second'] = second
    sb_df['lat'] = lat
    sb_df['lon'] = lon
    sb_df['wind'] = true_wind_speed
    sb_df['wdir'] = true_wind_dir
    sb_df['RelAz'] = RelAz
    sb_df['heading'] = heading
    sb_df['tilt'] = tilt
    sb_df['roll'] = roll
    sb_df['pitch'] = pitch
    
    filename = '/mnt/d/AMT31/Optics_all/Data/HSAS/HSAS_anc/' + "JamesCookAnc_" + yyyymmdd + '.sb'
    export_2_seabass(sb_header, sb_df, filename)
 
    # plot  functions
    plot_geometry()
    plot_windsonic()
    plot_windsurf()
    plot_att()
