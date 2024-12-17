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

import datetime
import scipy

import glob

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as cl

import h5py
import datetime
import subprocess


def run_fcheck(fnout, fcheckdir):
    'Seabass fil checking function - requires pearl script fcheck4 to run'
    print('running fcheck...')

    subprocess.run([fcheckdir + "/fcheck4/fcheck4.pl", fnout])

    print('...done')

    return

def read_hdf_L2(file, grp, var):
    'function to read invii file and open in pandas dataframe format - individual fiule'
    
    f = h5py.File(file, 'r')
    data = f[grp][var][()]
    
    return pd.DataFrame(data)


def read_hdf_L2_daily(hdf_files, grp, var):
    ' function to read hdf files and open in pandas dataframe format'
   
    # print(grp)
    # print(var)
    
    f = read_hdf_L2(hdf_files[0], grp, var) # used to test for spectra
    if len(f.iloc[0]) >100:
        
        for i in range(len(hdf_files)):  # loop over hdf files in day
              print(hdf_files[i])
              df_temp = pd.DataFrame(read_hdf_L2(hdf_files[i], grp, var))
              time = convert_datetime(df_temp['Datetag'].values , df_temp['Timetag2'].values) 
              df_temp.insert(i, 'time', time)
              #df_temp.drop(columns = ['Datetag' ,'Timetag2'])
              if i == 0: # append rows to dataframe
                   df_tot = df_temp
              else:
                   df_tot = pd.concat([df_tot, df_temp])
        
        df_tot = df_tot.drop(columns = ['Datetag' ,'Timetag2'])
        df_tot = df_tot.set_index(['time'])
   
    else:
        
        for i in range(len(hdf_files)):  # loop over hdf files in day
              df_temp = pd.DataFrame(read_hdf_L2(hdf_files[i], grp, var))
              if i == 0: # append rows to dataframe
                   df_tot = df_temp
              else:
                   df_tot = pd.concat([df_tot, df_temp])
        
    return df_tot



def read_anc_daily(hdf_files):
      ' Converts ancillary data into pandas dataframe for a days worth of files'
      
      for i in range(len(hdf_files)):  # loop over hdf files in day
          df_i = []
          for j in range(len(anc_var)): # loop over ancillary fields
              df_temp = pd.DataFrame(read_hdf_L2(hdf_files[i], 'ANCILLARY', anc_var[j]))
              if j == 0:# use first ancillary field to extract datetime index
                  time = convert_datetime(df_temp['Datetag'].values , df_temp['Timetag2'].values) 
                  df_i = pd.DataFrame(index=time)
                  df_i.insert(j, df_temp.keys()[2],df_temp.iloc[:,2].values) # add anc fields as columns
                  df_i.insert(j+1, df_temp.keys()[3],df_temp.iloc[:,3].values)
              else:
                  if len(df_temp.keys()) == 3:
                      df_i.insert(j, df_temp.keys()[2],df_temp.iloc[:,2].values)
                  elif len(df_temp.keys()) == 4:
                      df_i.insert(j, df_temp.keys()[2],df_temp.iloc[:,2].values)
                      df_i.insert(j+1, df_temp.keys()[3],df_temp.iloc[:,3].values)
          if i == 0: # append rows to dataframe
              df_tot = df_i
          else:
              df_tot = pd.concat([df_tot, df_i])
              
      return df_tot

def convert_datetime(date, time):
    ' converts HCP time to date-time format'
    # 'DATETAG'  - Vector of length y, YYYYDOY.0, where DOY is the UTC sequential day of the year
    # 'TIMETAG2' - Vector of length y, HHMMSSUUU UTC (hour, minute, second, millisecond)
    # breakpoint()
    
    # convert to ints
    int_date = []
    int_time = []
    for i in range(len(date)):
      int_date.append(int(date[i]))
      int_time.append(int(time[i]))
    
    # convert to time_str
    time_str =  [str(x)+str(y) for x, y in zip(int_date, int_time)]
    
    # convert to date-time vector
    datetimes =[]
    for i in range(len(date)):
        datetimes.append(datetime.datetime.strptime(time_str[i][0:13], "%Y%j%H%M%S"))
      
    return datetimes

def date_stringformat(time):
  
    date_nosep =[]
    for i in range(len(time)):
        date_nosep.append(str(time.date[i])[0:4]  + str(time.date[i])[5:7] + str(time.date[i])[8:10])
    return date_nosep

def loop_over_all_days(hdf_files_list, grp, var):
    'Loop pandas output over all days in deployment'
    
    for k in range(len(hdf_files_list)):
        df_tot_k = read_hdf_L2_daily(hdf_files_list[k], grp, var)  # Rrs that has
        if k == 0:
            df_tot = df_tot_k 
        else:
            df_tot = pd.concat([df_tot, df_tot_k])
        
    return df_tot


def loop_over_all_days_anc(hdf_files_list):
    'Loop pandas output for anc data over all days in deployment'
    
    for k in range(len(hdf_files_list)):
        df_tot_k = read_anc_daily(hdf_files_list[k])
        if k == 0:
            df_tot = df_tot_k 
        else:
            df_tot = pd.concat([df_tot, df_tot_k])
        
    return df_tot


def spec_plot(spec, name, units):
       
       wl = np.array([float(i) for i  in (list(spec.keys()))])
       

       plt.figure(figsize =(12,12))
       plt.rcParams.update({'font.size': 18})
       plt.title(name + ' : ' + str(spec.index[0].date()) + ' to '  + str(spec.index[-1].date()))
       plt.plot(wl, np.array(spec).T)
       plt.plot(wl, np.zeros(len(np.array(spec).T)),linestyle='dashed',color='k')
       plt.xlabel('Wavelength [nm]')
       plt.ylabel(units)

       return


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


if __name__ == '__main__':

    
    # directory to ship net CDF files #
    processed_dir = '/mnt/d/AMT31/Optics_all/Data/HSAS/HSAS_data/Processed/'
    jday_dir = glob.glob(processed_dir + '*3*')
    
    ############################################
    # 0. Extract lists of HDF files
    ###############################################
    hdf_files_list = []
    for k in range(8): # up to jday 349
        hdf_files_list.append(glob.glob(jday_dir[k] + '/L2/*hdf*'))

    hdf_files = hdf_files_list[0] # use 0th file to output hdf keys (grps and vars)
    f = h5py.File(hdf_files[0], 'r')

    ################################################################
    # 1. Deriving ancillary data  in daily pandas dataframes
    ################################################################
    
    anc_var = list(f['ANCILLARY'].keys()) # list anc lies
    print(anc_var)
    anc = loop_over_all_days_anc(hdf_files_list)
    
    #  key fields from anc data which go in SeaBASS file in SeaBASS-ready format
    time = anc.index
    date = date_stringformat(time)
    
    lat = np.array(anc['LATITUDE'])
    lon = np.array(anc['LONGITUDE'])
  
    solar_zen = np.array(anc['SZA'])
    solar_az = np.array(anc['SOLAR_AZ'])
    rel_az = np.array(anc['REL_AZ'])
    wind =  np.array(anc['WINDSPEED'])
    
    
    ################################################################
    # 2. Deriving spectral data in daily pandas dataframes
    ################################################################
    

    ref_var = list(f['REFLECTANCE'].keys())
    print(ref_var)
    Rrs = loop_over_all_days(hdf_files_list, 'REFLECTANCE', 'Rrs_HYPER')
    Rrs_M02 = loop_over_all_days(hdf_files_list, 'REFLECTANCE', 'Rrs_HYPER_M02') # this is best choice for Rrs in SeaBASS
    nlw_M02 = loop_over_all_days(hdf_files_list, 'REFLECTANCE', 'nLw_HYPER_M02')
    Rrs_unc = loop_over_all_days(hdf_files_list, 'REFLECTANCE', 'Rrs_HYPER_unc')
    bincount = loop_over_all_days(hdf_files_list, 'REFLECTANCE', 'Ensemble_N') # Rrs/nLW bin count
    bincount = bincount['N'].values.astype(int)
    
    spec_plot(Rrs,'Rrs','[sr$^{-1}$]')
    spec_plot(Rrs_M02,'Rrs_M02','[sr$^{-1}$]')
    spec_plot(Rrs_unc,'Rrs_unc: absolute', '[sr$^{-1}$]')
    


    ###################################
    #  Formatting fields in SeaBASS file
    ###################################

    start_time = str(Rrs_M02.index[0].time()) + '[GMT]' 
    end_time = str(Rrs_M02.index[-1].time()) + '[GMT]' 
    min_lat = str(min(lat))[0:7]
    max_lat = str(max(lat))[0:7]
    min_lon = str(min(lon))[0:7] 
    max_lon = str(max(lon))[0:7]
    time.astype(str)
    
    sb_header = {
    "/begin_header": "",
    "/investigators=": "Tom_Jordan,Gavin_Tilstone,Xavier_Warren",
    "/affiliations=": "Plymouth_Marine_Laboratory,UConn",
    "/contact=": "tjor@pml.ac.uk",
    "/experiment=": "AMT_HSAS",
    "/cruise=": "AMT_31",
    "/station=": "NA",
    "/data_file_name=": "AMT31_Rrs_" + str(date[0]) + "_" + str(date[-1]) + ".sb",
    "/documents=": "TheMagnaCarta_1215edition.doc", 
    "/calibration_files=": "HED2027A.cal,HLD0464A.cal,HLD2054A.cal,HSE2027A.cal,HSL0464A.cal,HSL2054A.cal",
    "/data_type=": "above_water",
    "/water_depth=":"-9999",
    "/measurement_depth=":"0",
    "/data_status=": "preliminary",
    "/start_date=":  str(date[0]),
    "/end_date=": str(date[-1]),
    "/start_time=":  start_time,
    "/end_time=": end_time,
    "/north_latitude=":  max_lat +'[DEG]',
    "/south_latitude=":  min_lat +'[DEG]',
    "/east_longitude=":  max_lon +'[DEG]',
    "/west_longitude=":  min_lon +'[DEG]',
    "/missing=": "-9999",
    "/delimiter=": "comma",
    "/nadir=": "40",
    "/rho_correction=": "Z17",
    "/NIR_residual_correction=": "MA95",
    "/BRDF_correction=": "M02",
    "/fields=":"date,time,lat,lon,RelAz,SZA,wind,bincount,Rrs351.0,Rrs354.3,Rrs357.6,Rrs360.9,Rrs364.2,Rrs367.5,Rrs370.8,Rrs374.1,Rrs377.4,Rrs380.7,Rrs384.0,Rrs387.3,Rrs390.6,Rrs393.9,Rrs397.2,Rrs400.5,Rrs403.8,Rrs407.1,Rrs410.4,Rrs413.7,Rrs417.0,Rrs420.3,Rrs423.6,Rrs426.9,Rrs430.2,Rrs433.5,Rrs436.8,Rrs440.1,Rrs443.4,Rrs446.7,Rrs450.0,Rrs453.3,Rrs456.6,Rrs459.9,Rrs463.2,Rrs466.5,Rrs469.8,Rrs473.1,Rrs476.4,Rrs479.7,Rrs483.0,Rrs486.3,Rrs489.6,Rrs492.9,Rrs496.2,Rrs499.5,Rrs502.8,Rrs506.1,Rrs509.4,Rrs512.7,Rrs516.0,Rrs519.3,Rrs522.6,Rrs525.9,Rrs529.2,Rrs532.5,Rrs535.8,Rrs539.1,Rrs542.4,Rrs545.7,Rrs549.0,Rrs552.3,Rrs555.6,Rrs558.9,Rrs562.2,Rrs565.5,Rrs568.8,Rrs572.1,Rrs575.4,Rrs578.7,Rrs582.0,Rrs585.3,Rrs588.6,Rrs591.9,Rrs595.2,Rrs598.5,Rrs601.8,Rrs605.1,Rrs608.4,Rrs611.7,Rrs615.0,Rrs618.3,Rrs621.6,Rrs624.9,Rrs628.2,Rrs631.5,Rrs634.8,Rrs638.1,Rrs641.4,Rrs644.7,Rrs648.0,Rrs651.3,Rrs654.6,Rrs657.9,Rrs661.2,Rrs664.5,Rrs667.8,Rrs671.1,Rrs674.4,Rrs677.7,Rrs681.0,Rrs684.3,Rrs687.6,Rrs690.9,Rrs694.2,Rrs697.5,Rrs700.8,Rrs704.1,Rrs707.4,Rrs710.7,Rrs714.0,Rrs717.3,Rrs720.6,Rrs723.9,Rrs727.2,Rrs730.5,Rrs733.8,Rrs737.1,Rrs740.4,Rrs743.7,Rrs747.0,Rrs351.0_unc,Rrs354.3_unc,Rrs357.6_unc,Rrs360.9_unc,Rrs364.2_unc,Rrs367.5_unc,Rrs370.8_unc,Rrs374.1_unc,Rrs377.4_unc,Rrs380.7_unc,Rrs384.0_unc,Rrs387.3_unc,Rrs390.6_unc,Rrs393.9_unc,Rrs397.2_unc,Rrs400.5_unc,Rrs403.8_unc,Rrs407.1_unc,Rrs410.4_unc,Rrs413.7_unc,Rrs417.0_unc,Rrs420.3_unc,Rrs423.6_unc,Rrs426.9_unc,Rrs430.2_unc,Rrs433.5_unc,Rrs436.8_unc,Rrs440.1_unc,Rrs443.4_unc,Rrs446.7_unc,Rrs450.0_unc,Rrs453.3_unc,Rrs456.6_unc,Rrs459.9_unc,Rrs463.2_unc,Rrs466.5_unc,Rrs469.8_unc,Rrs473.1_unc,Rrs476.4_unc,Rrs479.7_unc,Rrs483.0_unc,Rrs486.3_unc,Rrs489.6_unc,Rrs492.9_unc,Rrs496.2_unc,Rrs499.5_unc,Rrs502.8_unc,Rrs506.1_unc,Rrs509.4_unc,Rrs512.7_unc,Rrs516.0_unc,Rrs519.3_unc,Rrs522.6_unc,Rrs525.9_unc,Rrs529.2_unc,Rrs532.5_unc,Rrs535.8_unc,Rrs539.1_unc,Rrs542.4_unc,Rrs545.7_unc,Rrs549.0_unc,Rrs552.3_unc,Rrs555.6_unc,Rrs558.9_unc,Rrs562.2_unc,Rrs565.5_unc,Rrs568.8_unc,Rrs572.1_unc,Rrs575.4_unc,Rrs578.7_unc,Rrs582.0_unc,Rrs585.3_unc,Rrs588.6_unc,Rrs591.9_unc,Rrs595.2_unc,Rrs598.5_unc,Rrs601.8_unc,Rrs605.1_unc,Rrs608.4_unc,Rrs611.7_unc,Rrs615.0_unc,Rrs618.3_unc,Rrs621.6_unc,Rrs624.9_unc,Rrs628.2_unc,Rrs631.5_unc,Rrs634.8_unc,Rrs638.1_unc,Rrs641.4_unc,Rrs644.7_unc,Rrs648.0_unc,Rrs651.3_unc,Rrs654.6_unc,Rrs657.9_unc,Rrs661.2_unc,Rrs664.5_unc,Rrs667.8_unc,Rrs671.1_unc,Rrs674.4_unc,Rrs677.7_unc,Rrs681.0_unc,Rrs684.3_unc,Rrs687.6_unc,Rrs690.9_unc,Rrs694.2_unc,Rrs697.5_unc,Rrs700.8_unc,Rrs704.1_unc,Rrs707.4_unc,Rrs710.7_unc,Rrs714.0_unc,Rrs717.3_unc,Rrs720.6_unc,Rrs723.9_unc,Rrs727.2_unc,Rrs730.5_unc,Rrs733.8_unc,Rrs737.1_unc,Rrs740.4_unc,Rrs743.7_unc,Rrs747.0_unc",
    "/units=":"yyyymmdd,hh:mm:ss,degrees,degrees,degrees,degrees,m/s,none,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr,1/sr",
    "/end_header": "",
    }


    # format as dataframe
    sb_df = pd.DataFrame(index=time)
    sb_df['date'] = date
    sb_df['time'] = time.time.astype(str)
    sb_df['lat'] = lat
    sb_df['lon'] = lon
    sb_df['RelAz'] =  rel_az 
    sb_df['SZA'] =  solar_zen
    sb_df['wind'] = wind
    sb_df['bincount'] = bincount
    
    Rrs_M02 = Rrs_M02.iloc[:,0:121]
    Rrs_unc= Rrs_unc.iloc[:,0:121] # 120 wl bins matches default HCP seabass files
    
    sb_df = pd.concat([sb_df, Rrs_M02,Rrs_unc],axis=1,ignore_index=True)
    # sb_df.replace([np.inf, -np.inf], -int(9999) , inplace=True)
    
    filename = '/mnt/d/AMT31/Optics_all/Data/HSAS/HSAS_data/Processed/DailyRrsSb/' + "AMT31_Rrs_" + str(date[0]) + "_" + str(date[-1]) + ".sb"
    export_2_seabass(sb_header, sb_df, filename)
    
    fcheckdir = '/mnt/d/AMT31/Optics_all/Source_HyperCP_HSAS/Fcheck/'
    run_fcheck(filename, fcheckdir)