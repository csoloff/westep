import matplotlib.pyplot as plt

import csv

from astropy.io.votable import parse

from astropy.io.votable import parse_single_table

from astropy.io import ascii

import astropy.coordinates as coord

from astropy.coordinates import SkyCoord

import astropy.units as u

from astropy.io import ascii

import asciitable

import numpy as np

from prettytable import PrettyTable

import math

from datetime import datetime, timedelta

from suntime import Sun, SunTimeException

from astroquery.gaia import Gaia

import warnings

from tqdm import tqdm

import time

import pandas as pd

import os

warnings.filterwarnings("ignore")

def read(file, column):

    list=[]

    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                list.append(float(row[column]))
            except:
                pass
    return(list)


votable = parse_single_table("Master.xml")

results = votable.array


rot = 23.934472 # sidereal day length (hr)


today = datetime.now()

d = (today - datetime.fromisoformat('2021-09-23')).days

# Van Vleck coordinates
long = -72.6595
lat = 41.5556

sun = Sun(lat, long)

today_sr = sun.get_sunrise_time() - timedelta(hours=4)
today_ss = sun.get_sunset_time() - timedelta(hours=4)

night = today_ss - today_sr

UT7 = 2

lmst7 = (24.06570982)*d + (today_ss.replace(tzinfo=None) - datetime.combine(today, datetime.min.time()) ).seconds/86400*23.93446959



while lmst7 > 24:
    lmst7-=24



def UT(t): # given hours since 7:00 PM, ouputs UT in hrs
    universal = 2+t
    if universal >= 24:
        universal -= 24
    return universal

def Local(t): # given hours since 7:00 PM, ouputs local time in hrs
    l = 19+t
    if l >= 24:
        l -= 24
    return l

#print(Local(1))

def HA(t, RA, Dec): # given hours since 7:00 PM, ouputs hour angle in hrs
    HourAngle = lmst7+(t/24)*rot-RA/360*24
    while HourAngle<-12 or HourAngle>12:
        if HourAngle<-12:
            HourAngle+=24
        if HourAngle>12:
            HourAngle-=24
    return(HourAngle)

def Alt(t, RA, Dec): # given hours since 7:00 PM, ouputs altitude in degrees
    Altitude = (360/(2*math.pi))*math.asin(math.sin(2*math.pi * lat/360)*math.sin(2*math.pi * Dec/360)+math.cos(2*math.pi * lat/360)*math.cos(2*math.pi * Dec/360)*math.cos(2*math.pi * HA(t, RA, Dec)/24))
    return(Altitude)

def Az(t, RA, Dec):  # given hours since 7:00 PM, ouputs azimuth in degrees
    Azimuth = (360/(2*math.pi))*math.acos((math.sin(2*math.pi * Dec/360)-math.sin(2*math.pi * Alt(t, RA, Dec)/360)*math.sin(2*math.pi * lat/360))/(math.cos(2*math.pi * Alt(t, RA, Dec)/360)*math.cos(2*math.pi * lat/360)))
    return(Azimuth)

def X(t, RA, Dec): # given hours since 7:00 PM, outputs airmass
    z = 90 - Alt(t, RA, Dec)
    airmass = 1/(math.cos(2*math.pi * z/360))*(1-0.0012*(1/((math.cos(2*math.pi * z/360))**2)-1))
    return(airmass)

times=[]
hrs = round(night.seconds/(60*60))
for i in range(0, hrs+1):
    times.append(i)




def tlist(Gaia_run):
    acc=0
    acc2=-1
    name=[]
    sp=[]
    dec=[]
    ra=[]
    teff=[]
    vmag=[]
    distance=[]
    diameter=[]
    rot=[]
    bs=[]
    bs2=[]
    bs3=[]
    rise=[]

    ly=[]
    percent=0
    leng=len(results['MAIN_ID'])
    for x in tqdm(results['MAIN_ID']):
        width = u.Quantity(0.3983, u.deg)
        height = u.Quantity(0.3983, u.deg)
        d = float(results['DEC_d'][acc])
        r = float(results['RA_d'][acc])

        y=0
        b=True
        for z in range(0, 10*hrs+1):
            if X(z/10, r, d) < 1.5 and X(z/10, r, d) > 0:
                if b:
                    rise_i = today_ss + timedelta(hours=z/10)
                    b=False
                y+=1

        spectral = str(results['SP_TYPE'][acc])#[2:-1]
        if y>50 and float(results['Distance_distance'][acc]) >0 and 'V' in spectral and 'IV' not in spectral and results['Fe_H_Teff'][acc] > 6500 and results['Fe_H_Teff'][acc] < 10000: # and float(results['ROT_Vsini'][acc]) >0; d > 30 and r < RAMAX and r > RAMIN
            if Gaia_run:
                c = SkyCoord(results['RA_d'][acc], results['DEC_d'][acc], frame = 'icrs', unit = 'deg')

                result = Gaia.query_object(coordinate = c, width=width, height=height) # queries using Gaia

                bs.append(len(np.argwhere(result['phot_g_mean_mag'] < 15))-1)
                bs2.append(len(np.argwhere(result['phot_g_mean_mag'] < 12))-1)
                bs3.append(len(np.argwhere(result['phot_g_mean_mag'] < 9))-1)
                acc2+=1
            else:
                bs.append('NA')
                bs2.append('NA')
                bs3.append('NA')
                acc2+=1
            w=True
            if bs[acc2]==1:
                w = bs[acc2]==1 and bs2[acc2]>0 and bs3[acc2]>0
            if bs[acc2]!=0 and w:
                name.append(str(x))#[2:-1])
                sp.append(spectral)
                dec.append(results['DEC_d'][acc])
                ra.append(results['RA_d'][acc])
                teff.append(results['Fe_H_Teff'][acc])
                vmag.append(results['FLUX_V'][acc])
                distance.append(results['Distance_distance'][acc])
                diameter.append(results['Diameter_diameter'][acc])
                rot.append(results['ROT_Vsini'][acc])
                ly.append(y/10)
                rise.append(str(rise_i)[0:16])
            else:
                del bs[acc2]
                del bs2[acc2]
                del bs3[acc2]
                acc2-=1
        acc+=1


    master = [name, sp, dec, ra, teff, vmag, distance, diameter, rot]
    tab = PrettyTable()
    tab.field_names = ['Name', 'sp', 'ra', 'dec', 'teff', 'vmag', 'distance', 'diameter', 'rot', 'Hrs of X < 1.5', 'N<15', 'N<12', 'N<9', 'Start']
    acc=0
    for x in master[0]:
        tab.add_row([str(x), str(master[1][acc]), master[3][acc], master[2][acc], master[4][acc], master[5][acc], master[6][acc], master[7][acc], master[8][acc], ly[acc], bs[acc], bs2[acc], bs3[acc], rise[acc]])
        acc+=1
    tab.sortby = 'Start'#'rot'#'distance'
    print(tab.get_string())
    print(len(master[0]))
    master.extend([ly, bs, bs2, bs3, rise])
    return(master)


output = tlist(False)  # change boolean to toggle reference star search (reduces run time)
