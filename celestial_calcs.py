#!Python3
# celestial_calcs - a module to carry out astronomical time conversions

import datetime as dt
import dateutil as du
import numpy as np
import math
import re

def degrees_to_degreesMS(degrees):
    """need to round seconds up/down rather than fix"""
    fix_degrees = int(np.fix(degrees))
    minutes = int(np.fix((degrees - fix_degrees) * 60))
    seconds = round((degrees - fix_degrees - minutes / 60) * 3600)
    return str(fix_degrees)+'° '+str(minutes)+ "' " + str(seconds)+'"'

def degrees_to_HA(degrees):
    hour = int(np.fix(degrees/15))
    minutes = int(np.fix((degrees - hour * 15) * 4))
    seconds = round((degrees - hour * 15 - minutes * .25) * 3600/15)
    return str(hour)+'h '+str(minutes)+'\''+str(seconds)+'\"'

def HA_to_degrees(HA):
    hms_regex = re.compile(r'(\d\d?)(h?)(\s?)(\d\d?)(m?|\'?)(\s?)(\d\d?\.?\d*)(s?|\"?)(\s?)')
    hms = hms_regex.search(HA)
    hours = hms[1]
    minutes = hms[4]
    seconds = hms[7]
    degrees = float(hours) * 15 + float(minutes)/4 + float(seconds)/240
    return degrees

def dec_to_degrees(dec):
    degMS_regex = re.compile(r'(-?)(\d\d?\d?)(deg|°|\*)(\s?)(\d\d?)(m?|\'?)(\s?)(\d\d?\.?\d*)(s?|\"/)(\s?)')
    degMS = degMS_regex.search(dec)
    degrees = degMS[2]
    minutes = degMS[5]
    seconds = degMS[8]
    degrees = float(degrees) + float(minutes)/60 + float(seconds)/3600
    degrees = degrees*-1 if degMS[1] == '-' else degrees
    return degrees

def hms_to_decimal(datetime_object):
    hours = int(datetime_object.strftime('%H'))
    minutes = int(datetime_object.strftime('%M'))
    seconds = float(datetime_object.strftime('%S.%f'))
    decimal_time = hours + minutes/60 + seconds/3600
    return decimal_time

def decimal_to_hms(y, m, d, decimal_time):
    H = int(np.fix(decimal_time))
    M = int(np.fix((decimal_time % 1) * 60))
    S = int(np.fix((decimal_time * 3600) - (H * 3600) - (M * 60)))
    MS = int(np.fix(decimal_time * 3600 * 1000000 - H * 3600 * 1000000 - M * 60 * 1000000 - S * 1000000))
    hms_time = dt.datetime(y, m, d, H, M, S, MS)
    return hms_time

def parse_longitude(longitude_str):   
    regex_longitude_nbr = re.compile(r'\d+') #Define regex to parse longitude degrees
    longitude_nbr = float(regex_longitude_nbr.search(longitude_str).group(0)) #Parse longitude degrees
    regex_east_west = re.compile(r'W|w|E|e|West|west|East|east') #Define regex to parse east-west
    east_west = regex_east_west.search(longitude_str).group(0).lower() #Parse east-west
    return{'longitude':longitude_nbr, 'east_west':east_west}

def date_to_JD(ut_date):
    """convert UT date to Julian date. Date passed shoudl be datetime object"""
    year = int(ut_date.strftime('%Y'))
    month = int(ut_date.strftime('%m'))
    day = int(ut_date.strftime('%d'))
    
    #Step 1
    if month > 2:
        y = year
        m = month
    else:
        y = year - 1
        m = month + 12
    
    #Step 2
    if year < 0:
        T = 0.75
    else:
        T = 0
    
    #Step 3/4
    if ut_date.date() >= dt.date(1582, 10, 15):
        A = np.fix(y/100)
        B = 2 - A + np.fix(A/4)
    else:
        A = 0
        B = 0
    
    #Step 4
    JD = B + np.fix(365.25 * y - T) + np.fix(30.6001 * (m + 1)) + day + 1_720_994.5
    return JD
        

def LCT_to_UT(lct, longitude_str, dst):
    """function to convert LCT to UT"""
    
    #Convert lct to datetime object
    try:
        lct = dt.datetime.strptime(lct, '%m/%d/%Y %H:%M:%S')
    except ValueError:
        try:
            lct = dt.datetime.strptime(lct, '%m/%d/%Y')
        except ValueError:  
            print('Invalid date-time entered')
    
    longitude = parse_longitude(longitude_str)['longitude']
    east_west = parse_longitude(longitude_str)['east_west']
    
    #Determine whether to apply positive or negative adj to LCT
    adj_factor = 1 if east_west == 'e' else -1
    
    #Calculate adjustment to LCT to get to UT (accounting for position in timezone)
    adj_hours = longitude / 15 * adj_factor
    dst_adj = 1 if dst == 'y' else 0
    ut = lct - dt.timedelta(hours = adj_hours + dst_adj)
    
    return ut

def UT_to_GST(ut_date):
    """function to convert UT to GST"""
    #Step 1 - validated
    JD = date_to_JD(ut_date)
    
    #Step 2 - validated
    JD0 = date_to_JD(dt.datetime(int(ut_date.strftime('%Y')), 1, 1, 0, 0, 0)) - 1
    
    #Step 3 - validated
    Days = JD - JD0
    
    #Step 4 - validated
    T = (JD0 - 2_415_020.0) / 36_525.0
    
    #Step 5 - validated
    R = 6.6460656 + (2400.051262 * T) + (.00002581 * (T**2))
    
    #Step 6 - validated
    B = 24 - R + 24 * (int(ut_date.strftime('%Y')) - 1900)
    
    #Step 7 - validated
    T0 = (0.0657098 * Days) - B
    
    #Step 8 - validated
    UT = hms_to_decimal(ut_date)
    
    #Step 9 - validated
    GST = T0 + 1.002738 * UT
    
    #Step 10  validated
    if GST < 0:
        GST += 24
    elif GST > 24:
        GST -= 24
        
    #Convert to datetime and return
    year = int(ut_date.strftime('%Y'))
    month = int(ut_date.strftime('%m'))
    day = int(ut_date.strftime('%d'))
    GST = decimal_to_hms(year, month, day, GST)
    return GST
    
def GST_to_LST(GST_datetime_object, longitude_str):
    #Convert longitude to +/- depending on east vs. west vs GST
    longitude = parse_longitude(longitude_str)['longitude']
    east_west = parse_longitude(longitude_str)['east_west']
    longitude = longitude*-1 if east_west == 'w' else longitude
    
    #Step 1 - Convert GST to decimal format
    GST_decimal = hms_to_decimal(GST_datetime_object)
    
    #Step 2 - Time zone adjustment
    adjustment = longitude / 15
    
    #Step 3 = Calc LST decimal
    LST_decimal = GST_decimal + adjustment
    
    #Step 4 - if LST is negative or if >24 subtract 24h
    if LST_decimal < 0:
        LST_decimal += 24
    elif LST_decimal > 24:
        LST_decimal -= 24
    
    year = int(GST_datetime_object.strftime('%Y'))
    month = int(GST_datetime_object.strftime('%m'))
    day = int(GST_datetime_object.strftime('%d'))
    LST_hms = decimal_to_hms(year, month, day, LST_decimal)
    
    return LST_hms

def altaz_to_eq(alt, az, lat, north_south):
    """function that converts alt-azimuth coordinates to equatorial. Assume args in degrees"""
    
    #Convert degrees to radians
    altitude = alt*np.pi/180
    azimuth = az*np.pi/180
    latitude = lat*np.pi/180
    latitude = latitude if north_south.lower() == 'n' else latitude*-1
    
    #Calculate declination
    T0 = math.sin(altitude)*math.sin(latitude)+math.cos(altitude)*math.cos(latitude)*math.cos(azimuth)
    declination = math.asin(T0)
    
    #Calculate hour angle
    T1 = math.sin(altitude) - math.sin(latitude) * math.sin(declination)
    cosH = T1 / (math.cos(latitude)*math.cos(declination))
    H = math.acos(cosH)*180/np.pi
    hour_angle = 360-H if math.sin(azimuth) > 0 else H
    
    #Convert back to degreesMS format
    declination = declination*180/np.pi
    hour_angle = degrees_to_HA(hour_angle)
    declination = degrees_to_degreesMS(declination)
    
    return {'declination':declination, 'hour_angle':hour_angle}

def eq_to_altaz(HA_deg, dec_deg, lat, north_south):
    """function that converts equatrial coordinates to alt-az. Assume args in degrees"""
    
    #Convert degrees to radians
    HA = HA_deg*np.pi/180
    dec = dec_deg*np.pi/180
    latitude = lat*np.pi/180
    latitude = latitude if north_south.lower() == 'n' else latitude*-1
    
    T0 = math.sin(dec)*math.sin(latitude) + math.cos(dec)*math.cos(latitude)*math.cos(HA)
    h = math.asin(T0)  
    T1 = math.sin(dec) - math.sin(latitude) * math.sin(h)
    T2 = T1 / (math.cos(latitude) * math.cos(h))
    A = math.acos(T2)
    
    h = h*180/np.pi
    A = A*180/np.pi
    A = 360-A if math.sin(HA) > 0 else A
    
    h = degrees_to_degreesMS(h)
    A = degrees_to_degreesMS(A)
    
    altaz_coords = {'altitude': h, 'azimuth': A}
    return altaz_coords



# 1. convert ha and dec to degrees
ha = "7h 00m 00s"
dec = "49deg 54m 20s"
ha = HA_to_degrees(ha)
dec = dec_to_degrees(dec)

# 2. convert to altaz
eq_coords = eq_to_altaz(ha, dec, 80, 'S')
alt = eq_coords['altitude']
az = eq_coords['azimuth']

print(f'alt: {alt}')
print(f'azimuth: {az}')




