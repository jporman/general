#!Python3
# celestial_calcs - a module to carry out astronomical time conversions

from cmath import pi
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
    plus_minus = '-' if (fix_degrees == 0 and degrees < 0) else ''
    return plus_minus+str(fix_degrees)+'° '+str(abs(minutes))+ "' " + str(abs(seconds))+'"'

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

def degMS_to_degrees(dec):
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
    if ut_date >= dt.datetime(1582, 10, 15, 0, 0, 0):
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

def obliquity_of_ecliptic(year, month = 1, day = 1, hour = 0, minute = 0, second = 0):
    """Must pass at least the year of the epoch for which ecliptic is
    being calculated. Return decimal for ease of use in other calculations"""
    standard_epoch = dt.datetime(year, month, day, hour, minute, second)
    julian_day = date_to_JD(standard_epoch)
    julian_centuries = (julian_day - 2_451_545.0)/36_525 # Nbr julian centures since 1/0.5/2000
    De = 46.815*julian_centuries + 0.000_6*julian_centuries**2 \
        - 0.001_81*julian_centuries**3
    epsilon_zero = 23.439292 # equal to obliquity of ecliptic at 1/0.0/2000, or 23deg 26' 21.45"
    epsilon = epsilon_zero - De/3600
    return epsilon

def ecliptic_to_equatorial(epoch = 2000, ecliptic_lat = '0deg 0m 0s', ecliptic_long = '0deg 0m 0s'):
    epsilon = obliquity_of_ecliptic(epoch)
    beta_var = degMS_to_degrees(ecliptic_lat)
    lambda_var = degMS_to_degrees(ecliptic_long)

    # Convert to radians
    epsilon_rad = epsilon/360*(2*math.pi)
    beta_var_rad = beta_var/360*(2*math.pi)
    lambda_var_rad = lambda_var/360*(2*math.pi)
    
    # Calc declination
    T = math.sin(beta_var_rad)*math.cos(epsilon_rad) + math.cos(beta_var_rad)*math.sin(epsilon_rad)*math.sin(lambda_var_rad)
    declination = math.asin(T)
    declination = declination/(2*math.pi)*360

    # Calc right ascension
    y = math.sin(lambda_var_rad)*math.cos(epsilon_rad) - math.tan(beta_var_rad)*math.sin(epsilon_rad)
    x = math.cos(lambda_var_rad)
    right_ascension = math.atan(y/x)
    quad_adj = 0
    if y>0 and x>0: 
        quad_adj += 0 
    elif y>0 and x<0: 
        quad_adj += math.pi
    elif y<0 and x>0: 
        quad_adj += 2*math.pi 
    elif y<0 and x<0: 
        quad_adj += math.pi
    right_ascension += quad_adj
    right_ascension = right_ascension/(2*math.pi)*360

    # Convert to HA or degrees_MS format
    declination = degrees_to_degreesMS(declination)
    right_ascension = degrees_to_HA(right_ascension)
    equatorial_coords = {'declination': declination, 'right_ascension': right_ascension}

    return equatorial_coords

def equatorial_to_ecliptic(epoch = 2000, right_ascension = '0h0m0s', declination='0deg0m0s'):
    """Pass RA in hms format and declination in degMS"""
    epsilon = obliquity_of_ecliptic(epoch)
    epsilon_rads = epsilon/360*2*math.pi
    ra = HA_to_degrees(right_ascension)
    ra_rads = ra/360*2*math.pi
    dec = degMS_to_degrees(declination)
    dec_rads = dec/360*2*math.pi

    T_var = math.sin(dec_rads) * math.cos(epsilon_rads) - math.cos(dec_rads) * math.sin(epsilon_rads) * math.sin(ra_rads)
    ecliptic_long = math.asin(T_var)
    ecliptic_long = ecliptic_long/(2*math.pi)*360

    y = math.sin(ra_rads) * math.cos(epsilon_rads) + math.tan(dec_rads)*math.sin(epsilon_rads)
    x = math.cos(ra_rads)
    ecliptic_lat = math.atan(y / x)
    quad_adj = 0
    if y>0 and x>0: 
        quad_adj += 0 
    elif y>0 and x<0: 
        quad_adj += math.pi
    elif y<0 and x>0: 
        quad_adj += 2*math.pi 
    elif y<0 and x<0: 
        quad_adj += math.pi
    ecliptic_lat += quad_adj
    ecliptic_lat = ecliptic_lat/(2*math.pi)*360

    # Convert to degMS format
    ecliptic_long = degrees_to_degreesMS(ecliptic_long)
    ecliptic_lat = degrees_to_degreesMS(ecliptic_lat)

    ecliptic_coords = {'ecliptic longitude': ecliptic_long, 'ecliptic latitude': ecliptic_lat, 'epoch': epoch}

    return ecliptic_coords

def galactic_to_equatorial(galactic_lat = 0, galactic_long = 0, epoch= 1950):
    """Convertes 1950 or J2000 epoch galactic coords to equatorial coordinates"""
    if epoch == 1950:
        gnp_RA_rad = 192.25 / 360 * 2 * math.pi
        gnp_dec_rad = 27.4 / 360 * 2 * math.pi
        N_zero = 33 / 360 * 2 * math.pi # Galactic north pole ascending longitude
    elif epoch == 2000:
        gnp_RA_rad = 192.8598 / 360 * 2 * math.pi
        gnp_dec_rad = 27.128027 / 360 * 2 * math.pi
        N_zero = 32.9319 / 360 * 2 * math.pi # Galactic north pole ascending longitude
    else:
        print('Invalid epoch entered, using 1950 epoch')
        gnp_RA_rad = 192.25 / 360 * 2 * math.pi
        gnp_dec_rad = 27.4 / 360 * 2 * math.pi
        N_zero = 33 / 360 * 2 * math.pi # Galactic north pole ascending longitude
    
    
    galactic_lat = degMS_to_degrees(galactic_lat)
    galactic_long = degMS_to_degrees(galactic_long)
    galactic_lat_rad = galactic_lat / 360 * 2 * math.pi
    galactic_long_rad = galactic_long / 360 * 2 * math.pi

    T = math.cos(galactic_lat_rad)*math.cos(gnp_dec_rad)*math.sin(galactic_long_rad - N_zero) \
        + math.sin(galactic_lat_rad)*math.sin(gnp_dec_rad)
    equatorial_dec_rad = math.asin(T)

    y = math.cos(galactic_lat_rad)*math.cos(galactic_long_rad - N_zero)
    x = math.sin(galactic_lat_rad)*math.cos(gnp_dec_rad) - \
        math.cos(galactic_lat_rad)*math.sin(gnp_dec_rad)*math.sin(galactic_long_rad - N_zero)
    equatorial_RA_rad = math.atan(y/x)
    quad_adj = 0
    if y>0 and x>0: 
        quad_adj += 0 
    elif y>0 and x<0: 
        quad_adj += math.pi
    elif y<0 and x>0: 
        quad_adj += 2*math.pi 
    elif y<0 and x<0: 
        quad_adj += math.pi
    equatorial_RA_rad += (quad_adj + gnp_RA_rad)

    if equatorial_RA_rad > 2 * math.pi:
        equatorial_RA_rad -= 2 * math.pi
    
    equatorial_RA = equatorial_RA_rad * 360 / (2 * math.pi)
    equatorial_RA = degrees_to_HA(equatorial_RA)
    equatorial_dec = equatorial_dec_rad * 360 / (2 * math.pi)
    equatorial_dec = degrees_to_degreesMS(equatorial_dec)
    equatorial_coords = {'right_ascension': equatorial_RA, 'declination': equatorial_dec}
    
    return equatorial_coords

def equatorial_to_galactic(equatorial_RA = 0, equatorial_dec = 0, epoch = 1950):

    equatorial_RA = HA_to_degrees(equatorial_RA)
    equatorial_dec = degMS_to_degrees(equatorial_dec)
    eq_RA_rad = equatorial_RA / 360 * 2 * math.pi
    eq_dec_rad = equatorial_dec / 360 * 2 * math.pi

    if epoch == 1950:
        gnp_RA_rad = 192.25 / 360 * 2 * math.pi
        gnp_dec_rad = 27.4 / 360 * 2 * math.pi
        N_zero = 33 / 360 * 2 * math.pi # Galactic north pole ascending longitude
    elif epoch == 2000:
        gnp_RA_rad = 192.8598 / 360 * 2 * math.pi
        gnp_dec_rad = 27.128027 / 360 * 2 * math.pi
        N_zero = 32.9319 / 360 * 2 * math.pi # Galactic north pole ascending longitude
    else:
        print('Invalid epoch entered, using 1950 epoch')
        gnp_RA_rad = 192.25 / 360 * 2 * math.pi
        gnp_dec_rad = 27.4 / 360 * 2 * math.pi
        N_zero = 33 / 360 * 2 * math.pi # Galactic north pole ascending longitude

    T_zero = math.cos(eq_dec_rad)*math.cos(gnp_dec_rad)*math.cos(eq_RA_rad - gnp_RA_rad) +\
        math.sin(eq_dec_rad)*math.sin(gnp_dec_rad) 
    galactic_lat_rad = math.asin(T_zero)

    y = math.sin(eq_dec_rad)-math.sin(galactic_lat_rad)*math.sin(gnp_dec_rad)
    x = math.cos(eq_dec_rad)*math.sin(eq_RA_rad - gnp_RA_rad)*math.cos(gnp_dec_rad)
    galactic_long_rad = math.atan(y/x)
    quad_adj = 0
    if y>0 and x>0: 
        quad_adj += 0 
    elif y>0 and x<0: 
        quad_adj += math.pi
    elif y<0 and x>0: 
        quad_adj += 2*math.pi 
    elif y<0 and x<0: 
        quad_adj += math.pi
    galactic_long_rad += (quad_adj + N_zero)

    if galactic_long_rad > 2 * math.pi:
        galactic_long_rad -= 2 * math.pi
    
    galactic_long = galactic_long_rad * 360 / (2 * math.pi)
    galactic_long = degrees_to_degreesMS(galactic_long)
    galactic_lat = galactic_lat_rad * 360 / (2 * math.pi)
    galactic_lat = degrees_to_degreesMS(galactic_lat)
    galactic_coords = {'galactic_longitude': galactic_long, 'galactic_latitude': galactic_lat}
    
    return galactic_coords

def precession_corrections(RA_uncorrected = None, dec_uncorrected = None, epoch_from = 1950, epoch_to = 1950):
    """returns correction to RA in hours and dec in degrees, as well as corrected RA and Dec"""
    RA_uncorrected = HA_to_degrees(RA_uncorrected)
    dec_uncorrected = degMS_to_degrees(dec_uncorrected)
    T_var = (epoch_to - 1900) / 100
    M_var = 3.07234 + .00186*T_var # Seconds of time
    N_d = 20.0468 - 0.0085*T_var # Arcseconds
    N_t = N_d / 15 # Converts arcseconds to seconds of time
    D_var = epoch_to - epoch_from

    # delta_RA = correction to RA_uncorrected in degrees
    delta_RA = (M_var + N_t*math.sin(RA_uncorrected / 360 * 2 * math.pi) * \
        math.tan(dec_uncorrected / 360 * 2 * math.pi)) * D_var
    delta_RA = delta_RA / 3600 # Conversion to hours
    RA_corrected = RA_uncorrected + delta_RA * 15 # Conversion to degrees...do it elsewhere
    
    # delta_dec = correction to dec_uncorrected in degrees
    delta_dec = (N_d*math.cos(RA_uncorrected / 360 * 2 * math.pi)) * D_var # In arcseconds
    delta_dec = delta_dec / 3600
    dec_corrected = dec_uncorrected + delta_dec

    RA_corrected = degrees_to_HA(RA_corrected)
    dec_corrected = degrees_to_degreesMS(dec_corrected)
    corrections_degrees = {'delta_RA': delta_RA, 'delta_dec': delta_dec, \
        'RA_corrected': RA_corrected, 'dec_corrected': dec_corrected}

    return corrections_degrees

def solve_keppler(orbital_eccentricity = None, mean_anomaly = None, \
    eccentric_anomaly = None, solve_method = 'Simple', termination_criteria = .0001):
    """Return mean anomaly if eccentric is given, or eccentric if mean is given. One can select numerical method"""
    while not orbital_eccentricity:
        orbital_eccentricity = input('Please enter an orbital eccentricity: ')

    if eccentric_anomaly:
        E_rads = eccentric_anomaly /360 * 2 * math.pi
        mean_anomaly = e_rads - orbital_eccentricity * math.sin(e_rads)
        mean_anomaly = mean_anomaly / (2*math.pi) * 360
        return mean_anomaly
    elif mean_anomaly:
        M_rads = mean_anomaly / 360 * 2 * math.pi
        if solve_method == 'Simple':
            E_rads_prior = M_rads
            E_rads_current = M_rads + orbital_eccentricity * math.sin(E_rads_prior)
            iteration_count = 1
            while abs(E_rads_current - E_rads_prior) > termination_criteria:
                E_rads_prior = E_rads_current
                E_rads_current = M_rads + orbital_eccentricity * math.sin(E_rads_prior)
                iteration_count += 1
        elif solve_method == 'Newton/Raphson':
            E_rads_prior = M_rads
            E_rads_current = E_rads_prior - \
                (E_rads_prior - orbital_eccentricity * math.sin(E_rads_prior) - M_rads) / \
                (1 - orbital_eccentricity * math.cos(E_rads_prior))
            iteration_count = 1
            while abs(E_rads_current - E_rads_prior) > termination_criteria:
                E_rads_prior = E_rads_current
                E_rads_current = E_rads_prior - \
                    (E_rads_prior - orbital_eccentricity * math.sin(E_rads_prior) - M_rads) / \
                    (1 - orbital_eccentricity * math.cos(E_rads_prior))
                iteration_count += 1
        else:
            print('Numerical method not recognized')
            return None
        print(f'Keppler solution found in {iteration_count} iterations using {solve_method} numerical method')
        eccentric_anomaly = E_rads_current * 360 / (2 * math.pi)
        return eccentric_anomaly
    else:
        print('Missing anomaly input. Cannot compute')
        return None

    
###########################################
#---------------MAIN SCRIPT---------------#
###########################################

# ecliptic_coords = equatorial_to_ecliptic(epoch = 2000, \
#     right_ascension = '11h 10m 13s', declination = '30deg 05m 40s')

# print('Ecliptic longitude: ', ecliptic_coords['ecliptic longitude'], '\n'\
#     'Ecliptic latitude: ', ecliptic_coords['ecliptic latitude'])

# eq_coords = ecliptic_to_equatorial(epoch = 2000, \
#     ecliptic_long = '120deg 30m 30s', ecliptic_lat = '0deg 00m 00s')

# print('Right ascension: ', eq_coords['right_ascension'], '\n'\
#     'Declination: ', eq_coords['declination'])

# eq_coords = galactic_to_equatorial(galactic_lat = '30deg 25m 40s', \
#     galactic_long = '120deg 00m 00s', epoch = 2000)
# print('Right Ascension: ', eq_coords['right_ascension'], \
#     '\nDeclination: ', eq_coords['declination'])

# galactic_coords = equatorial_to_galactic(equatorial_RA = '11h 10m 13s', \
#     equatorial_dec = '30deg 05m 40s', epoch = 2000)
# print('Galactic latitude: ', galactic_coords['galactic_latitude'], \
#     '\nGalactic longitude: ', galactic_coords['galactic_longitude'])

# corrections = precession_corrections(RA_uncorrected = '12h 34m 34s', \
#     dec_uncorrected = '29deg 49m 08s', epoch_from = 2000, epoch_to = 2015.0)
# print('RA correction: ', corrections['delta_RA'], \
#     '\nDec correction: ', corrections['delta_dec'], \
#     '\nRA corrected: ', corrections['RA_corrected'], \
#     '\nDec corrected: ', corrections['dec_corrected'])

print(solve_keppler(orbital_eccentricity = .850000, \
    mean_anomaly = 5.498078, \
    solve_method = 'Newton/Raphson', \
    termination_criteria = .000_002))


