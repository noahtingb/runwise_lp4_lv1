import numpy as np
import coden.Solweig as so
import coden.Clearnessindex as ci
import coden.Pet_Calc as p
import coden.Sun_Position as sp
import datetime

def Solweig_2015a_metdata_noload(year,doy,hour,minu, location, UTC=0):
    """
    This function is used to process the input meteorological file.
    It also calculates Sun position based on the time specified in the met-file

    :param inputdata:
    :param location:
    :param UTC:
    :return:
    """
    dectime = doy+hour / 24 + minu / (60*24.)
    halftimestepdec = 0
    sunmaximum = 0.
    
    sunmax = dict()

    YMD = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1)
    # Finding maximum altitude in 15 min intervals (20141027)
    fifteen = 0.
    sunmaximum = -90.
    sunmax={'zenith': 90.}
    while sunmaximum <= 90. - sunmax['zenith']:
                sunmaximum = 90. - sunmax['zenith']
                fifteen = fifteen + 15. / 1440.
                YMDHM = YMD+datetime.timedelta(days=(60*10)/1440.0 + fifteen)
                time = {"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC}
                sunmax = sp.sun_position(time,location)
    altmax = sunmaximum

    YMDHM = YMD + datetime.timedelta(hours=hour) + datetime.timedelta(minutes=minu) - datetime.timedelta(days=halftimestepdec)

    time = {"year":YMDHM.year,'month': YMDHM.month,'day':YMDHM.day,'hour': YMDHM.hour,'min':YMDHM.minute,'sec':0,"UTC":UTC}
    sun = sp.sun_position(time, location)
    altitude = 90. - sun['zenith']
    azimuth = sun['azimuth']
    zen = sun['zenith'] * (3.141592653589793238/180.)
        # day of year and check for leap year
    if (time['year'] % 4) == 0 and ( ((time['year'] % 100)==0 and (time['year'] % 400) == 0) or ((time['year'] % 100)!=0)):
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    jday = sum(dayspermonth[0:time['month'] - 1]) + time['day']

    return None, altitude, azimuth, zen, jday, None, dectime, altmax
def indexflask(form):
        month = int(form["month"])
        day = int(form["day"])
        hour = int(form["hour"])
        year = int(form["year"])
        minu = 30
        Ta = float(form["Ta"])        
        RH = float(form["RH"])
        Ws = float(form["Ws"])
        location = form["loc"]

        if month > 12 or month < 0:
            print("petresult.html","Incorrect month filled in")
        if day > 31 or day < 0:
            print("petresult.html","Incorrect day filled in")
        if hour > 23 or hour < 0:
            print("petresult.html","Incorrect hour filled in")
        if Ta > 60 or Ta < -75:
            print("petresult.html", "Unreasonable air temperature filled in",Ta)
        if RH > 100 or RH < 0:
            print("petresult.html", "Unreasonable relative humidity filled in")
        if Ws > 100 or Ws < 0:
            print("petresult.html", "Unreasonable Wind speed filled in")
        #day of year
        if (year % 4) == 0 and ( ((year % 100)==0 and (year % 400) == 0) or ((year % 100)!=0)):
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        doy= sum(dayspermonth[0:month - 1]) + day
        # Currently looked to England
        UTC = 0
        # Radiation
        P = -999. 
        _, altitude, _, zen, jday, _, _, _ = Solweig_2015a_metdata_noload(year,doy,hour,minu, location, UTC)
        #print(altitude,zen,jday)
        if altitude > 0.:
            I0, _, _= ci.clearnessindex_2013b(zen, jday, Ta, RH / 100., 40, location, P)
            radG = I0[0]*(1-form["c1"]*1e-2)*(1-form["c2"]*1e-2)*(1-form["c3"]*1e-2)
        else:
            radG = 0.
        # Main calculation
        if Ta is not None and RH is not None and Ws is not None and radG is not None:
            _, resultPET, _ = petcalc(Ta, RH, Ws, radG, year, month, day, hour, minu,location)
        return resultPET

def petcalc(Ta, RH, Ws, radG, year, month, day, hour, minu,location):
#    sh = 1.  # 0 if shadowed by building
#    vegsh = 1.  # 0 if shadowed by tree

    # Location and time settings. Should be moved out later on
    UTC = 0

    # Human parameter data. Should maybe be move out later on
    absK = 0.70
    absL = 0.95
    pos = 0
    mbody = 75.
    ht = 180 / 100.
    clo = 0.9
    age = 35
    activity = 80.
    sex = 1

    if pos == 0:
        Fside = 0.22
        Fup = 0.06
        Fcyl = 0.28
    else:
        Fside = 0.166666
        Fup = 0.166666
        Fcyl = 0.2


    # Environmental data. Should maybe bo moved out later on.
    albedo_b = 0.2
    albedo_g = 0.15
    ewall = 0.9
    eground = 0.95
#    svf = 0.6

    # Meteorological data, Should maybe be move out later on.
    onlyglobal = 1

    doy = day_of_year(year, month, day)

    radD = -999.
    radI = -999.
    P   = -999.
    _, altitude, azimuth, zen, jday, _, dectime, altmax = Solweig_2015a_metdata_noload(year,doy,hour, minu, location, UTC)

    svfalfa = np.arcsin(np.exp((np.log((1.-0.6))/2.)))

    # %Creating vectors from meteorological input
    
    TgK = 0.37
    Tstart = -3.41
    TmaxLST = 15
    TgK_wall = 0.58
    Tstart_wall = -3.41
    TmaxLST_wall = 15

    # If metfile starts at night
    CI = 1.

    skyvaultalt = np.atleast_2d([])
    skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
    skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
    for j in range(7):
        for k in range(1, int(360/skyvaultaziint[j]) + 1):
            skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])

    skyvaultalt = np.append(skyvaultalt, 90)

    diffsh = np.zeros((145))
    svfalfadeg = svfalfa / (np.pi / 180.)
    for k in range(0, 145):
        if skyvaultalt[k] > svfalfadeg:
            diffsh[k] = 1

    # main loop
    # Nocturnal cloudfraction from Offerle et al. 2003
    #cant have this funktion with this data

    Tmrt,CI = so.Solweig1D_2020a_calc( albedo_b, absK, absL, ewall,
                                                            Fside, Fup, Fcyl,
                                                            altitude, azimuth, zen, jday,
                                                            onlyglobal, location, dectime, altmax,
                                                            Ta, RH, radG, radD, radI, P, TgK, Tstart, albedo_g, eground, TgK_wall, Tstart_wall,
                                                            TmaxLST, TmaxLST_wall, svfalfa, CI, 1, diffsh)

    # Recalculating wind speed based on pwerlaw
    WsPET = (1.1 / 10) ** 0.2 * Ws
    resultPET = p._PET(Ta, RH, Tmrt, WsPET, mbody, age, ht, activity, clo, sex)

    return Tmrt, resultPET, None

def day_of_year(yyyy, month, day):
    if (yyyy % 4) == 0 and ( ((yyyy % 100)==0 and (yyyy % 400) == 0) or ((yyyy % 100)!=0)):
        dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return sum(dayspermonth[0:month - 1]) + day
