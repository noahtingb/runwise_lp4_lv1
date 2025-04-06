import numpy as np
from coden.Clearnessindex import clearnessindex_2013b
import numpy as np

def Perez_v3(zen, azimuth, radD, radI, jday):
    """
    This function calculates distribution of luminance on the skyvault based on
    Perez luminince distribution model.
    
    Created by:
    Fredrik Lindberg 20120527, fredrikl@gvc.gu.se
    Gothenburg University, Sweden
    Urban Climte Group
    
    Input parameters:
     - zen:     Zenith angle of the Sun (in degrees)
     - azimuth: Azimuth angle of the Sun (in degrees)
     - radD:    Horizontal diffuse radiation (W m-2)
     - radI:    Direct radiation perpendicular to the Sun beam (W m-2)
     - jday:    Day of year
    
    Output parameters:
     - lv:   Relative luminance map (same dimensions as theta. gamma)

    :param zen:
    :param azimuth:
    :param radD:
    :param radI:
    :param jday:
    :return:
    """
    acoeff=np.array([[1.3525,-0.2576,-0.269,-1.4366],[-1.2219,-0.773,1.4148,1.1016],[-1.1,-0.2515,0.8952,0.0156],[-0.5484,-0.6654,-0.2672,0.7117],[-0.6,-0.3566,-2.5,2.325],[-1.0156,-0.367,1.0078,1.4051],[-1.,0.0211,0.5025,-0.5119],[-1.05,0.0289,0.426,0.359]])    
    bcoeff=np.array([[-7.6700e-01,7.0000e-04,1.2734e+00,-1.2330e-01],[-2.0540e-01,3.6700e-02,-3.9128e+00,9.1560e-01],[2.7820e-01,-1.8120e-01,-4.5000e+00,1.1766e+00],[7.2340e-01,-6.2190e-01,-5.6812e+00,2.6297e+00],[2.9370e-01,4.9600e-02,-5.6812e+00,1.8415e+00],[2.8750e-01,-5.3280e-01,-3.8500e+00,3.3750e+00],[-3.0000e-01,1.9220e-01,7.0230e-01,-1.6317e+00],[-3.2500e-01,1.1560e-01,7.7810e-01,2.5000e-03]])
    ccoeff=np.array([[2.8,0.6004,1.2375,1.],[6.975,0.1774,6.4477,-0.1239],[24.7219,-13.0812,-37.7,34.8438],[33.3389,-18.3,-62.25,52.0781],[21.,-4.7656,-21.5906,7.2492],[14.,-0.9999,-7.1406,7.5469],[19.,-5.,1.2438,-1.9094],[31.0625,-14.5,-46.1148,55.375]])    
    dcoeff=np.array([[1.8734e+00,6.2970e-01,9.7380e-01,2.8090e-01],[-1.5798e+00,-5.0810e-01,-1.7812e+00,1.0800e-01],[-5.0000e+00,1.5218e+00,3.9229e+00,-2.6204e+00],[-3.5000e+00,1.6000e-03,1.1477e+00,1.0620e-01],[-3.5000e+00,-1.5540e-01,1.4062e+00,3.9880e-01],[-3.4000e+00,-1.0780e-01,-1.0750e+00,1.5702e+00],[-4.0000e+00,2.5000e-02,3.8440e-01,2.6560e-01],[-7.2312e+00,4.0500e-01,1.3350e+01,6.2340e-01]])    
    ecoeff=np.array([[0.0356,-0.1246,-0.5718,0.9938],[0.2624,0.0672,-0.219,-0.4285],[-0.0156,0.1597,0.4199,-0.5562],[0.4659,-0.3296,-0.0876,-0.0329],[0.0032,0.0766,-0.0656,-0.1294],[-0.0672,0.4016,0.3017,-0.4844],[1.0468,-0.3788,-2.4517,1.4656],[1.5,-0.6426,1.8564,0.5636]])
    
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    altitude = 90-zen
    zen = zen * deg2rad
    azimuth = azimuth * deg2rad
    altitude = altitude * deg2rad
    Idh = radD
    # Ibh = radI/sin(altitude)
    Ibn = radI

    # Skyclearness
    PerezClearness = ((Idh+Ibn)/(Idh+1.041*np.power(zen, 3)))/(1+1.041*np.power(zen, 3))
    # Extra terrestrial radiation
    day_angle = jday*2*np.pi/365
    #I0=1367*(1+0.033*np.cos((2*np.pi*jday)/365))
    I0 = 1367*(1.00011+0.034221*np.cos(day_angle) + 0.00128*np.sin(day_angle)+0.000719 *
               np.cos(2*day_angle)+0.000077*np.sin(2*day_angle))    # New from robinsson

    # Optical air mass
    # m=1/altitude; old
    if altitude >= 10*deg2rad:
        AirMass = 1/np.sin(altitude)
    elif altitude < 0:   # below equation becomes complex
        AirMass = 1/np.sin(altitude)+0.50572*np.power(180*complex(altitude)/np.pi+6.07995, -1.6364)
    else:
        AirMass = 1/np.sin(altitude)+0.50572*np.power(180*altitude/np.pi+6.07995, -1.6364)

    # Skybrightness
    # if altitude*rad2deg+6.07995>=0
    PerezBrightness = (AirMass*radD)/I0
    if Idh <= 10:
        # m_a=0;m_b=0;m_c=0;m_d=0;m_e=0;
        PerezBrightness = 0
    if altitude < 0:
        print("Airmass")
        print(AirMass)
        print(PerezBrightness)
    # sky clearness bins
    if PerezClearness < 1.065:
        intClearness = 0
    elif PerezClearness < 1.230:
        intClearness = 1
    elif PerezClearness < 1.500:
        intClearness = 2
    elif PerezClearness < 1.950:
        intClearness = 3
    elif PerezClearness < 2.800:
        intClearness = 4
    elif PerezClearness < 4.500:
        intClearness = 5
    elif PerezClearness < 6.200:
        intClearness = 6
    elif PerezClearness > 6.200:
        intClearness = 7
    else:
        raise ValueError('No valid PerezClearness, are inputs NaN?')

    m_a = acoeff[intClearness,  0] + acoeff[intClearness,  1] * zen + PerezBrightness * (acoeff[intClearness,  2] + acoeff[intClearness,  3] * zen)
    m_b = bcoeff[intClearness,  0] + bcoeff[intClearness,  1] * zen + PerezBrightness * (bcoeff[intClearness,  2] + bcoeff[intClearness,  3] * zen)
    m_e = ecoeff[intClearness,  0] + ecoeff[intClearness,  1] * zen + PerezBrightness * (ecoeff[intClearness,  2] + ecoeff[intClearness,  3] * zen)

    if intClearness > 0:
        m_c = ccoeff[intClearness, 0] + ccoeff[intClearness, 1] * zen + PerezBrightness * (ccoeff[intClearness, 2] + ccoeff[intClearness, 3] * zen)
        m_d = dcoeff[intClearness, 0] + dcoeff[intClearness, 1] * zen + PerezBrightness * (dcoeff[intClearness, 2] + dcoeff[intClearness, 3] * zen)
    else:
        # different equations for c & d in clearness bin no. 1,  from Robinsson
        m_c = np.exp(np.power(PerezBrightness * (ccoeff[intClearness, 0] + ccoeff[intClearness, 1] * zen), ccoeff[intClearness, 2]))-1
        m_d = -np.exp(PerezBrightness * (dcoeff[intClearness, 0] + dcoeff[intClearness, 1] * zen)) + dcoeff[intClearness, 2] + \
            PerezBrightness * dcoeff[intClearness, 3] * PerezBrightness

    skyvaultalt = np.atleast_2d([])
    skyvaultazi = np.atleast_2d([])

    if True:
        # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
        skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
        skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
        for j in range(7):
            for k in range(1, int(360/skyvaultaziint[j]) + 1):
                skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])
                skyvaultazi = np.append(skyvaultazi, k*skyvaultaziint[j])

        skyvaultalt = np.append(skyvaultalt, 90)
        skyvaultazi = np.append(skyvaultazi, 360)

    skyvaultzen = (90 - skyvaultalt) * deg2rad
    skyvaultalt = skyvaultalt * deg2rad
    skyvaultazi = skyvaultazi * deg2rad

    # Angular distance from the sun from Robinsson
    cosSkySunAngle = np.sin(skyvaultalt) * np.sin(altitude) + \
                     np.cos(altitude) * np.cos(skyvaultalt) * np.cos(np.abs(skyvaultazi-azimuth))

    # Main equation
    lv = (1 + m_a * np.exp(m_b / np.cos(skyvaultzen))) * ((1 + m_c * np.exp(m_d * np.arccos(cosSkySunAngle)) +
                                                           m_e * cosSkySunAngle * cosSkySunAngle))

    # Normalisation
    lv = lv / np.sum(lv)

    #x = np.atleast_2d([])
    #lv = np.transpose(np.append(np.append(np.append(x, skyvaultalt*rad2deg), skyvaultazi*rad2deg), lv))
    x = np.transpose(np.atleast_2d(skyvaultalt*rad2deg))
    y = np.transpose(np.atleast_2d(skyvaultazi*rad2deg))
    z = np.transpose(np.atleast_2d(lv))
    lv = np.append(np.append(x, y, axis=1), z, axis=1)
    return lv, PerezClearness, PerezBrightness

def Lside_veg_v2020a(azimuth,altitude,Ta,Tw,SBC,ewall,Ldown,esky,t,F_sh,CI,Lup):
    
    # This m-file is the current one that estimates L from the four cardinal points 20100414
    #Building height angle from svf
    svfalfa=np.arcsin(np.exp((np.log(1-0.6))/2))
    
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t
    
    F_sh = 2*F_sh-1  #(cylindric_wedge scaled 0-1)
    
    c=1-CI
    Lsky_allsky = esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4)
    
    viktveg, viktwall, viktsky, viktrefl = -1.1102230246251565e-16, 0.719056214891864, 0.28094378510813617, 0.7190562148918639
    alfaB=np.arctan(svfalfa)
    betaB=np.arctan(np.tan((svfalfa)*F_sh))
    betasun=((alfaB-betaB)/2)+betaB
    Lsky = (0.6 * Lsky_allsky) * viktsky * 0.5
    Lrefl = (Ldown + Lup) * (viktrefl) * (1 - ewall) * 0.5
    Lground = Lup * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5

    if altitude > 0:  # daytime
        # betasun = np.arctan(0.5*np.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180-t))  and  (azimuth <= (360-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziE*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    
    Lsum = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    if altitude>0: # daytime
        if (azimuth <= (90-t))  or  (azimuth > (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziS*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum+= Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (360-t))  or  (azimuth <= (180-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziW*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        if (azimuth > (90-t))  and  (azimuth <= (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziN*(np.pi/180)))**4)*viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        
        Lsum += Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

        return Lsum
    return Lsum*4

def Kside_veg_v2019a(radI, radD, radG,altitude, albedo, F_sh, Kup, lv, diffsh):
    # New reflection equation 2012-05-25
    
    deg2rad = np.pi / 180
    KsideD = 0

    # Direct radiation ###
    # Kside with cylinder ###
    KsideI = radI * np.cos(altitude * deg2rad)

    ### Diffuse and reflected radiation ###
#    viktveg, viktwall = 4.440892098500626e-16, 0.7190562148918634
    svfviktbuveg = 0.7190562148918634# + (4.440892098500626e-16) * (1 - 1)

    ### Anisotropic Diffuse Radiation after Perez et al. 1993 ###

    aniAlt = lv[0][:, 0]
    aniLum = lv[0][:, 2]

    phiVar = np.zeros((145, 1))
    radTot = np.zeros(1)

    for ix in range(0, 145):  # Azimuth delta
        if ix < 60:
            aziDel = 12
        elif ix >= 60 and ix < 108:
            aziDel = 15
        elif ix >= 108 and ix < 126:
            aziDel = 20
        elif ix >= 126 and ix < 138:
            aziDel = 30
        elif ix >= 138 and ix < 144:
            aziDel = 60
        elif ix == 144:
            aziDel = 360

        phiVar[ix] = (aziDel * deg2rad) * (np.sin((aniAlt[ix] + 6) * deg2rad) - np.sin((aniAlt[ix] - 6) * deg2rad))  # Solid angle / Steradian
        radTot = radTot + (aniLum[ix] * phiVar[ix] * np.sin(aniAlt[ix] * deg2rad))  # Radiance fraction normalization
    lumChi = (aniLum * radD) / radTot  # Radiance fraction normalization
    for idx in range(0, 145):
            anglIncC = np.cos(aniAlt[idx] * deg2rad)  # Angle of incidence, np.cos(0) because cylinder - always perpendicular
            KsideD += diffsh[idx] * lumChi[idx] * anglIncC * phiVar[idx]  # Diffuse vertical radiation
    Keast = (albedo * (svfviktbuveg * (radG * (1 - F_sh) + radD * F_sh)) + Kup) * 0.5
    return Keast, KsideI, KsideD
def cylindric_wedge(zen, svfalfa):

    # Fraction of sunlit walls based on sun altitude and svf wieghted building angles
    # input: 
    # sun zenith angle "beta"
    # svf related angle "alfa"

    beta=zen
    alfa= svfalfa
    
    # measure the size of the image
    # sizex=size(svfalfa,2)
    # sizey=size(svfalfa,1)
    
    xa=1-2./(np.tan(alfa)*np.tan(beta))
    ha=2./(np.tan(alfa)*np.tan(beta))
    ba=(1./np.tan(alfa))
    hkil=2.*ba*ha
    if xa<0:
        qa=np.tan(beta)/2
        Za=((ba**2)-((qa**2)/4))**0.5
        phi=np.arctan(Za/qa)
        A=(np.sin(phi)-phi*np.cos(phi))/(1-np.cos(phi))
        ukil=2*ba*xa*A
    else:
        qa,Za,phi,A,ukil=0,0,0,0,0    
    
    Ssurf=hkil+ukil
    
    F_sh=(2*np.pi*ba-Ssurf)/(2*np.pi*ba)#Xa
    
    return F_sh
def diffusefraction(radG,altitude,Kt,Ta,RH):
    """
    This function estimates diffuse and directbeam radiation according to
    Reindl et al (1990), Solar Energy 45:1

    :param radG:
    :param altitude:
    :param Kt:
    :param Ta:
    :param RH:
    :return:
    """

    alfa = altitude*(np.pi/180)

    if Ta <= -999.00 or RH <= -999.00 or np.isnan(Ta) or np.isnan(RH):
        if Kt <= 0.3:
            radD = radG*(1.020-0.248*Kt)
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG*(1.45-1.67*Kt)
        else:
            radD = radG*0.147
    else:
        RH = RH/100
        if Kt <= 0.3:
            radD = radG*(1 - 0.232 * Kt + 0.0239 * np.sin(alfa) - 0.000682 * Ta + 0.0195 * RH)
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG*(1.329- 1.716 * Kt + 0.267 * np.sin(alfa) - 0.00357 * Ta + 0.106 * RH)
        else:
            radD = radG*(0.426 * Kt - 0.256 * np.sin(alfa) + 0.00349 * Ta + 0.0734 * RH)

    radI = (radG - radD)/(np.sin(alfa))

    # Corrections for low sun altitudes (20130307)
    if radI < 0:
        radI = 0
    if altitude < 1 and radI > radG:
        radI=radG
    if radD > radG:
        radD = radG
    return radI, radD

def daylen(doy, xlat):
    # Calculation of declination of sun (Eqn. 16). Amplitude= +/-23.45
    # deg. Minimum = DOY 355 (DEC 21), maximum = DOY 172.5 (JUN 21/22).
    # Sun angles.  SOC limited for latitudes above polar circles.
    # Calculate daylength, sunrise and sunset (Eqn. 17)
    #doy is day of year, xlat is the latitude in degrees
    rad=np.pi/180.0
    DEC = -23.45 * np.cos(2.0*np.pi*(doy+10.0)/365.0)
    SOC = np.tan(rad*DEC) * np.tan(rad*xlat)
    SOC = min(max(SOC,-1.0),1.0) #correction for polar circles
    # SOC=alt
    day_len = 12.0 + 24.0*np.arcsin(SOC)/np.pi
    sun_upp = 12.0 - day_len/2.0 
    sun_down = 12.0 + day_len/2.0
    return day_len, DEC, sun_down, sun_upp

def Solweig1D_2020a_calc(albedo_b, absK, absL, ewall, Fside, Fup, Fcyl, altitude, azimuth, zen, jday,
                         onlyglobal, location, dectime, altmax, Ta, RH, radG, radD, radI, P, TgK, Tstart, albedo_g, eground, TgK_wall, Tstart_wall, TmaxLST, TmaxLST_wall,
                         svfalfa,  svfbuveg, CI, diffsh):
    # This is the core function of the SOLWEIG1D model, 2019-Jun-21
    # Fredrik Lindberg, fredrikl@gvc.gu.se, Goteborg Urban Climate Group, Gothenburg University, Sweden
    # Instrument offset in degrees
    t = 0.
    # Stefan Bolzmans Constant
    SBC = 5.67051e-8
    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location['latitude'])
    # Vapor pressure
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH / 100.)
    # Determination of clear - sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (1 - (1 + msteg) * np.exp(-((1.2 + 3.0 * msteg) ** 0.5))) + 0  # -0.04 old error from Jonsson et al.2006

    if altitude > 0: # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt = clearnessindex_2013b(zen, jday, Ta, RH / 100., radG, location, P)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if onlyglobal == 1:
            I0, CI, Kt = clearnessindex_2013b(zen, jday, Ta, RH / 100., radG, location, P)
            if (CI > 1) or (CI == np.inf):
                CI = 1

            radI, radD = diffusefraction(radG, altitude, Kt, Ta, RH)

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        if True:
            zenDeg = zen*(180/np.pi)
            lv = Perez_v3(zenDeg, azimuth, radD, radI, jday)   # Relative luminance

            aniLum = 0.
            for idx in range(0, 145):
                aniLum = aniLum + diffsh[idx] * lv[0][idx][2]     # Total relative luminance from sky into each cell

            dRad = aniLum * radD   # Total diffuse radiation from sky into each cell


        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of  dfm
        Tgamp = (TgK * altmax - Tstart) + Tstart
        Tgampwall = (TgK_wall * altmax - (Tstart_wall)) + (Tstart_wall)
        Tg = Tgamp * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST / 24 - SNUP / 24)) * np.pi / 2) + Tstart # 2015 a, based on max sun altitude
        Tgwall = Tgampwall * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST_wall / 24 - SNUP / 24)) * np.pi / 2) + (Tstart_wall) # 2015a, based on max sun altitude

        if Tgwall < 0:  # temporary for removing low Tg during morning 20130205
            # Tg = 0
            Tgwall = 0

        # New estimation of Tg reduction for non - clear situation based on Reindl et al.1990
        radI0, _ = diffusefraction(I0, altitude, 1., Ta, RH)
        corr = 0.1473 * np.log(90 - (zen / np.pi * 180)) + 0.3454  # 20070329 correction of lat, Lindberg et al. 2008
        CI_Tg = (radI / radI0) + (1 - corr)
        if (CI_Tg > 1) or (CI_Tg == np.inf):
            CI_Tg = 1
        Tg = Tg * CI_Tg  # new estimation
        Tgwall = Tgwall * CI_Tg

        if Tg < 0.:
            Tg = 0.
        #Tg[Tg < 0] = 0  # temporary for removing low Tg during morning 20130205

        gvf = 1

        Lup = SBC * eground * ((gvf + Ta + Tg + 273.15) ** 4)

        # Building height angle from svfs
        F_sh = float(cylindric_wedge(zen, svfalfa))  # Fraction shadow on building walls based on sun alt and svf
        #F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * np.sin(altitude * (np.pi / 180)) + dRad + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))
        
        Kup = albedo_g * (radI * np.sin(altitude * (np.pi / 180.))) + dRad + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)

        Kesnw, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, altitude, albedo_b, F_sh, Kup, lv, diffsh)
    else:  # # # # # # # NIGHTTIME # # # # # # # #
        Tgwall = 0
        # Nocturnal K fluxes set to 0
        Knight = 0.
        Kdown = 0.
        Kup = 0.
        Kesnw = 0.#Keast, Kwest, Knorth, Ksouthö
        KsideI = 0.
        KsideD = 0.
        F_sh = 0.
        Tg = 0.

        # # # # Lup # # # #
        Lup = SBC * eground * ((Knight + Ta + Tg + 273.15) ** 4)
        I0 = 0

    # # # # Ldown # # # # 
    Ldown = (0.6 + 1 - 1) * esky * SBC * ((Ta + 273.15) ** 4) + \
            (2 - 1 - 1) * ewall * SBC * ((Ta + 273.15) ** 4) + \
            (1 - 0.6) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + \
            (2 - 0.6 - 1) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)

    if CI < 0.95:  # non - clear conditions
        c = 1 - CI
        Ldown = Ldown * (1 - c) + \
                c * (0.6 * SBC * ((Ta + 273.15) ** 4) + (1 - 0.6) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) +
                     (1 - 0.6) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4))

    # # # # Lside # # # # FIXAD
    Lsum = Lside_veg_v2020a(azimuth, altitude, Ta, Tgwall, SBC, ewall, Ldown,esky, t, F_sh, CI, Lup)

    # # # # Calculation of radiant flux density and Tmrt # # # #
    # Human body considered as a cylinder with Perez et al. (1993)
    Sstr = absK * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (4*Kesnw) * Fside) + absL * ((Ldown + Lup) * Fup + (Lsum) * Fside) #4 directions

    Tmrt = float(np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.15)
    return Tmrt, CI
