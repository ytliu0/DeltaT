import numpy as np

# Constant arrays for cubic spline polynomials, 
# see http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt
# The coefficients after 2013 have been modified to include data after 2019.
y0 = np.array([-720, -100, 400, 1000, 1150, 1300, 1500, 1600, 1650, 1720, 1800,
             1810, 1820, 1830, 1840, 1850, 1855, 1860, 1865, 1870, 1875, 1880,
             1885, 1890, 1895, 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935,
             1940, 1945, 1950, 1953, 1956, 1959, 1962, 1965, 1968, 1971, 1974,
             1977, 1980, 1983, 1986, 1989, 1992, 1995, 1998, 2001, 2004, 2007,
             2010, 2013, 2016, 2019])
y1 = np.array([-100, 400, 1000, 1150, 1300, 1500, 1600, 1650, 1720, 1800, 1810,
             1820, 1830, 1840, 1850, 1855, 1860, 1865, 1870, 1875, 1880, 1885,
             1890, 1895, 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940,
             1945, 1950, 1953, 1956, 1959, 1962, 1965, 1968, 1971, 1974, 1977,
             1980, 1983, 1986, 1989, 1992, 1995, 1998, 2001, 2004, 2007, 2010,
             2013, 2016, 2019, 2022])
a0 = np.array([20371.848, 11557.668, 6535.116, 1650.393, 1056.647, 681.149, 292.343,
             109.127, 43.952, 12.068, 18.367, 15.678, 16.516, 10.804, 7.634,
             9.338, 10.357, 9.04, 8.255,
             2.371, -1.126, -3.21, -4.388, -3.884, -5.017, -1.977, 4.923, 11.142,
             17.479, 21.617, 23.789, 24.418, 24.164, 24.426, 27.05, 28.932,
             30.002, 30.76, 32.652, 33.621, 35.093, 37.956, 40.951, 44.244,
             47.291, 50.361, 52.936, 54.984, 56.373, 58.453, 60.678, 62.898,
             64.083, 64.553, 65.197, 66.061, 66.919, 68.128, 69.248])
a1 = np.array([-9999.586, -5822.27, -5671.519, -753.21, -459.628, -421.345,
             -192.841, -78.697, -68.089, 2.507, -3.481, 0.021, -2.157, -6.018,
             -0.416, 1.642, -0.486, -0.591, -3.456, -5.593, -2.314, -1.893, 0.101,
             -0.531, 0.134, 5.715, 6.828, 6.33, 5.518, 3.02, 1.333, 0.052, -0.419,
             1.645, 2.499, 1.127, 0.737, 1.409, 1.577, 0.868, 2.275, 3.035, 3.157,
             3.199, 3.069, 2.878, 2.354, 1.577, 1.648, 2.235, 2.324, 1.804, 0.674,
             0.466, 0.804, 0.839, 1.005, 1.341, 0.620])
a2 = np.array([776.247, 1303.151, -298.291, 184.811, 108.771, 61.953, -6.572,
             10.505, 38.333, 41.731, -1.126, 4.629, -6.806, 2.944, 2.658, 0.261,
             -2.389, 2.284, -5.148, 3.011, 0.269, 0.152, 1.842, -2.474, 3.138,
             2.443, -1.329, 0.831, -1.643, -0.856, -0.831, -0.449, -0.022, 2.086,
             -1.232, 0.22, -0.61, 1.282, -1.115, 0.406, 1.002, -0.242, 0.364,
             -0.323, 0.193, -0.384, -0.14, -0.637, 0.708, -0.121, 0.21, -0.729,
             -0.402, 0.194, 0.144, -0.109, 0.275, 0.061, -0.782])
a3 = np.array([409.16, -503.433, 1085.087, -25.346, -24.641, -29.414, 16.197, 3.018,
             -2.127, -37.939, 1.918, -3.812, 3.25, -0.096, -0.539, -0.883, 1.558,
             -2.477, 2.72, -0.914, -0.039, 0.563, -1.438, 1.871,
             -0.232, -1.257, 0.72, -0.825, 0.262, 0.008, 0.127, 0.142, 0.702,
             -1.106, 0.614, -0.277, 0.631, -0.799, 0.507, 0.199, -0.414, 0.202,
             -0.229, 0.172, -0.192, 0.081, -0.165, 0.448, -0.276, 0.11, -0.313,
             0.109, 0.199, -0.017, -0.084, 0.128, -0.071, -0.281, 0.193])

# Integration constants for year < -720 and > 2022
c1 = 1.007739546148514 # chosen to make DeltaT continuous at y = -720
c2 = -150.3150351029286 # chosen to make DeltaT continuous at y = 2022

# Table for estimating the errors in Delta T for years in [-2000,2500] based on 
# http://astro.ukho.gov.uk/nao/lvm/
ytab = np.array([-2000, -1600, -900, -720, -700, -600, -500, -400, -300, -200, -100, 0,
        100, 200, 300, 400, 500, 700, 800, 900, 1000, 1620, 1660, 1670, 1680,
        1730, 1770, 1800, 1802, 1805, 1809, 1831, 1870, 2022.5, 2023.5, 2030, 2040,
        2050, 2100, 2200, 2300, 2400, 2500])
eps_tab = np.array([1080, 720, 360, 180, 170, 160, 150, 130, 120, 110, 100, 90, 80, 70,
           60, 50, 40, 30, 25, 20, 15, 20, 15, 10, 5, 2, 1, 0.5, 0.4, 0.3, 0.2,
           0.1, 0.05, 0.1, 1, 2, 4, 6, 10, 20, 30, 50, 100])
nytab = len(ytab)

def spline(y, y0=y0, y1=y1, a0=a0, a1=a1, a2=a2, a3=a3):
    """
    Calculate Delta T by cubic spline polynomial. 
    y can be a scalar or a 1D numpy array.
    """
    i = np.searchsorted(y0, y, 'right')-1
    t = (y - y0[i])/(y1[i]-y0[i])
    return a0[i] + t*(a1[i] + t*(a2[i] + t*a3[i]))

def integrated_lod(y, C):
    """
    Integrated lod (deviation of mean solar day from 86400s) equation from 
    http://astro.ukho.gov.uk/nao/lvm/:
    lod = 1.72 t − 3.5 sin(2*pi*(t+0.75)/14) in ms/day, where t = (y - 1825)/100
    Using 1ms = 1e-3s and 1 Julian year = 365.25 days,
    lod = 0.62823*t - 1.278375*sin(2*pi/14*(t + 0.75) in s/year
    Integrate the equation gives
    C + 31.4115*t^2 + 894.8625/pi*cos(2*pi/14*(t + 0.75))
    in seconds. C is the integration constant.
    y can be a scalar or a 1D numpy array.
    """
    t = 0.01*(y - 1825)
    return C + 31.4115*t*t + 284.8435805251424*np.cos(0.4487989505128276*(t + 0.75))

def initialize_IERS_variables(s):
    data = np.loadtxt('DeltaT_IERS.csv', delimiter=',', skiprows=1)
    if s > 1:
       # reduce data 
       n = data.shape[0]
       data = data[range(0,n,s),:]
    jd_beg = 2400000.5 + data[0,3]
    jd_end = 2400000.5 + data[-1,3]
    y = 2000 + (data[-1,3] - 51544)/365.2425
    c = data[-1,6] - integrated_lod(y,0)
    return jd_beg, jd_end, c, data[:,6], data[:,5]

# Global variables for interpolating the IERS data file
s_IERS = 1   # Number of days between two data points (if > 1, data will be reduced)
jd_beg, jd_end, c_IERS, DeltaT_IERS, sigma = initialize_IERS_variables(s_IERS)
nIERS = len(sigma)

def DeltaT(y):
    """
    Compute Delta T using the fitting and extrapolation formulae by 
    Stephenson et al (2016) and Morrison et al (2021). See 
    http://astro.ukho.gov.uk/nao/lvm/
    The input y can be a scalar or a 1D array.
    Return Delta T in seconds. If y is a 1D array, Delta T is a 1D numpy array.
    """
    if isinstance(y, (int,float)):
       if y < -720:
           return integrated_lod(y, c1)
       if y > 2022:
           return integrated_lod(y, c2)
       return spline(y)
    else:
       y = np.array(y)
       return np.where(y < -720, integrated_lod(y, c1), 
                       np.where(y <= 2022, spline(y), integrated_lod(y, c2)) )

def DeltaT_error_estimate(y):
    """
    Estimate the error of Delta T based on the tables in http://astro.ukho.gov.uk/nao/lvm/
    The table only gives the error estimate for y in [-2000,2500]. The error outside this 
    range is estimated by quadratic functions, but they are probably not reliable.
    The input y can be a scalar or a 1D array.
    Return the result in seconds. If y is a 1D array, the result is a 1D numpy array.
    """
    k1 = 0.74e-4
    k2 = 2.2e-4

    if isinstance(y, (int,float)):
       if y < ytab[0]:
           return k1*(y-1825)**2
       if y >= ytab[nytab-1]:
           return k2*(y-1825)**2
       return eps_tab[np.searchsorted(ytab, y, 'right')-1]
    else:
       y = np.array(y)
       return np.where(y < ytab[0], k1*(y-1875)**2, 
                       np.where(y < ytab[nytab-1], eps_tab[np.searchsorted(ytab, y, 'right')-1],
                         k2*(y-1875)**2) )

def DeltaT_with_error_estimate(y):
    """
    Compute Delta T using the fitting and extrapolation formulae by 
    Stephenson et al (2016) and Morrison et al (2021) and provides an error estimate.
    The input y can be a scalar or a 1D array.
    Return a string if y is a scalar, and a 1D array of strings if y is a 1D array
    """
    if isinstance(y, (int,float)):
       dT = DeltaT(y)
       eps = DeltaT_error_estimate(y)
       if eps > 10:
           eps = round(eps)
           dT = round(dT)
       elif eps < 0.09:
           dT = round(dT, 2)
       elif eps < 0.9:
           dT = round(dT, 1)
       else:
           dT = round(dT)
       return str(dT)+u' \u00B1 '+str(eps)+' seconds'
    else:
       out = ['']*len(y)
       y = np.array(y)
       dT = DeltaT(y)
       eps = DeltaT_error_estimate(y)
       for i,e in enumerate(eps):
           if e > 10:
             out[i] = str(round(dT[i]))+u' \u00B1 '+str(round(e))+' seconds'
           elif e < 0.09:
             out[i] = str(round(dT[i],2))+u' \u00B1 '+str(e)+' seconds'
           elif e < 0.9:
             out[i] = str(round(dT[i],1))+u' \u00B1 '+str(e)+' seconds'
           else:
             out[i] = str(round(dT[i]))+u' \u00B1 '+str(e)+' seconds'
       return out
    
def DeltaT_interpolate_IERS(jd):
   """
   Calculate Delta T by linear interpolating the IERS Delta T data.
   If jd is out of range of the DeltaT_IERS array, return -99999.
   """
   if isinstance(jd, (int,float)):
      d = (jd - jd_beg)/s_IERS
      i = np.int_(d)
      d -= np.floor(d)
      dT = -99999
      if jd >= jd_beg and jd < jd_end:
        dT = DeltaT_IERS[i] + d*(DeltaT_IERS[i+1] - DeltaT_IERS[i])
      return dT
   jd = np.array(jd)
   d = (jd - jd_beg)/s_IERS
   i = np.int_(d)
   d -= np.floor(d)
   out_of_range = (jd < jd_beg) | (jd >= jd_end)
   i = np.where(out_of_range, 0, i)
   dT = DeltaT_IERS[i] + d*(DeltaT_IERS[i+1] - DeltaT_IERS[i])
   return np.where(out_of_range, -99999, dT)

def DeltaT_interpolate_IERS_error_estimate(jd):
   """
   Estimate the error of Delta T computed from DeltaT_interpolate_IERS(). 
   If jd is out of range of the DeltaT_IERS array, return -99999.
   """
   if isinstance(jd, (int,float)):
      i = np.int_((jd - jd_beg)/s_IERS)
      eps = -99999
      if jd >= jd_beg and jd < jd_end:
        eps = max(0.5*np.abs(DeltaT_IERS[i+1] - DeltaT_IERS[i]), sigma[i])
      return eps
   jd = np.array(jd)
   i = np.int_((jd - jd_beg)/s_IERS)
   out_of_range = (jd < jd_beg) | (jd >= jd_end)
   i = np.where(out_of_range, 0, i)
   hd = 0.5*np.abs(DeltaT_IERS[i+1] - DeltaT_IERS[i])
   eps = np.where(sigma[i] > hd, sigma[i], hd)
   return np.where(out_of_range, -99999, eps)

def jdy(jd):
   """
   Convert jd to y
   """
   if isinstance(jd, (int,float)):
      return 2000 + (jd - 2451544.5)/365.2425 if jd >= 2299160.5 else (jd+0.5)/365.25 - 4712
   return np.where(jd >= 2299160.5, 2000 + (jd - 2451544.5)/365.2425, (jd+0.5)/365.25 - 4712)

def DeltaT_hybrid(jd):
   """
   Compute Delta T at Julian date jd using a hybrid method:
   If jd is in the range of IERS Delta T table, compute Delta T by linear interpolating the IERS Delta T data.
   If jd is outside the range of the IERS data, compute Delta T by cubic spline interpolation or extrapolation by the integrated lod.
   """
   if isinstance(jd, (int,float)):
      if jd >= jd_end:
         return integrated_lod(jdy(jd), c_IERS)
      elif jd >= jd_beg:
         return DeltaT_interpolate_IERS(jd)
      else:
         return DeltaT(jdy(jd))
   else:
      jd = np.array(jd)
      y = jdy(jd)
      return np.where(jd >= jd_end, integrated_lod(y, c_IERS), 
               np.where(jd >= jd_beg, DeltaT_interpolate_IERS(jd), DeltaT(y)) )
   
def DeltaT_hybrid_error_estimate(jd):
   if isinstance(jd, (int,float)):
      if jd >= jd_end or jd < jd_beg:
        return DeltaT_error_estimate(jdy(jd))
      else:
        return DeltaT_interpolate_IERS_error_estimate(jd)
   else:
    jd = np.array(jd)
    y = jdy(jd)
    return np.where((jd >= jd_end) | (jd < jd_beg),  DeltaT_error_estimate(y), 
                    DeltaT_interpolate_IERS_error_estimate(jd))
   
def DeltaT_hybrid_with_error_estimate(jd):
   dT = DeltaT_hybrid(jd)
   eps = DeltaT_hybrid_error_estimate(jd)
   if isinstance(jd, (int,float)):
    if eps > 1:
        dT = round(dT)
    elif eps > 0.1:
        dT = round(dT, 1)
    elif eps > 0.01:
        dT = round(dT, 2)
    else:
       dT = round(dT,6)
    return str(dT)+u' \u00B1 '+np.format_float_positional(eps, 2, fractional=False)+' seconds'
   else:
      dT = np.where(eps > 1, np.round(dT), np.where(eps > 0.1, np.round(dT,1), np.where(eps > 0.01, np.round(dT,2), np.round(dT,6))))
      return [str(dT[i])+u' \u00B1 '+np.format_float_positional(eps[i],2, fractional=False) for i in range(len(dT))]
   