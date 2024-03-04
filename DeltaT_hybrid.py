from DeltaT import *

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
