import numpy as np
import DeltaT as DT

def regression_cubic_spline(y, DeltaT, a_in, knots):
    """
    Fit a regression cubic spline on DeltaT ~ y with cubic coefficients a_in in the first interval.
    Input:
      y: 1D numpy array for years.
      DeltaT: 1D numpy array for DeltaT. y and DeltaT must have the same length.
      a_in: 1D numpy array of size 4 for the coefficients of cubic polynomials in the first interval 
      of y in [knots[0], knots[1]].
      knots: 1D numpy array storing the knots for the regression spline. Note that all values of y must be between knots[0] and knots[-1].
    Return:
      Tuple: spline, residuals, c2 
      See explanation in the doc string of extend_Morrison_etal_cubic_spline().
    """
    def h(t, m):
        """
        Function h(t, m) = 0 if t < m and h(t,m) = (t-m)^3 if t >= m
        """
        return np.where(t < m, 0, (t-m)**3)
    
    def hd(t,m, n):
        """
        Calculate the n-th derivative of h.
        """
        if n==1:
            return np.where(t < m, 0, 3*(t-m)**2)
        elif n==2:
            return np.where(t < m, 0, 6*(t-m))
        elif n==3:
            return np.where(t < m, 0, 6)
        else:
            return 0

    dy = knots[1] - knots[0]
    t = (y - knots[0])/dy # the rescaled time variable
    # Follow the method in Section 7.4 of *An Introduction to Statistical Learning* 
    # by James, Witten, Hastie and Tibshirani (https://www.statlearning.com/) to fit 
    # the data with the equation
    # DeltaT = a_in[0] + a_in[1]*t + a_in[2]*t^2 + a_in[3]*t^3 + sum_{m=1}^n c_m*h(t, m)
    # h(t, m) = 0 if t < m and (t-m)^3 if t >= m and values of m are determined by 
    # m = (knots[1:] - knots[0])/dy.
    # It's clear that the fitted DeltaT and its first two derivatives will be continuous, but 
    # its third derivative will be discontinuous at the knots.
    # Define z = DeltaT - a_in[0] + a_in[1]*t + a_in[2]*t^2 + a_in[3]*t^3. Then 
    # z = sum_{m=1}^n c_m*h(t, m)
    # Coefficients c_m can be computed by linear regression with 0 intercept.
    # Note that any values of c_m will not change the fitted z in the first interval, so we can 
    # remove data from the first interval.
    t = t[y > knots[1]]
    z = DeltaT[y > knots[1]] - a_in[0] - a_in[1]*t - a_in[2]*t*t - a_in[3]*t**3
    nc = len(knots)-2
    # Set up the features
    X = np.zeros((len(t), nc))
    m = (knots[1:(nc+1)] - knots[0])/dy
    for i in range(nc):
        X[:,i] = h(t, m[i])
    Xmax = np.max(X, axis=0)
    reg_coef = np.linalg.lstsq(X/Xmax, z, rcond=None)[0]/Xmax
    # Now construct the coefficients of the cubic polynomials in the intervals 
    # between knots from the regression coefficients.
    # If f(x) = c0 + c1*(x-x0) + c2*(x-x0)^2 + c3*(x-x0)^3, then 
    # c0 = f(x0), c1 = f'(x0), c2 = f''(x0)/2, c3 = f'''(x0)/6, where ' = d/dx.
    a = np.zeros((4,nc+1))
    a[:,0] = a_in
    for i in range(nc):
        fac = (knots[i+2] - knots[i+1])/dy  # take into account of variable rescaling
        e = 1e-12    # a tiny number compared to m[i]
        a[0,i+1] = a_in[0] + a_in[1]*m[i] + a_in[2]*m[i]**2 + a_in[3]*m[i]**3 + sum(reg_coef*h(m[i]+e, m))
        a[1,i+1] = fac*(a_in[1] + 2*a_in[2]*m[i] + 3*a_in[3]*m[i]**2 + sum(reg_coef*hd(m[i]+e, m, 1)))
        a[2,i+1] = fac*fac*(a_in[2] + 3*a_in[3]*m[i] + 0.5*sum(reg_coef*hd(m[i]+e, m, 2)))
        a[3,i+1] = fac*fac*fac*(a_in[3] + sum(reg_coef*hd(m[i]+e, m, 3))/6)
    y0 = knots[0:(nc+1)]
    y1 = knots[1:(nc+2)]
    a = np.round(a, 3)  # round a to 3 decimal places

    # Calculate the residuals of the fit
    res = DeltaT - DT.spline(y, y0, y1, a[0,:], a[1,:], a[2,:], a[3,:])
    # Calculate the root mean square error of the cubic spline fit in the intervals between knots
    eps = np.zeros(nc+1)
    for i in range(nc+1):
        r = res[(y >= knots[i]) & (y < knots[i+1])]
        eps[i] = np.sqrt(np.mean(r*r))
    # Delta T at the last knot
    DeltaT_last = sum(a[:,-1])
    # integration constant c2 for extrapolation using the integrated lod function for y > knots[-1]
    c2 = DeltaT_last - DT.integrated_lod(knots[-1], 0)
    spline = {'y0':y0, 'y1':y1, 'a0':a[0,:], 'a1':a[1,:], 'a2':a[2,:], 'a3':a[3,:], 'epsilon':eps}
    residuals = {'y':y, 'res':res}
    return spline, residuals, c2

def extend_Morrison_etal_cubic_spline(knots_in):
    """
    Fit a new cubic spline polynomials in the intervals specified by the array knots_in using the IERS data in 'DeltaT_IERS.csv'.
    The cubic spline fit before knots_in[0] are not altered, but the cubic cofficients between y0[i] and knots_in[0] may be altered because of possible variable rescaling. Here y0[i] is the largest year in y0 less than knots_in[0] in the Morrison et al's cubic spline knots.
    Return:
      Tuple: spline, residuals, c2.
      spline is a dictionary with keys 'y0', 'y1', 'a0', 'a1', 'a2', 'a3', and 'epsilon'.
      Here y0 = knots[0:(nc+1)], y1 = knots[1:(nc+2)], where nc = len(knots)-2; 
      knots = np.append(y0[i], knots_in);
      a0, a1, a2, a3 are the coefficients of cubic spline polynomials so that 
      fitted DeltaT = a0 + a1*t + a2*t^2 + a3*t^3 
      with t = (y - y0)/(y1 - y0).
      Epsilon is the root mean square error of the fitting formula in the interval [y0, y1].
      The fitting coefficients are for the intervals [knots[i], knots[i+1]] with i=0, 1, 2,... nc.
      residuals is a dictionary with keys 'y' and 'res'. res is a 1D array storing the residuals of the fit at y.
      c2 is the integration constant for extrapolation using the integrated lod function for y > knots[-1].
    """
    # Search for i so that y0[i] is the largest year in the Morrison et al's cubic spline knots less than knots_in[0]
    i = np.searchsorted(DT.y0, knots_in[0], 'left')-1
    # Add y0[i] to knots_in
    knots = np.append(DT.y0[i], np.array(knots_in, dtype=np.float64))
    # The time variable in the interval [y0[i], knots_in[0]] is t = (y - y0[i])/(knots_in[0]- y0[i]), 
    # which is in general not the same as the original t = (y - y0[i])/(y1[i] - y0[i]). So the original 
    # cubic coefficients should be modified by a1[i] -> fac*a1[i], a2[i] -> fac^2*a2[i], a3[i] -> fac^3*a3[i], 
    # where fac = (knots_in[0]- y0[i])/(y1[i] - y0[i]).
    fac = (knots_in[0] - DT.y0[i])/(DT.y1[i] - DT.y0[i])
    a_init = np.array([DT.a0[i], DT.a1[i]*fac, DT.a2[i]*fac**2, DT.a3[i]*fac**3])
    data = np.loadtxt('DeltaT_IERS.csv', delimiter=',', skiprows=1)
    y = (data[:,3] - 51544)/365.2425 + 2000
    yrange= (y >= knots[0]) & (y < knots[-1])
    # filter y and DeltaT to the range where y is in [knots[0], knots[-1])
    y = y[yrange]
    DeltaT = data[yrange, 6]
    return regression_cubic_spline(y, DeltaT, a_init, knots)
