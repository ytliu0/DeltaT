# &Delta;T Calculator

&Delta;T is the difference between the terrestrial time (TT) and universal time (UT1). See my [dynamical time](http://ytliu.epizy.com/eclipse/dynamical_time.html) page for a brief explanation of these terms.

This package provides a python code to implement the &Delta;T calculation according to the fitting and extrapolation formulae in [Stephenson et al (2016)](https://royalsocietypublishing.org/doi/10.1098/rspa.2016.0404) and [Morrison et al (2021)](https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776). Specifically, values of &Delta;T from -720 to 2022 are computed using their [spline fit cubic polynomials](http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt), but the coefficients after 2013 have been modified to extend the fitting formulae to 2022 (their original formulae are fitted up to 2019). Outside this range &Delta;T is extrapolated by integrating their long-term lod function. The uncertainty in &Delta;T is estimated based on [their tables](http://astro.ukho.gov.uk/nao/lvm/) for years between -2000 to +2500. Outside this range, the error is estimated using quadratic functions. However, they are probably not reliable. 

The cubic spline polynomials are fitted to the &Delta;T data from -720 to 2022. Daily values of &Delta;T in 2022 and later can be constructed from the data in [IERS bulletins](https://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html), which can be used to calculate the error of the extrapolation formula after 2022 and to extend the fitting formula. The csv file `DeltaT_IERS.csv` contains the &Delta;T data from 1972 to the recent time computed from the IERS data.

All functions that compute &Delta;T are in `DeltaT.py`. `DeltaT(y)` computes &Delta;T (in seconds) as a function of year y using the fitting and extrapolation formulae. `DeltaT_error_estimate(y)` provides an error estimate of &Delta;T. `DeltaT_with_error_estimate(y)` combines the two functions and returns &Delta;T(y) with an error estimate. `DeltaT_hybrid(jd)` calculates &Delta;T by linear interpolation from the IERS data if the Julian date `jd` is in inside the range of the IERS data. Outside the data range, the function computes &Delta;T by the cubic spline polynomials or by the integrated long-term lod function. `DeltaT_hybrid_error_estimate(jd)` estimates the error of the computed &Delta;T. Inside the IERS data range, the error is estimated by either the measurement error of UT1 or half of the difference between two values of &Delta;T used in the interpolation, whichever is larger.  The input `y` and `jd` can be a scalar or a 1D array. 

The functions in `regression_cubic_spline.py` can be used to update the current fitting formulae.

The notebook `DeltaT_examples.ipynb` shows examples of using these functions.

A javascript version of the code is used in my [eclipse website](http://ytliu.epizy.com/eclipse/).

## References

- F.R. Stephenson, L.V. Morrison, and C.Y. Hohenkerk, [Measurement of the Earth's rotation: 720 BC to AD 2015](https://royalsocietypublishing.org/doi/10.1098/rspa.2016.0404), *Proc. R. Soc. A.*, **472**:20160404 (2016).

- L.V. Morrison, F.R. Stephenson, C.Y. Hohenkerk, and M. Zawilski, [Addendum 2020 to 'Measurement of the Earth's rotation: 720 BC to AD 2015'](https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776), *Proc. R. Soc. A.*, **477**:20200776 (2021).
