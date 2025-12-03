# BiomedicalSignalProcessing-Project1
Processing of ECG signals.


# Algorithm notes

* 1) NOISE FILTER
- eliminates baseline drift and powerline interference

- y(k) = 1/(K*K) * sigma\[m = k - K + 1, k\] sigma\[n = m - K + 1, m\] x(n) - 1/(L*L) * sigma\[m = k - L + 1, k\] sigma\[n = m - L + 1, m\] x(n) 
- K and L are constants


* 2) DIFFERENTIATOR
- accounts for steep slopes, calcs position of R peak
- band limited differentiator
- h_d = 1/3 (-1, -2, 0, 2, 1)



* 3) ENERGY COLLECTOR
- MULTIPLICATOR - squarer
- ya(k) = {xa(k)}^2
- convolution in freq domain

- INTEGRATOR
- yb(k) = sigma\[n=k - N + 1, k\] xb(n)
- N is number of samples to integrate

* 4) ADAPTIVE MIN DIST CLASSIFIER
- muQRS(k) = 0.9muQRS(k-1) + 0.1yb
- muNQRS(k) = 0.9muNQRS(k-1) + 0.1yb



5) MINMAX SEARCHER
- inputs DIFFERENTIATOR and MIN DIST CLASSIFIER




# FIX find_energy_peaks