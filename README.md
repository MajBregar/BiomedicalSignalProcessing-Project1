# BiomedicalSignalProcessing-Project1
Processing of ECG signals.


# Algorithm notes


# 1) preprocess, get bandpass filtered signal and integrated signal
# 2) DETECTION


# a - detect a fiducial mark (candidate for R wave)
* R wave at same location as rising edge of integrated signal (max slope or peak of rwave)


# b - threshold adjusting
* pI - peak value of integrated signal
* pF - peak value of bandpass filtered signal

* pI process
    * SPKI - estimate of the pI signal peaks 
    * NPKI - estimate of the pI noise peaks (any peak thats not a QRS)

    * when detected new peak -> classify as NOISE or SIGNAL peak
    * peak is SIGNAL <=>  (peak > THR_I1) OR (peak > THR_I2 in case of searchback) 
    * peak is NOISE <=> !SIGNAL

    * if pI is SIGNAL => SPKI = 0.125 pI + 0.875 SPKI
    * if pI is NOISE => NPKI = 0.125 PEAKI + 0.875 NPKI
    * THRESHOLD_I1 = NPKI + 0.25 (SPKI - NPKI)
    * THRESHOLD I2 = 0.5 THRESHOLD I1
    * if peak is SIGNAL AND searchback was required => SPKI = 0.25 pI + 0.75 SPKI

* pF process    
    * SPKF - estimate of the pF signal peaks
    * NPKF - estimate of the pF noise peaks (any peak thats not a QRS)
 
    * if pF is SIGNAL => SPKF = 0.125 pF + 0.875 SPKF
    * if pF is NOISE => NPKF = 0.125 pF + 0.875 NPKF
    * THRESHOLD_F1 = NPKF + 0.25 (SPKF - NPKF)
    * THRESHOLD_F2 = 0.5 THRESHOLD_F1
    * if peak is SIGNAL AND searchback was required => SPKF = 0.25 pF + 0.75 SPKF.

* irregular heart rates 
    * THRESHOLD_I1 *= 0.5
    * THRESHOLD_F1 *= 0.5


* peak is QRS (SIGNAL) <=> pI is SIGNAL AND pF is SIGNAL


# c - AVG_RR interval and rate limits
* RR_AVERAGE1 = 0.125 (RR(i), RR(i-1), ..., RR(i-7)) //avg of 8 latest beats
* RR_AVERAGE2 = 0.125 (RRS(i), RRS(i-1), ..., RRS(i-7)) //avg of 8 latest beats that fell into rate limts
* rate limits:
    * RR_LOW = 0.92 RR_AVERAGE2
    * RR_HIGH = 1.16 RR_AVERAGE2
    * RR_MISSED = 1.66 RR_AVERAGE2
* if no QRS found within interval RR_MISSED => QRS candidate = max_peak(thr_1, thr_2)
* if (all i: RR_LOW <= RR(i) <= RR_HIGH) => RR_AVERAGE2 = RR_AVERGE1 (regular rate)


# d - T-wave identification
* if (RR(i) < 360 ms) AND (RR(i) >= 200 ms) AND (pF < LatestQRSF) => T Wave




# PIPELINE 
## 1. - preprocess to aquire BP filtered signal and MWI integrated signal
## 2. - for 2 seconds stay in LEARNING PHASE (set thresholds)
## 3. EXECTUTE FOR EVERY SAMPLE - ONLINE
#   determine pF and pI from signals (just that they are maximums)
#   classfiy pF and pI as SIGNAL or NOISE
#   update the estimates (NPKI/NPKF, SPKI/SPKF)
#   if both are SIGNAL -> QRS candidate
#       set refractory period
#       if 200ms < RR < 360ms:
#           chekcking Twave (if slope < 0.5 prev_qrs_slope) TWave
#       update RR avgs and limits

#    no qrs complex: look for searchback
    


# TODO 
* make a learning phase that initializes the state THEN run it again
* see where the missing 2 delay is
* double check everything if its truly online
* evaluate everything and REMEMBER TO LOOK AT SLIDES TO SEE WHAT THE FORUMULAS ARE
