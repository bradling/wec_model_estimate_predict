# wec_model_estimate_predict
This software, developed for my thesis, is matlab software that models WECs, allows development of estimators and disturbance predictors

# Overview


# Use Instructions
I am still working on uploading example use files demonstrating how to use this software.


## Wec modeling
A WEC model has been developed for a heaving only body, and planar motion (heave, surge, and pitch)
The WEC modeling framework extracts the hydrodynamic parameters from ANSYS Aqwa .LIS output files.
Alternatively, hydrodynamic parameters can also be formatted in a Matlab structure if they were calculated using other software.

I am currently working on uploading example use files for the WEC Model.

## Disturbance estimation
A Kalman filter and an Extended Kalman filter are applied to estimate an unknown disturbance to a state space model.
A disturbance model is augmented to the original system, with the choice of a harmonic or persistant disturbance model.
A tool for tuning the process noise covariance matrix is also included in this work.

## Prediction
A recursive, autoregressive, time-series prediciton method is used in this work.


# For more theroretical background
See:
1. My thesis: to be published soon

2. My METS extended abstract: http://s36.a2zinc.net/clients/pennwell/nha2015/Public/Calendar.aspx?TrackId=1381,1406,1407&View=Calendar

3. My OMAE conference paper: 

