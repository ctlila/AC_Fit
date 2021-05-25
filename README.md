# AC_Fit
Software to fit AC magnetic susceptibility measurements

--------------------------------------------------------------------------------------------
!!! MATLAB'S OPTIMIZATION AND STATISTICS TOOLBOX NEEDED !!!
you can download it in matlab -> Home -> Add-Ons -> search for 'Optimization toolbox' and 'Statistics and Machine Learning Toolbox' -> Install
--------------------------------------------------------------------------------------------


Open .m file in matlab and either run directly (default settings) or change parameters of optional features.

The default values of all prompted variables are set as : 

    MW = 838808 (mg/mol)
    
    mass = 11.3 (mg) for sample
    
    mass_teflon = 7 (mg)
    
    nb_points = 21
    
    first_row = 1      #This is the row index of your first point. Generally 1.
    
    smallest_step_T =  0.5 (K)
    smallest_step_H =  100 (Oe)

You may change these default values in lines 58 to 117.

INPUT :
.dat file of AC measurements from PPMS ; Molecular weight of the compound, mass of the sample, mass of teflon.

OUTPUT : 
2 .png files for chi' and chi" fitted curves, and 2 .txt files with the fit parameters (chiT-chiS, tau, alpha, chiS) and with chi' and chi". 


The default values are suited to correspond to the data files provided in the repository :
FSS2_60D_ScanT_0.25T.dat
FSS2_60D_ScanT_0.1T.xlsx
FSS2_60D_ScanT_0T.dat
FSS2_60D_ScanH20K.xlsx

Consequently you may press enter for each prompted field if using one of these file.


Please report any error/bug.


Enjoy !
