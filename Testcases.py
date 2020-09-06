from Annular_Combustor_Model import *
import os
##from checking import *

ProgramPath = os.getcwd()
timetracker = datetime.now()

# Pressure is in Pascals
# Temperature is in kelvin
# mdot_A is in kg/s

## args = (Parameter_1_start, Parameter_1_end, Parameter_2_start, Parameter_2_end, Parameter_3_start, Parameter_3_end, i_end, j_end, k_end,
##         Paremeter_1_increment, Paremeter_2_increment, Paremeter_3_increment, percentage of primary jet air going in to primary zone, Evaporation Switch, SMD
##         Recirculation Button)

HHV_Propane=50368.0
HHV_JetA=48020.0

##Evaporation = raw_input('Evaporation On or Off [y/n]? ')
##
##if Evaporation == 'y':
##    Evaporation = 'on'
##    SMD = raw_input('Enter SMD to start with at Idle = ')
##    SMDStopping = raw_input('Enter SMD to end with at Idle = ')
##
##    SMD30 = raw_input('Enter SMD to start with at 30% = ')
##    SMDStopping30 = raw_input('Enter SMD to end with at 30% = ')
##else:
##    Evaporation = 'off'
##    SMD = SMDStopping = SMD30 = SMDStopping30 = 0
##    
##Recirculation = raw_input('Recirculation On or Off [y/n]? ')
##if Recirculation == 'y':
##    Recirculation = 'on'
##    re = raw_input('Recirculation Percentage? ')
##    re = float(re)
##    re = re/100
##else:
##    Recirculation = 'off'
##    re = 0

Evaporation = 'off'
Recirculation = 'off'
ydir = 10

while ydir < 20:
    re = 0.3
    print '# of PSRs in Primary Zone = ', ydir

##    ##------------------------------------------------------
##    ##                        Idle
##    ##------------------------------------------------------
    while re <= 0.31:
        args2 = (0.05, 0.05, 0.05, 0.05, 0.07, 0.07, 5, ydir, 1,
            0.05, 0.05, 0.03, 0.25, Evaporation, 0, Recirculation,re)
        args3 = (0.05, 0.05, 0.05, 0.05, 0.07, 0.07, 5, ydir, 1,
            0.05, 0.05, 0.03, 0.25, Evaporation, 0, Recirculation,re)
        args1 = (0.05, 0.05, 0.05, 0.05, 0.07, 0.07, 5, ydir, 1,
                 0.05, 0.05, 0.03, 0.25, Evaporation, 0, Recirculation,re)
        print(str('Idle'))
        P3 = 315000
        #T3 = 435.0
        T3 = 430.0
        phi_overall = 0.19332
        mdot_A = 2.580
        #mdot_F = 0.04
        mdot_F = 0.034
        mdot_F=HHV_JetA/HHV_Propane*mdot_F
        Procedure(P3, T3, phi_overall, mdot_A, mdot_F, 7, 'Simplex', *args1)
        print ''
        os.chdir(ProgramPath)

    ##------------------------------------------------------
    ##                        30%
    ##------------------------------------------------------
        print(str(30))
        P3 = 760000
        #T3 = 615.0
        T3 = 612.0
        phi_overall = 0.21029
        mdot_A = 5.790
        mdot_F = 0.083
        mdot_F=HHV_JetA/HHV_Propane*mdot_F
        Procedure(P3, T3, phi_overall, mdot_A, mdot_F, 30, 'Simplex', *args1)
        print ''
        os.chdir(ProgramPath)
    ##------------------------------------------------------
    ##                    85% and 100%
    ##------------------------------------------------------
        print(str(85))
        P3 = 1600000
        #T3 = 736.0
        T3 = 734.0
        phi_overall = 0.31407
        mdot_A = 10.0
        mdot_F = 0.2141
        mdot_F=HHV_JetA/HHV_Propane*mdot_F
        Procedure(P3, T3, phi_overall, mdot_A, mdot_F, 85, 'Simplex', *args3)
        os.chdir(ProgramPath)

        print(str('TakeOff'))
        P3 = 1820000
        #T3 = 764.0
        T3 = 762.0        
        phi_overall = 0.34921
        mdot_A = 10.67
        mdot_F = 0.254
        mdot_F=HHV_JetA/HHV_Propane*mdot_F
        Procedure(P3, T3, phi_overall, mdot_A, mdot_F, 100, 'Simplex', *args3)
        os.chdir(ProgramPath)

        re += 0.05

    ydir += 1


print(datetime.now()-timetracker)
raw_input("Please press enter to exit the program")

