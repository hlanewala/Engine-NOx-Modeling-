import sys
sys.path.append('C:\Users\user\Documents\Master thesis\Cantera Simulation\PW305')
from Cantera import *
from Cantera.Reactor import *
from ReactorClass import *
from EmissionsCalculator import Emissions
from datetime import datetime
from numpy import numarray
import numpy
import math
import os
import csv

'''
This modeling procedure has following characterisitcs
1) This model consists of 10-19 PSRs in the radial direction in primary zone. This distribution of PSRs is based on an assumed normal distribution of equivalence ratio in
the PZ. Several values of unmixedness parameter were tried to observe the relationship between fuel/air unmixedness and the accuracy in prediction of emissions.
2) For all values of unmixedness, it was noticed that NOx and CO trended well across various power settings but the predictions lacked in accuracy.
3) More knowledge about this procedure can be obtained if higher values of unmixedness do not result in un-physical solutions which is currently a big issue.
Currently, the unmixedness value of greater than 0.28 results in some of the fuel and air burning at -ve equivalence ratios.
4) In this procedure, there is one primary zone, one intermediate zone and three dilution zones. Their air splits can be seen in the air flows section.
'''

def Procedure(P3, T3, phi_overall, mdot_A, mdot_F, y, NozzleType, *args):
    ##______________________________________________________
    ## Overall combustor data:
    ##------------------------------------------------------
    ## Pecentage Air flows through the liner
    ##------------------------------------------------------
    #Inlet Air Splits################################
    #Percentage Air Through Combustor Inlet Holes.
    #Each flow split corresponds to a percentage of flow passing through a jet,
    #louver or swirler. The percent flows came from P&WC
    NozzleFaceFlow_1 = 3.0/100
    NozzleFaceFlow_2= 3.3/100
    NozzleFaceFlow_3= 3.0/100
    NozzleFaceFlow_4= 7.6/100
    inlet5 = 11.3/100    
    inlet6 = 3.4/100   
    inlet7 = 11.4/100
    inlet8 = 29.0/100
    inlet9 = 3.0/100
    inlet10 = 2.4/100   
    inlet11 = 0.6/100    
    inlet12 = 1.8/100
    inlet13 = 4.0/100   
    inlet14 = 3.9/100
    inlet15 = 3.9/100
    inlet16 = 4.9/100
    inlet17 = 2.4/100
    inlet18 = 1.1/100
        
    ##-----------------------------------------------------
    ## Operating data
    ##-----------------------------------------------------

    ## mass flow in each segment (1/12 of the combustor)
    mdot_A = mdot_A/24
    mdot_F = mdot_F/24

    AirMass = 0
    Evaporation = args[13]
    SMD = args[14]
    Recirculation_Button = args[15]
    Parameter_3 = args[4]
    i_end = args[6]
    j_end = args[7]
    k_end = args[8]
    re = args[16]
    UnmixednessParameterValue = args[16]
    NumberOfPSRsInPZ = args[7]
    
    Init_Prod = numpy.empty([args[8],args[7],args[6]], dtype=object)
    Product = numpy.empty([args[8],args[7],args[6]], dtype=object)
    Leftoverfuel = numpy.empty([args[8],args[7],args[6]], dtype=object)
    a = numpy.empty([args[8],args[7],args[6]], dtype=object)
    mdot = numpy.zeros([args[8],args[7],args[6]])
    tres = numpy.zeros([args[8],args[7],args[6]])
    vol = numpy.zeros([args[8],args[7],args[6]])
    R_Height = numpy.zeros([args[6],args[7]])
    R_Length = numpy.zeros(args[6])
    #R_Width = numpy.zeros(args[8])
    SimulationRun = 0
    PZLocations = 0

    NReactors = 25
    Number_of_Nozzles = 24
    
    KineticMech = 'gri30.cti'
    fuel = importPhase(KineticMech)
    fuel.set(T = 300.0, P = P3, X = 'C3H8:1')
    #air = importPhase('air_capital.cti')
    air = importPhase(KineticMech)
    #air.set(T = T3, P = P3)
    air.set(T = T3, P = P3, Y = 'O2:0.23049799047, N2:0.74950874863, AR:0.013702863416, H2O:0.00629')

    ## stoichiometric fuel/air of JP4
    C, H = CH_Atoms(fuel)
    MolesOfAir=2.0*16.0+3.76*14.0*2.0
    f_stoich=(C*12.0+1.0*H)/(5.0*MolesOfAir)
    f_stoich = 0.063379

    ## total mass flow rate
    total_mdot = mdot_A + mdot_F

    ##-----------------------------------------------------    
    #               Reactor Volumes (m^3)
    ##-----------------------------------------------------
    # In python, the arrays are declared as: array([@ z = 0[@ y = 0[x],@ y = 1[x]], @ z = 1[@ y = 0[x],@ y = 1[x]]])
    # and they are accessed as: array[z,y,x] 
    # volumes = [@ z = 0[@ y = 0[volumes of reactors at a particular x-location], @y = 1[volumes of reactors at a particular x-location]],
    #            @ z = 1[@ y = 0[volumes of reactors at a particular x-location], @y = 1[volumes of reactors at a particular x-location]]]

    #Reactor Dimensions####
    #The Radius from the center line of the combustor to the center line of the engine
    Radius=(((8.7+5.5)/2.0)*0.0254)/2.0   #m
    #The total width from the middle of the nozzle to halfway inbetween the two nozzles
    ################################
    #                      #      ##
    #                      #     # #
    #                      #    #  #
    #                      #   #   #
    #                      #  #    #
    #                      #  #    #
    #                      #   #   #
    #                      #    #  #
    #                      #     # #
    #                      #      ##
    ################################
    Width =2.0*math.pi*Radius/Number_of_Nozzles     #m
    #Widths##########################
    #In the line of the Nozzle
    R_Width = 3.3*0.0254
    #Lengths#########################
    #R_Length[x]
    #IZ
    R_Length[1]=0.95*0.0254
    #DZ
    R_Length[2]=0.685*0.0254
    R_Length[3]=0.89*0.0254
    R_Length[4]=0.8*0.0254
    #Heights#########################
    #R_Height[x,y]
    #All heights are specific to the given reactor
  
    R_Height[1,0]=(1.85/16.5*8.7)*0.0254 + (1.85/16.5*8.7)*0.0254 + (1.5/16.5*8.7)*0.0254

    R_Height[2,0]=(1.65/16.5*8.7)*0.0254 + (1.5/16.5*8.7)*0.0254 + (1.0/16.5*8.7)*0.0254

    R_Height[3,0]=(1.6/16.5*8.7)*0.0254 + (1.2/16.5*8.7)*0.0254 + (0.9/16.5*8.7)*0.0254
    
    R_Height[4,0]=(2.8/16.5*8.7)*0.0254

    R_PZVolume = 0.00036102614196870998

    PSR = FlowData(UnmixednessParameterValue, NumberOfPSRsInPZ)
    Phis = UnmixednessData(y, UnmixednessParameterValue, NumberOfPSRsInPZ)

    for r in range(0, len(PSR)):
        vol[0,r,0] = (R_PZVolume/(len(PSR)))
    
    vol[0,0,1] = R_Length[1] * R_Height[1,0] * R_Width

    vol[0,0,2] = R_Length[2] * R_Height[2,0] * R_Width

    vol[0,0,3] = R_Length[3] * R_Height[3,0] * R_Width

    vol[0,0,4] = R_Length[4] * R_Height[4,0] * R_Width
    
    ##-----------------------------------------------------
    #                      Air Splits 
    ##-----------------------------------------------------
        
    airpercent = NozzleFaceFlow_1 + NozzleFaceFlow_2 + NozzleFaceFlow_3 + NozzleFaceFlow_4 + inlet5 + inlet17 + inlet16 + inlet15
    for r in range(0,len(PSR)):
        a[0,r,0] = 'airpercent*PSR[' + str(r) + ']'
    a[0,0,1] = 'inlet6+0.125*(inlet7+inlet8)'
    a[0,0,2] = '0.875*(inlet7+inlet8)+inlet14+inlet13'
    a[0,0,3] = 'inlet9+inlet12'
    a[0,0,4] = 'inlet10+inlet11+inlet18' 
    
    ##-----------------------------------------------------
    #                    Product Splits 
    ##-----------------------------------------------------

    Init_Prod[0,0,1] = 'FlowReactors('
    #Reactor 1 - 17
    for r in range(0,len(PSR)):
        Init_Prod[0,r,0] = 'FlowReactors(fuel, Phis[' + str(r) + ']*f_stoich*air_mdot, air, air_mdot, mechanism = KineticMech)'
        Init_Prod[0,0,1] += 'Product[0,' + str(r) + ',0], mdot[0,' + str(r) + ',0], '
    #Reactor 2
    Init_Prod[0,0,1] += 'air, air_mdot, mechanism = KineticMech)'
    
    #Reactor 3
    Init_Prod[0,0,2] = 'FlowReactors(Product[0,0,1], mdot[0,0,1], air, air_mdot, mechanism = KineticMech)'
    #Reactor 4
    Init_Prod[0,0,3] = 'FlowReactors(Product[0,0,2], mdot[0,0,2], air, air_mdot, mechanism = KineticMech)'
    #Reactor 5
    Init_Prod[0,0,4] = 'FlowReactors(Product[0,0,3], mdot[0,0,3], air, air_mdot, mechanism = KineticMech)'
    
    def PSRdetector(k,j,i):
        if (i == 0):
            return 1
        else:
            return 0
    ##-----------------------------------------------------
    #       No Input Required from this point on
    ##-----------------------------------------------------    

    # Primary, Intermediate and Dilution Zone Length
    Total_Length, PZ_Length, IZ_Length, DZ_Length = 0.1586, 0.0416, 0.0624, 0.0546     

    #Parameter_3 is inlet 2,3 and 4 air distribution parameter
    #Parameter_2 is inlet 10 (dilution) air distribution parameter
    #Parameter_1 is nozzle plane mixture parameter    
    while Parameter_3 <= args[5]:
        Parameter_2 = args[2]
        data = []
        while Parameter_2 <= args[3]:
            print('Current Unmixedness Parameter Value = ' + str(UnmixednessParameterValue))
            print ' '
            Parameter_1 = args[0]            
            data.append([Parameter_2, ' ', ' ', ' ', ' ', ' ', ' ', ' '])
            while Parameter_1 <= args[1]:
                data.append([Parameter_1, 'C3H8', 'EINOx', 'EICO', 'Phi', 'Effective Phi', 'Temperature', 'mdot', 'AirMass', 'tres', 'FuelMass', 'leftoverfuel', 'Reacted Fuel (%)'])
                ## Change the previous meaning of Parameter_3 to 0.3 in Reactor 1,2,3 and 4
                Counter = 0
                air_mdot = 0
                AirMass = 0                
                for i in range(0,i_end):
                    for k in range(0,k_end):
                        for j in range(0,j_end):                                                             
                            try:
                                #This ratio statement works only if there are same number of reactors in the plane of the nozzle and the plane between the
                                #fuel nozzles or if there are no reactors in the plane between the fuel nozzle
                                ratio = vol[k,j,i]/vol.sum(0)[j,i]
                                #The subsequent instructions are executed only if ratio is a positive finite value b/w 0 and 1 including the value 1
                                if(ratio == 0 or math.isnan(ratio)):
                                    break
                            except:
                                break
                            AirMass += air_mdot
                            try:
                                air_mdot = mdot_A*(eval(a[k,j,i]))
                            except TypeError:
                                air_mdot = 0
##                            print 'i = ', i, 'j = ', j, 'k = ', k
                            try:
                                prod = eval(Init_Prod[k,j,i])
                            except TypeError:
                                prod = 0
                            if PSRdetector(k,j,i):                                
                                if (Recirculation_Button == 'on' and i == 0):
                                    if k == 0:                                        
                                        stuff1 = Recirculation()
                                        Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot, leftoverfuel, RemainingSMD = stuff1[0], stuff1[1], stuff1[2], stuff1[3], stuff1[4], stuff1[5]                                        
                                        Recirculation_Button = 'on'
                                        s1 = 0
                                        try:
                                            SimulationRun = stuff1[6]
                                            SimulationRun = 1
                                        except:
                                            SimulationRun = 0
                                    elif k == 1:
                                        stuff2 = Recirculation()
                                        Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot, leftoverfuel, RemainingSMD = stuff2[0], stuff2[1], stuff2[2], stuff2[3], stuff2[4], stuff2[5]
                                        Recirculation_Button = 'off'
                                        s2 = 0
                                        try:
                                            SimulationRun = stuff2[6]
                                            SimulationRun = 1
                                        except:
                                            SimulationRun = 0
                                    
                                elif SimulationRun:                                    
                                    if k == 0:
                                        s1 += 6
                                        Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot, leftoverfuel, RemainingSMD = stuff1[s1], stuff1[s1+1], stuff1[s1+2], stuff1[s1+3], stuff1[s1+4], stuff1[s1+5]
                                        try:
                                            SimulationRun = stuff1[s1+6]
                                            SimulationRun = 1
                                            if (j == (j_end-1) and i == 0):
                                                Recirculation_Button = 'on'
                                        except:
                                            SimulationRun = 1
                                    elif k == 1:
                                        s2 += 6
                                        Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot, leftoverfuel, RemainingSMD = stuff2[s2], stuff2[s2+1], stuff2[s2+2], stuff2[s2+3], stuff2[s2+4], stuff2[s2+5]
                                        try:
                                            SimulationRun = stuff2[s2+6]
                                            SimulationRun = 1
                                        except:
                                            SimulationRun = 0
                                        
                                else:
                                    Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot, leftoverfuel, RemainingSMD = prod.PSR(vol[k,j,i], PZ_Length, Evaporation)

                                # Modification to contents of downstream reactors due to unreacted liquid fuel in the upstream PSR
                                if (leftoverfuel and i == 0):
                                    if k < 1 and Recirculation == 'off':
                                        Init_Prod[1,j,i] = ContentModifier(Init_Prod[1,j,i], leftoverfuel, k, j, i)
####                                        print 'Init_Prod[1,'+str(j)+',0] = ', Init_Prod[1,j,i]
                                    for Ydirection in range(0, j_end):
                                        try:
                                            if k == 0:
                                                leftoverfuel = leftoverfuel - leftoverfuel*Parameter_1
                                            Init_Prod[0,Ydirection,i+1] = ContentModifier(Init_Prod[0,Ydirection,i+1], leftoverfuel, k, j, i)
##                                            print 'Init_Prod[0,'+str(Ydirection)+',1] = ', Init_Prod[0,Ydirection,i+1]
                                        # This except statement applies to a situation where there is more upstream than downstream reactors in the
                                        # y-direction. For instance, in a model where there's two reactors in the PZ due to stoich reactor but only one
                                        # reactor in the IZ.
                                        except AttributeError:
                                            print ('Warning: Only Reactor[' + str(k) + ',' + str(Ydirection-1) + ',' + str(i+1) + '] ' + 'string was modified')
                                            pass
                                        
                                Counter += 1
                            else:                                
                                Product[k,j,i], mdot[k,j,i], tres[k,j,i], fuel_mdot = prod.PFR(vol[k,j,i], NReactors)
                                leftoverfuel = 0
                                Counter += 1
                                
                            Leftoverfuel[k,j,i] = leftoverfuel   
                            Reactedfuel = (fuel_mdot-leftoverfuel)/fuel_mdot*100
                            Phi = fuel_mdot-leftoverfuel
                            Air_mdot = mdot[k,j,i]-Phi
                            EffectivePhi = (Phi/(mdot[k,j,i]-Phi))/f_stoich 
                            Phi = (fuel_mdot/(mdot[k,j,i]-Phi))/f_stoich                            
                            print('Product['+str(k)+','+str(j)+','+str(i)+'] = ' + str(Product[k,j,i].temperature()))
                            Fuel, EINOx, EICO = Emissions(Product[k,j,i], mdot[k,j,i], fuel_mdot, fuel)
                            #temperature = (Product[k,j,i].temperature()-273)*1.8+32
                            temperature = Product[k,j,i].temperature()
##                            arguments = (EINOx, EICO, Phi, Product[k,j,i].temperature(), mdot[k,j,i], tres[k,j,i]*1000, Reactedfuel)
                            arguments = (EINOx, EICO, Phi, temperature, mdot[k,j,i], tres[k,j,i]*1000, Reactedfuel, EffectivePhi)
                            EINOx, EICO, Phi, Temperature, TotalMass, Tres, Reactedfuel, EffectivePhi = DigitRounder(*arguments)
                            ReactorString = ReactorIdentity(k,j,i,PZLocations)                                
                            data.append([ReactorString, Fuel, EINOx, EICO, Phi, EffectivePhi, Temperature, TotalMass, Air_mdot, Tres, fuel_mdot, leftoverfuel, Reactedfuel])
                        
                data.append([' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '])
                print(' ')    
                Parameter_1 += args[9]
                
            Parameter_2 += args[10]
            
        SMD = SMD/1E-6
        path = os.getcwd()
        path = PathString(path)
        os.chdir(path)
        s = 'PW305_'+str(y)+'_'+str(UnmixednessParameterValue)+'_'+str(NumberOfPSRsInPZ)+'_Analysis.csv'
        ##s = 'Fuel_Split_Analysis3.csv'
        csv_file = open(str(s),'wb')
        csv_writer = csv.writer(csv_file)
        for values in data:
            csv_writer.writerow(values)
        csv_file.close()
        
        Parameter_3 += args[11]
