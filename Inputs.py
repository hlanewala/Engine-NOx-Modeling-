import math
import csv
import numpy

def TestCases(i,fuel,TestCase):
###############################################################################################################
#                  ICAO Engine Exhaust Emissions Data Bank - PW305                                #
##############################################################################################################
#                                  Single  Torodial                                                           #
###############################################################################################################      
    
    if TestCase[i]!=None:
        P3  = float(TestCase[i][2])*1000.0         #Pa
        T3  = float(TestCase[i][1])                #K
        #Total mass flow rate of air
        mdot_A = float(TestCase[i][5])             #(kg/s)
        #Total mass flow rate of fuel
        mdot_F = float(TestCase[i][6])             #(kg/s)
        #Emissions Data
        EICO_Actual=float(TestCase[i][9])          #(g/kgfuel)
        EINOx_Actual=float(TestCase[i][7])         #(g/kgfuel)
        EIUHC_Actual=float(TestCase[i][8])         #(g/kgfuel)
        Condition = TestCase[i][0]
        humidity = float(TestCase[i][4])
        Nozzle='Pressure'
    else:
        print 'Error - No More Cases'
        os.system('pause')
        os.system('exit')
        
    if fuel=='Propane':
        HHV_Propane=50368.0
        HHV_JetA=48020.0
        mdot_F=HHV_JetA/HHV_Propane*mdot_F
        C=3
        H=8
        Air=2.0*16.0+3.76*14.0*2.0
        FAR_Stoch=(C*12.0+1.0*H)/(5.0*Air)
    elif fuel=='Methane':
        HHV_Methane=55528.0
        HHV_JetA=48020.0
        mdot_F=HHV_JetA/HHV_Methane*mdot_F
        C=1
        H=4
        Air=2.0*16.0+3.76*14.0*2.0
        FAR_Stoch=(C*12.0+1.0*H)/(2.0*Air)
    elif fuel=='Jet':
        C=10.0
        H=22.0
        Air=2.0*16.0+3.76*14.0*2.0
        FAR_Stoch=(C*12.0+1.0*H)/(15.0*Air)
    else:
        print 'Error - Unknown fuel'
        os.system('pause')
    #Sends test cases back to original program    
    return [P3,T3,mdot_A,mdot_F,FAR_Stoch,Condition,Nozzle,humidity,EICO_Actual,EINOx_Actual,EIUHC_Actual]

def Max_Geometry():
###############################################################################################################
#                                                  Max Geometry
###############################################################################################################

    #General Combustor Geometry#######
    #Maximum Number of Reactors in each plane####
    ymax=5
    xmax=6
    zmax=1

    #Total number of nozzles in the entire combustor
    Number_of_Nozzles = 24.0
    
    return[zmax,xmax,ymax,Number_of_Nozzles]

def Geometry(R_Width,R_Length,R_Height,R_Volume,Reactors,Ratio_Width,Number_of_Nozzles,xmax,ymax,zmax):
###############################################################################################################
#                                             Reactor Volumes                                                 #
###############################################################################################################
    #Reactor Dimensions####
    #The Radius from the center line of the combustor to the center line of the engine
    Radius=(((8.7+5.5)/2.0)*0.0254)/2.0   #m
    #The total width from the middle of the nozzle to halfway inbetween the two nozzles
    ############
    #  #      ##
    #  #     # #
    #  #    #  #
    #  #   #   #
    #  #  #    #
    #  #  #    #
    #  #   #   #
    #  #    #  #
    #  #     # #
    #  #      ##
    ############
    Width =2.0*math.pi*Radius/Number_of_Nozzles     #m
    #After Calculating the width and the size of the fuel spray zone the offline reactor
    #would have a width of .1 inches compared to the 1.65 inches of the nozzle reactor.
    #Therefore we will be assuming that thtere is no offline reactor 
    #Widths##########################
    #In the line of the Nozzle
    R_Width[0] = 3.3*0.0254/2.0
    
    #Lengths#########################
    #PZ
    R_Length[0]=0.0
    R_Length[1]=0.74*0.0254
    #IZ
    R_Length[2]=0.95*0.0254
    #DZ
    R_Length[3]=0.685*0.0254
    R_Length[4]=0.89*0.0254
    R_Length[5]=0.8*0.0254
    #Heights#########################
    #All heights are specific to the given reactor
    R_Height[1,0]=(0.4/16.5*8.7)*0.0254
    R_Height[1,1]=(1.8/16.5*8.7)*0.0254
    R_Height[1,2]=(3.55/16.5*8.7)*0.0254
    R_Height[1,3]=(0.4/16.5*8.7)*0.0254

    R_Height[2,0]=(1.85/16.5*8.7)*0.0254
    R_Height[2,1]=(1.85/16.5*8.7)*0.0254
    R_Height[2,2]=(1.5/16.5*8.7)*0.0254

    R_Height[3,0]=(1.65/16.5*8.7)*0.0254
    R_Height[3,1]=(1.5/16.5*8.7)*0.0254
    R_Height[3,2]=(1.0/16.5*8.7)*0.0254

    R_Height[4,0]=(1.6/16.5*8.7)*0.0254
    R_Height[4,1]=(1.2/16.5*8.7)*0.0254
    R_Height[4,2]=(0.9/16.5*8.7)*0.0254

    R_Height[5,0]=(2.8/16.5*8.7)*0.0254

    

    #Volumes##########################
    R_Volume[0,0,0]=(math.pi*((3.3/2.0*0.0254)**2.0-(2.7/2.0*0.0254)**2.0))/2.0*0.34*R_Width[0]
    R_Volume[0,0,1]=(math.pi*((2.7/2.0*0.0254)**2.0))/2.0*0.34*R_Width[0]
    R_Volume[0,0,2]=(math.pi*((2.7/2.0*0.0254)**2.0))/2.0*0.66*R_Width[0]
    R_Volume[0,0,3]=(math.pi*((3.3/2.0*0.0254)**2.0-(2.7/2.0*0.0254)**2.0))/2.0*0.66*R_Width[0]

    for z in range(0,zmax):
        for x in range(1,xmax):
            for y in range(0,ymax):
                R_Volume[z,x,y]=R_Height[x,y]*R_Length[x]*R_Width[z]

    #Reactor Names####################
    Reactors[0,0,0] = 'PZ#1 - Near Wall Top'
    Reactors[0,0,1] = 'PZ#1 - Middle Top'
    Reactors[0,0,3] = 'PZ#1 - Near Wall Bottom'
    Reactors[0,0,2] = 'PZ#1 - Middle Bottom'

    Reactors[0,1,0] = 'PZ#2 - Near Wall Top'
    Reactors[0,1,1] = 'PZ#2 - Middle Top'
    Reactors[0,1,2] = 'PZ#2 - Middle Bottom'
    Reactors[0,1,3] = 'PZ#2 - Near Wall Bottom'

    Reactors[0,2,0] = 'IZ - Top'
    Reactors[0,2,1] = 'IZ - Middle'
    Reactors[0,2,2] = 'IZ - Bottom'

    Reactors[0,3,0] = 'DZ#1 - Top'
    Reactors[0,3,1] = 'DZ#1 - Middle'
    Reactors[0,3,2] = 'DZ#1 - Bottom'

    Reactors[0,4,0] = 'DZ#2 - Top'
    Reactors[0,4,1] = 'DZ#2 - Middle'
    Reactors[0,4,2] = 'DZ#2 - Bottom'

    Reactors[0,5,0] = 'Final Reactor'

    #Ratio of one reactor with to another is used to determine the amount of air
    #passing through inlet holes is comming into one side compared to the other
    return[R_Width,R_Length,R_Height,R_Volume,Reactors,Ratio_Width]

def FlowSplits(Air_Rate,Fuel_Rate,Flow_Splits,Ratio_Width,re):
###############################################################################################################
#                                                  Flow Splits                                                #
###############################################################################################################

    #Inlet Air Splits################################
    #Percentage Air Through Combustor Inlet Holes.
    #Each flow split corresponds to a percentage of flow passing through a jet,
    #louver or swirler. The percent flows came from P&WC
    NozzleFaceFlow = 3.0
    NozzleFaceFlow_1= 7.6
    NozzleFaceFlow_2= 3.3
    NozzleFaceFlow_3= 3.0
    Inlet1_1 = 4.3    
    Inlet1_2 = 2.4  
    Inlet2_1 = 4.3  
    Inlet2_2 = 3.4  
    Inlet3 = 4.2    
    Inlet4 = 10.3   
    Inlet5_1 = 5.8/2
    Inlet5_2 = 5.8/2  
    Inlet6 = 17.1
    Inlet7 = 8.6   
    Inlet8 = 2.8    
    Inlet9 = 1.9
    Inlet10 = 2.5   
    Inlet11 = 1.9     
    
    
    #Inlet air to each reactor########################
    
    Air_Rate [0,0,0] = 0.0
    Air_Rate [0,0,1] = NozzleFaceFlow+NozzleFaceFlow_1+NozzleFaceFlow_2*0.4
    Air_Rate [0,0,2] = NozzleFaceFlow_3+NozzleFaceFlow_2*0.6
    Air_Rate [0,0,3] = Inlet1_1+Inlet1_2

    Air_Rate [0,1,0] = Inlet2_1+Inlet2_2
    Air_Rate [0,1,1] = 0.0
    Air_Rate [0,1,2] = Inlet3*0.75
    Air_Rate [0,1,3] = Inlet3*0.25

    Air_Rate [0,2,0] = 0.0
    Air_Rate [0,2,1] = 0.0
    Air_Rate [0,2,2] = Inlet5_1

    Air_Rate [0,3,0] = (Inlet4+Inlet6)*0.3
    Air_Rate [0,3,1] = (Inlet4+Inlet6)*0.7+Inlet7*0.3
    Air_Rate [0,3,2] = Inlet7*0.7+Inlet5_2

    Air_Rate [0,4,0] = Inlet8
    Air_Rate [0,4,1] = 0.0
    Air_Rate [0,4,2] = Inlet9

    Air_Rate [0,5,0] = Inlet10+Inlet11
    
    #Inlet Fuel Splits######################
    #Set the amount of fuel entering each reactor
    Fuel_Rate[0,0,1] = 0.5
    Fuel_Rate[0,0,2] = 0.5
    
    #Product Flow#####################################
    #Parameters
    Parameter_1=0.15    #Primary Zone Near Wall reactor content (bottom of combustor)
    Parameter_2=0.25    #Primary Zone Near Wall reactor content (top of combustor)
    Parameter_3=0.2     #PZ Reactor Split
    Parameter_4=0.6     #IZ Near Wall Air Split
    Parameter_5=0.30    #Middle intermediate zone contents moving up towards wall  
    
    #Primary Zone #1 Reactors
    #Nozzle
    Flow_Splits[0,0,0,0]=[1.0,0,1,0]
    Flow_Splits[0,0,0,1]=[1.0,0,1,1]
    Flow_Splits[0,0,0,2]=[1.0,0,0,1]
    Flow_Splits[0,0,0,3]=[1.0,0,0,0]
    Flow_Splits[1,0,0,3]=[Parameter_1,0,0,2]
    
    #Primary Zone #2 Reactors
    #Nozzle
    Flow_Splits[0,0,1,0]=[Parameter_2,0,2,0]
    Flow_Splits[0,0,1,1]=[1.0-Parameter_2,0,2,0]
    Flow_Splits[0,0,1,2]=[1.0-Parameter_1,0,0,2]
    Flow_Splits[0,0,1,3]=[1.0,0,0,3]

    #Intermediate Zone Reactors
    Flow_Splits[0,0,2,0]=[re,0,2,1]
    Flow_Splits[0,0,2,1]=[(1.0-Parameter_4),0,2,2]
    Flow_Splits[1,0,2,1]=[1.0-Parameter_3,0,1,2]
    Flow_Splits[0,0,2,2]=[Parameter_3,0,1,2]
    Flow_Splits[1,0,2,2]=[1.0,0,1,3]

    #Dilution Zone Reactors #1
    Flow_Splits[0,0,3,0]=[(1.0-re)*(Parameter_5),0,2,1]
    Flow_Splits[0,0,3,1]=[(1.0-re)*(1.0-Parameter_5),0,2,1]
    Flow_Splits[0,0,3,2]=[Parameter_4,0,2,2]

    #Dilution Zone Reactors #2
    Flow_Splits[0,0,4,0]=[1.0,0,3,0]
    Flow_Splits[0,0,4,1]=[1.0,0,3,1]
    Flow_Splits[0,0,4,2]=[1.0,0,3,2]
    
    #Final Reactor

    Flow_Splits[0,0,5,0]=[1.0,0,4,0]
    Flow_Splits[1,0,5,0]=[1.0,0,4,1]
    Flow_Splits[2,0,5,0]=[1.0,0,4,2]

    return [Flow_Splits,Air_Rate,Fuel_Rate]

def Set_Reactor_Evaporation(Reactor_Evaporation,zmax,xmax,ymax):
###############################################################################################################
#                                  Reactors with Evaporation                                                  #
###############################################################################################################
    z=0
    x=0
    for y in range(1,3):
        Reactor_Evaporation[z][x][y]=1
    Reactor_Evaporation[0][1][2]=1
    return Reactor_Evaporation
