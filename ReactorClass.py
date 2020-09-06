from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from datetime import datetime
import string
import numpy


class FlowReactors:

    def __init__(self, *args, **kwargs):

        self._code = len(args)/2
                
        if (self._code < 1):
            raise AttributeError('Please give atleast one gas as an input')                         

        self._tupledargs = args
        a = self.ArgumentArray(args)
        
        self._argslength = len(args)
        self._mechanism = kwargs['mechanism']
        try:
            self._SMD = kwargs['SMD']
        except:
            self._SMD = 0
        try:
            self._Tpz = kwargs['Tpz']
            self._tres = kwargs['tres']
        except:
            self._Tpz = 0
            self._tres = 0
        #self._Tpz = 0
        self._gas1 = a[0]
        self._args = a
        self._PureAirMass = self._args[len(a)-1]
        self._RemainingSMD = 0
        self._UnreactedFuel = 0
        self._EvapCounter = 0
        self._f = 0
           
    def ArgumentArray(self, args):
        a = numpy.empty(len(args),dtype=object)
        for i in range(0,len(args)):
            a[i] = args[i]

        return a

    def InitialGas(self):                    
        N2 = self.CharAnalyzer(self._gas1)
        s = str(N2)+str(':1')
        initial_gas = importPhase(self._mechanism)
        initial_gas.set(T = 2000, P = self._gas1.pressure(), X = s)
        return (initial_gas)     

    def CharAnalyzer(self, gas):
        if (gas.speciesName(0) == 'H2'):
            N2 = 'N2'
        else:
            N2 = 'n2'

        return N2   

    def FuelMassAnalyzer(self, gas, mass_flow):
        FuelSpeciesContainer = gas.speciesNames()
        MW_Carbon = gas.atomicWeights('C')[0]
        MW_Hydrogen = gas.atomicWeights('H')[0]
        f = 0

        LoopCounter = 0    
        for item in FuelSpeciesContainer:
            C_Atom = gas.nAtoms(item, 'C') ## It has been checked that 'C' or 'H' is valid for mechanisms of both upper and lower case types
            H_Atom = gas.nAtoms(item, 'H')
            MW_Species = gas.molecularWeights()
            MW_Species = MW_Species[LoopCounter]
            Y_species = gas.massFraction(LoopCounter)
            f += Y_species*(C_Atom*MW_Carbon+H_Atom*MW_Hydrogen)/MW_Species
            LoopCounter += 1

        f = mass_flow*f
        return f            

    def LiquidFuelDetector(self):
        self._args = self.ArgumentArray(self._tupledargs)
        Counter = 0
        LiquidFuelMass = numpy.zeros(self._argslength/2)
        self._LiquidFuelMassLocation = []
        for code in range(0, self._argslength, 2):
            gas = self._args[code]
            Fuel_Species = gas.massFractions().tolist()
            try:
                Fuel_Species = Fuel_Species.index(1.0)
                LiquidFuelMass[Counter] = float(self._args[code+1])
                self._LiquidFuelMassLocation.append(code+1)
                Counter += 1
            except ValueError:
                pass
        LiquidFuelMass = LiquidFuelMass.sum(0)
        return LiquidFuelMass

    def ReactingFuelMass(self):
        Locations = len(self._LiquidFuelMassLocation)
        for i in range(0, Locations):
            if i == 0:
                self._args[self._LiquidFuelMassLocation[i]] = self._f
            else:
                self._args[self._LiquidFuelMassLocation[i]] = 0                
        
        
    def Evaporation(self, gas, volume, length, tres): # Michael's Evaporation Model
        
        d2_law = 'Y'
        SMD = self._SMD
        Tpz = self._Tpz
        P3 = gas.pressure()
        #T3 = gas.temperature()
        ma = self._airmass
        #ma = self._PureAirMass #to be consistent with Mike
        mf = self.LiquidFuelDetector()
        A = volume/length        
        f = mf/ma
        
        T3_fuel = 298.0
        T_boil_fuel = 300.0+273.15
        if Tpz < T_boil_fuel:
            Tpz = Tpz+1000
        T = (Tpz + T_boil_fuel)/2.0
        
        #Jet-A
        #Fuel Properties####################################
        #All values are taken from 
        Cp_fuel                     = 1968.0                                            #(J/(kg*K))
        Thermal_conductivity_fuel   = -0.0002*T3_fuel+0.1667                            #(W/(m*K))
        rho_fuel                    = 797.88                                            #(kg/m^3)
        heat_of_vaporization_fuel   = -1.6218*T3_fuel**2.0+705.12*T3_fuel+320000.0      #(J/kg)
        kinetic_viscosity_fuel      = 1.65*10**-6                                       #(m^2/s)
        dynamic_viscosity_fuel      = rho_fuel * kinetic_viscosity_fuel                 #(N*s/m^2)
        surface_tension_fuel        = 25.35/1000.0                                      #(N/m)
        FAR_stoch                   = 0.067
        T_boil_fuel                 = 300.0+273.15        

        #Air Properties#####################################
        air = Air()
        air.set(T = (Tpz + T_boil_fuel)/2.0, P = P3)
        rho_air                     = air.density()
        Thermal_conductivity_air    = air.thermalConductivity()             #(W/(m*K))
        Cp_air                      = air.cp_mass()                         #(J/(kg*K))
        dynamic_viscosity_air       = air.viscosity()                       #(N*s/m^2)
        velocity_air                = ma/rho_air/A                          #(m/s) 
        velocity_fuel               = mf/rho_fuel/A 
        
        #Evaporation Constants##############################
        #Mass Transfer Unit
        print 'Tpz = ', Tpz 
        B   = (Tpz - T_boil_fuel)/heat_of_vaporization_fuel*Cp_fuel
        if B<=0.0:
            B=0.0
        #print 'Mass Transfer = '+str(B)
        #Reynolds Number
        Re  = SMD*abs(velocity_air-velocity_fuel)*rho_air/dynamic_viscosity_air
        #Re  = 0.0
        #print 'Reynolds Number = '+str(Re) 
        #Prandtl Number
        Pr  = Cp_air*dynamic_viscosity_air/Thermal_conductivity_air
        #print 'Prandtl Number = '+str(Pr) 
        #Thermal Conductivity - Ref: Introduction to Combustion pg. 377
        Thermal_conductivity = 0.4*Thermal_conductivity_fuel + 0.6*Thermal_conductivity_air
        #print 'Thermal Conductivity = '+str(Thermal_conductivity)
        #Number of Droplets - Gas Turbine Combustion pg. 54
        n=6.0/math.pi*rho_air/rho_fuel*volume/SMD**3.0*f
        
        #Liquid --> Gas Flow Rate###########################
        #Evaporation Constant
        Evaporation_Constant_Priniples  = 8.0*(Thermal_conductivity/Cp_fuel/rho_fuel)*math.log1p(B)*(1.0+0.30*Re**(0.5)*Pr**(0.33))
        #print 'Evaporation_Constant_Priniples = ', Evaporation_Constant_Priniples
##        print 'tres =', tres
        #Spray Evaporation - Gas Turbine Combustion pg. 54
        if d2_law=='N':
            gas  = rho_fuel*Evaporation_Constant_Priniples*SMD*math.pi/6.0*n
            print gas
            if gas>=mf:
                SMD_new=0.0
                gas=mf
                mf_new=0.0
            elif gas<=0.0:
                SMD_new=SMD
                gas=0.0
                mf_new=mf
            else:
                #Check what values this equation gives
                SMD_new = (SMD**3.0-6.0/1.0*gas/math.pi/rho_fuel/n*t)**(1/3)
                mf_new=mf-gas
                
        #D2 Law - Introduction to Combustion pg 377
        else:
            try:
                self._RemainingSMD = math.sqrt(SMD**2-Evaporation_Constant_Priniples*tres)
            except:
                self._RemainingSMD = 0
                
            if self._RemainingSMD>SMD:
                raise ValueError('Remaining SMD is calculated to be greater than the original SMD')
            else:
                #self._f is the reacted portion of the liquid fuel mass (mf)
                self._f = mf*(1-((self._RemainingSMD/2.0)**3.0/(SMD/2.0)**3.0))  
                self._UnreactedFuel = mf-self._f
                self.ReactingFuelMass()
##                print 'Remaining SMD = ', self._RemainingSMD

    def HasnainEvaporation(self, gas, volume, length, tres):
        # Variables for Calculation
        SMD = self._SMD
        Tpz = self._Tpz
        P3 = gas.pressure()        
        ma_total = self._airmass #Total Air in the reactor
        mf_total = self._fuelmass #Total Fuel in the reactor
        ma = self._PureAirMass #Air at T3 coming directly from the compressor in to this reactor
        mf = self.LiquidFuelDetector() #Fuel at 300 K coming directly from the fuel pump as well as unevaporated fuel coming from those reactors whose outputs 
                                       #are inputs to this reactor
        A = volume/length        
        f = mf/ma
        
        T3_fuel = 298.0
        T_boil_fuel = 300.0+273.15
        if Tpz < T_boil_fuel:
            Tpz = Tpz+1000
        T = (Tpz + T_boil_fuel)/2.0
        
        #Jet-A
        #Fuel Properties####################################
        #All values are taken from
        if T<=700:
            Cp_fuel                 = 1968.0#2.407178E2+5.09965*T-6.29026E-4*T**2-1.07155E-6*T**3                                            #(J/(kg*K))
        else:
            Cp_fuel                 = 1968.0#-1.35345889E4+9.14879E1*T-2.207E-1*T**2+2.91406E-4*T**3-2.153074E-7*T**4+8.386E-11*T**5-1.34404E-14*T**6 #(J/(kg*K))    
        Thermal_conductivity_fuel   = -0.0002*T_boil_fuel+0.1667                                  #(W/(m*K))
        rho_fuel                    = 797.88                                            #(kg/m^3)
        heat_of_vaporization_fuel   = -1.6218*T_boil_fuel**2.0+705.12*T_boil_fuel+320000.0                  #(J/kg)
        kinetic_viscosity_fuel      = 1.65*10**-6                                       #(m^2/s)
        dynamic_viscosity_fuel      = rho_fuel * kinetic_viscosity_fuel                 #(N*s/m^2)
        surface_tension_fuel        = 25.35/1000.0                                      #(N/m)      

        #Air Properties#####################################
        air = Air()        
        air.set(T = T, P = P3)
        #The reason that the volume of air is determined by the density of the input air below is that the input air is at T3 while the air created above is at
        #an average temperature T as can be seen. Volume of air is not really needed for this evaporation model anyways.
        volume_air                  = ma_total*tres/self._args[self._argslength-2].density()
        rho_air                     = air.density()
        Thermal_conductivity_air    = air.thermalConductivity()             #(W/(m*K))
        Cp_air                      = air.cp_mass()                         #(J/(kg*K))
        dynamic_viscosity_air       = air.viscosity()                       #(N*s/m^2)
        #For velocity of air, ma_total is used because total reacted fuel amount is negligible compared to the total air in a reactor
        velocity_air                = ma_total/rho_air/A                    #(m/s)
        #For velocity of fuel, mf is used because the evaporation is only applied to liquid portion of the fuel
        velocity_fuel               = mf/rho_fuel/A 
        
        #Evaporation Constants##############################
        #Mass Transfer Unit
        print 'Tpz = ', Tpz
##        print 'T_boil_fuel = ', T_boil_fuel
##        print 'heat_of_vaporization = ', heat_of_vaporization_fuel
##        print 'Cp_fuel = ', Cp_fuel
        B   = (Tpz - T_boil_fuel)/heat_of_vaporization_fuel*Cp_fuel
##        print 'Mass Transfer = '+str(B)
        
        #Reynolds Number
        Re  = SMD*abs(velocity_air-velocity_fuel)*rho_air/dynamic_viscosity_air
        #print 'Reynolds Number = '+str(Re) 
        #Prandtl Number
        Pr  = Cp_air*dynamic_viscosity_air/Thermal_conductivity_air
        #print 'Prandtl Number = '+str(Pr) 
        #Thermal Conductivity - Ref: Introduction to Combustion pg. 377
        Thermal_conductivity = 0.4*Thermal_conductivity_fuel + 0.6*Thermal_conductivity_air
        #print 'Thermal Conductivity = '+str(Thermal_conductivity)
        #Number of Droplets - Gas Turbine Combustion pg. 54
        #n = 6.0/math.pi*rho_air/rho_fuel*volume_air/SMD**3.0*f
        
        #Liquid --> Gas Flow Rate###########################
        #Evaporation Constant
        Evaporation_Constant_Priniples  = 8.0*(Thermal_conductivity/Cp_fuel/rho_fuel)*math.log1p(B)*(1.0+0.30*Re**(0.5)*Pr**(0.33))

        # Evaporation time calculation
        te = SMD**2/Evaporation_Constant_Priniples
##        print 't = ', tres
##        print 'te = ', te
        # Correlation for effective fuel-air ratio
        f = f*(-0.0295*te/tres+0.45)
        self._f = f*ma
        self._UnreactedFuel = mf-self._f
        self.ReactingFuelMass()
        print 'Correlation modified Phi = ', f/0.06378
####        print 'f = ', f

        if f < 0 and self._EvapCounter > 0:
            raise ValueError('Calculated Effective Fuel/Air is Negative')
        self._EvapCounter += 1

        try:
            self._RemainingSMD = math.sqrt(SMD**2-Evaporation_Constant_Priniples*tres)
        except:
            self._RemainingSMD = 0

    def PSRCalc(self, volume, tfinal):
        initial_gas = self.InitialGas()

        gas = self._gas1
        NoGas = 0 ## Ignition isn't provided if NoGas = 0
        mass_flow = 0

        upstreams = numpy.empty(self._code, dtype=object)
        ms = numpy.empty(self._code, dtype=object)

        reactor = Reactor(initial_gas, volume=volume, energy='on')

        ControllerCount = 0
        for code in range(0, self._argslength, 2):
            upstreams[ControllerCount] = Reservoir(self._args[code])
            ms[ControllerCount] = MassFlowController()
            ms[ControllerCount].install(upstreams[ControllerCount], reactor)
            ms[ControllerCount].set(self._args[code+1])
            mass_flow += self._args[code+1]
            ControllerCount += 1
                   
        exhaust = Reservoir(gas)
        
        v = Valve()
        v.install(reactor, exhaust)
        v.setValveCoeff(Kv=0.5)#Change made from 1.0 to 0.5
        
        sim = ReactorNet([reactor])
        
        tnow = 0.0
        
        tracker = datetime.now()
        LoopCounter = 0
        while (tnow < tfinal):
            LoopCounter += 1                     
            tnow = sim.step(tfinal)
            tres = reactor.mass()/mass_flow 
            currenttime = datetime.now()
            d = reactor.massFractions()

            IndexCounter = 0
            for item in d:
                if item > 1:
                    badguy = 1                    
                    baditem = item
                    break
                else:
                    badguy = 0
                IndexCounter += 1
            if badguy:
                break
            
            b = (currenttime.time().minute - tracker.time().minute)
            if(b > 2):
                break
         
        if(IndexCounter >= gas.nSpecies()):
            badSpecie = 'No Bad Species present'
            baditem = 'None'
        else:
            badSpecie = gas.speciesName(IndexCounter)
            
        tres = reactor.mass()/v.massFlowRate()
        T = reactor.temperature()
        P = reactor.pressure()
        reactor = Reactor(initial_gas)
        x = reactor.contents().moleFractions()
        initial_gas.setState_TPX(T,P,x)

        return initial_gas, mass_flow, tres


    def PFR(self, volume, NReactors):
        initial_gas = self.InitialGas()
        gas = self._gas1
        T = gas.temperature()
        P = gas.pressure()
        x = gas.moleFractions()
        initial_gas.setState_TPX(T,P,x)

        upstreams = numpy.empty(self._code-1, dtype=object)
        ms = numpy.empty(self._code-1, dtype=object)
        
        mass_flow = self._args[1]
        
        TOL = 1.0E-10
        Niter = 20

        nsp = gas.nSpecies()
  
        wdot = ['']*nsp
        wold = ['']*nsp

        volume_n = volume / NReactors

        tres = 0.0
               
        for i in range(0,NReactors):
            reactor = Reactor(initial_gas, volume=volume_n, energy='on')

            upstream = Reservoir(initial_gas)
            downstream = Reservoir(initial_gas)

            m = MassFlowController()
            m.install(upstream, reactor)
            m.set(mass_flow)

            if (i == 0):
                ControllerCount = 0
                for code in range(2, self._argslength, 2):
                    upstreams[ControllerCount] = Reservoir(self._args[code])
                    ms[ControllerCount] = MassFlowController()
                    ms[ControllerCount].install(upstreams[ControllerCount], reactor)
                    ms[ControllerCount].set(self._args[code+1])
                    mass_flow += self._args[code+1]
                    ControllerCount += 1
            
            v = Valve()
            v.install(reactor, downstream)
            v.setValveCoeff(Kv=0.1)

            sim = ReactorNet([reactor])
            
            dt = reactor.mass()/mass_flow

            tnow = 0.0
            wold = initial_gas.netProductionRates()
            
            while(tnow < Niter*dt):
                tnow += dt
                sim.advance(tnow)

                max_change = 0.0
                wdot = initial_gas.netProductionRates()

                for k in range(0,nsp):
                    max_change = max(math.fabs(wdot[k]-wold[k]), max_change)
                    wold[k] = wdot[k]

                if (max_change < TOL):
                    break
                
            tres +=reactor.mass()/mass_flow

            T = reactor.temperature()
            P = reactor.pressure()
            reactor = Reactor(initial_gas)
            x = reactor.contents().moleFractions()
            initial_gas.setState_TPX(T,P,x)

        f = self.FuelMassAnalyzer(initial_gas, mass_flow)
            
        return initial_gas, mass_flow, tres, f    

    def PSR(self, volume, length, Evaporation = 'off', tfinal = 100.0):
        
        initial_gas, mass_flow, tres = self.PSRCalc(volume, tfinal)
        self._fuelmass = self.FuelMassAnalyzer(initial_gas, mass_flow)
        self._airmass = mass_flow-self._fuelmass
        #The code between Evaporation process start and end only has significance if evaporation is 'on'
        #Evaporation process start
        if self._Tpz == 0 and Evaporation == 'on':
            T = initial_gas.temperature()
        elif Evaporation == 'on':
            T = 0
            tres = self._tres
        else:
            pass
        
        if (Evaporation == 'on' and self._SMD != 0):
            LoopCounter =  0
            while (abs(self._Tpz-T)>1):
                if LoopCounter == 0 and T != 0:
                    self._Tpz = T
                elif LoopCounter == 0 and T == 0:
                    T = self._Tpz
                    
                self._Tpz = self._Tpz + 1.0*(T-self._Tpz)
                self.Evaporation(initial_gas, volume, length, tres)
                initial_gas, mass_flow, tres = self.PSRCalc(volume, tfinal)
                T = initial_gas.temperature()
                LoopCounter += 1
                if LoopCounter > 14:
                    self._args[self._LiquidFuelMassLocation[0]] = self._f + self._UnreactedFuel
                    initial_gas, mass_flow, tres = self.PSRCalc(volume, tfinal)
                    print ('Warning: Evaporation Model did not converge')
                    self._UnreactedFuel = 0
                    break
        #Evaporation process end
        else:
            self._UnreactedFuel = 0
      
        f = self._fuelmass
     
        leftoverfuel = self._UnreactedFuel
        
        if (leftoverfuel < 1E-10):
            leftoverfuel = 0

        SMD = self._RemainingSMD
        
        #The mass flow returned to the main program does not include the mass of the left over fuel. 

        return initial_gas, mass_flow, tres, f, leftoverfuel, SMD

def DigitRounder(*args):
    EINOx = round(args[0],2)
    EICO = round(args[1],2)
    Phi = round(args[2],3)
    Temperature = int(round(args[3],0))
    MassFlow = round(args[4], 10)
    Tres = round(args[5],1)
    Reactedfuel = round(args[6],1)
    EffectivePhi = round(args[7],3)
    return EINOx, EICO, Phi, Temperature, MassFlow, Tres, Reactedfuel, EffectivePhi

def ContentModifier(s, lof, k, j, i):
    str_sep = 'mdot['+str(k)+','+str(j)+','+str(i)+']'
    Req_String = string.split(s, str_sep)
    if len(Req_String) == 1:
        return s
    else:
        x = Req_String[0].rsplit('],',1)[1]
        x = [x, str(lof)]
        x = string.join(x,'')
        
        s1 = [Req_String[0],str_sep]
        s1 = string.join(s1,'')
        s2 = ', fuel,'+x
        s2 = [s2,Req_String[1]]
        s2 = string.join(s2,'')

        s = [s1,s2]
        s = string.join(s,'')
        return s

def ReactorIdentity(k,j,i,x):
    if (i == 0 or i == x):
        ReactorString = 'PZ_'
    elif (i == x+1):
        ReactorString = 'IZ_'
    else:
        ReactorString = 'DZ_'

    if k == 0:
        ReactorString = ReactorString + 'nozzle_'
    else:
        ReactorString = ReactorString + 'between_'

    ReactorString = ReactorString + str(j)
    return ReactorString

def CH_Atoms(fuel):
    FuelIndex = fuel.massFractions()
    FuelIndex = FuelIndex.tolist()
    FuelIndex = FuelIndex.index(1)
    C_Atoms = fuel.nAtoms(FuelIndex, 'C')
    H_Atoms = fuel.nAtoms(FuelIndex, 'H')
    if (C_Atoms == 0 and H_Atoms == 0):
        C_Atoms = fuel.nAtoms(FuelIndex, 'c')
        H_Atoms = fuel.nAtoms(FuelIndex, 'h')
    return C_Atoms, H_Atoms

def PathString(s):
    d = string.split(s, '\\')
    try:
        d.index('Parametric Study Tests')
        return s
    except ValueError:
        s = s+'\Parametric Study Tests'
        return s

def SMDCalculator(*args):
    length = len(args)
    AverageSMD = 0
    Total = 0
    for i in range(0, length, 2):
        AverageSMD += float(args[i]*args[i+1])
        Total += args[i+1]

    if Total == 0:
        return AverageSMD
    else:
        AverageSMD = AverageSMD/Total
        return AverageSMD        

def UnmixednessData(condition, UnmixednessParameterValue, NReactors = ''):
    s = os.getcwd()
    try:
        s = s+'\PhisandFlows_'+str(UnmixednessParameterValue)+'_'+str(NReactors)+'.csv'
    except:
        s = s+'\PhisandFlows_'+str(UnmixednessParameterValue)+'.csv'
        
    with open(str(s), 'r') as f:
        PhisandFlows = f.readlines()
        f.close

    if condition == 'Idle':
        condition = 7
        
    x = PhisandFlows[0]
    p = string.split(x, sep=',')

    try:
        indexer = p.index(str(condition))
    except ValueError:
        indexer = 4
    
    y = []
    for i in range(1,len(PhisandFlows)):
        x = PhisandFlows[i]
        p = string.split(x, sep=',')
        y.append(float(p[indexer]))
        
    return y

def FlowData(UnmixednessParameterValue, NReactors = ''):
    s = os.getcwd()
    try:
        s = s+'\PhisandFlows_'+str(UnmixednessParameterValue)+'_'+str(NReactors)+'.csv'
    except:
        s = s+'\PhisandFlows_'+str(UnmixednessParameterValue)+'.csv'
        
    with open(str(s), 'r') as f:
        PhisandFlows = f.readlines()
        f.close
        
    y = []   
    for i in range(1,len(PhisandFlows)):
        x = PhisandFlows[i]
        p = string.split(x, sep=',')
        y.append(float(p[0]))
    
    return y
