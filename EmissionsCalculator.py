def Emissions(prod, mdot, fuel_mass, fuel):
    FuelIndex = fuel.massFractions()
    FuelIndex = FuelIndex.tolist()
    FuelIndex = FuelIndex.index(1)
    FuelMassFrac = prod.massFraction(FuelIndex)
    FuelMassFrac = FuelMassFrac*mdot*1000/fuel_mass
    if (prod.speciesName(0) == 'H2'):
        EINOx = prod.massFraction('NO')+prod.massFraction('NO2')
        EINOx = EINOx*(mdot)*1000/fuel_mass
        EICO = prod.massFraction('CO')
        EICO = EICO*(mdot)*1000/fuel_mass
    else:
        EINOx = prod.massFraction('no')+prod.massFraction('no2')
        EINOx = EINOx*(mdot)*1000/fuel_mass
        EICO = prod.massFraction('co')
        EICO = EICO*(mdot)*1000/fuel_mass

    return FuelMassFrac, EINOx, EICO
