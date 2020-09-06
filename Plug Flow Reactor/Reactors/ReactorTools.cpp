/////////////////////////////////////////////////////////////////////
///
/// \file Reactor.h
/// 
/// \author H. Lanewala
/// 
/// \brief This file contains the function definitions and 
///        class members for the plug flow reactor utilities.
/// 
/// This header file contains the function definitions and class 
/// members for the plug flow reactor utilities.
///
/////////////////////////////////////////////////////////////////////
#include "ReactorTools.h"

/**
 * Determine the mixture composition corresponding to the
 * specified equivalence ratio.
 *
 * \param gas An ideal gas mixture object.
 * \param fuel The fuel species name.
 * \param phi The equivalence ratio.
 * \param x Array of species mole fractions.
 */
void Composition( const IdealGasMix &gas, 
		  const string &fuel, 
		  const double phi, 
		  double *x )
{

  // the number of species
  int nsp = gas.nSpecies();

  // first, compute the stoichiometric fuel air ratio
  double C_atoms=gas.nAtoms(gas.speciesIndex(fuel), gas.elementIndex("C"));
  double H_atoms=gas.nAtoms(gas.speciesIndex(fuel), gas.elementIndex("H"));
  double ax=C_atoms+H_atoms/4.0;
  double fa_stoic=1.0/(4.76*ax);
  
  
  // determine the inlet composition
  for(int k=0;k<nsp;k++){
    if(k==gas.speciesIndex(fuel)){      x[k]=1.0; }
    else if(k==gas.speciesIndex("O2")){ x[k]=0.21/phi/fa_stoic; }
    else if(k==gas.speciesIndex("N2")){ x[k]=0.79/phi/fa_stoic; }
    else{ x[k]=0.0; }
  }


}


/**
 * Integrate a Well Stirred Reactor from the given initial gas
 * state to the time 'tfinal'.  Return the gas contents.
 *
 * \param gas An ideal gas mixture object.
 * \param tfinal The final time in [s].
 * \param mass_flow The mass flow rate in the reactor [kg/s].
 * \param volume The volume of the reactor [m^3].
 * \return T The final temperature [K].
 * \return p The final pressure [Pa].
 * \return x Array of species mole fractions.
 */
void WSR( IdealGasMix &gas, 
	  const double tfinal,
	  const double mass_flow,
	  const double volume,
	  double &T,
	  double &p,
	  double *x )
{

  //------------------------------------------------
  // setup reactor network
  //------------------------------------------------
  // the number of species
  int nsp = gas.nSpecies();

  // create a reservoir for the fuel/air inlet
  Reservoir inlet;
  inlet.insert(gas);

  // create the WSR, and fill it in initially with unburned mixture
  Reactor reactor;
  reactor.insert(gas);
  reactor.setInitialVolume(volume);

  // create a reservoir for the exhaust. The initial composition
  // doesn't matter.
  Reservoir exhaust;
  exhaust.insert(gas);

  // create and install the mass flow controllers. Controllers
  // m1 and m2 provide constant mass flow rates, and m3 provides
  // a short Gaussian pulse only to ignite the mixture
  MassFlowController m1;
  m1.install(inlet, reactor);
  m1.setMassFlowRate(mass_flow);
  
  // to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
  // The steady-state behavior is independent of how we do this, so we'll
  // just use a stream of pure atomic hydrogen.
  gas.setState_TPX(T, p, "H:1.0");
  Reservoir igniter;
  igniter.insert(gas);

  // The igniter will use a Guassiam 'functor' object to specify the
  // time-dependent igniter mass flow rate.
  double A = 0.1;
  double FWHM = 0.2;
  double t0 = 1.0;
  Gaussian igniter_mdot(A, t0, FWHM);
  
  MassFlowController m2;
  m2.install(igniter, reactor);
  m2.setFunction(&igniter_mdot);

  // put a valve on the exhaust line to regulate the pressure
  Valve v;
  v.install(reactor, exhaust);
  double Kv = 1.0;
  v.setParameters(1, &Kv);
  
  // the simulation only contains one reactor
  ReactorNet sim;
  sim.addReactor(&reactor);


  //------------------------------------------------
  // solve the WSR to steady state
  //------------------------------------------------
  //
  // take single steps to steady state, writing the results to a CSV file
  // for later plotting.
  //

  // output file object
  ofstream fout("WSR.dat");
  
  // the header
  fout << "VARIABLES = \"t (s)\" \"T (K)\" \"P (Pa)\" \"tres (s)\"" << endl;
  for(int k=0;k<nsp;k++) fout << " \"x_" << gas.speciesName(k) << "\"";
  fout << endl;
  fout << "ZONE" << endl;

  //
  // integrate to specified time
  //
  double tnow = 0.0;
  double tres;
  while (tnow < tfinal) {

    // take a step
    tnow = sim.step(tfinal);
    tres = reactor.mass()/v.massFlowRate();

    // output
    fout << tnow << ", "
	 << reactor.temperature() << ", "
	 << reactor.pressure() << ", "
	 << tres << ", ";
    ThermoPhase& c = reactor.contents();
    for (int k = 0; k < nsp; k++)
      fout << c.moleFraction(k) << ", ";
    fout << endl;

  } // endwile
  
  // close file
  fout.close();


  // set the state of the burned mixture
  reactor.contents().getMoleFractions(x);
  T = reactor.temperature();
  p = reactor.pressure();

}

/**
 * Integrate a Constant Volume Reactor from the given initial gas
 * state to the time 'tfinal'.  Return the gas contents.
 *
 * \param gas An ideal gas mixture object.
 * \param tfinal The final time in [s].
 * \param volume The volume of the reactor [m^3].
 * \return T The final temperature [K].
 * \return p The final pressure [Pa].
 * \return x Array of species mole fractions.
 */
void CVR( IdealGasMix &gas, 
	  const double tfinal,
	  const double volume,
	  double &T,
	  double &p,
	  double *x )
{

  //------------------------------------------------
  // setup reactor network
  //------------------------------------------------
  // the number of species
  int nsp = gas.nSpecies();

  // create the CVR, and fill it in initially with unburned mixture
  Reactor reactor;
  reactor.insert(gas);
  reactor.setInitialVolume(volume);
  
  // the simulation only contains one reactor
  ReactorNet sim;
  sim.addReactor(&reactor);


  //------------------------------------------------
  // solve the CVR to steady state
  //------------------------------------------------
  //
  // take single steps to steady state, writing the results to a CSV file
  // for later plotting.
  //

  // output file object
  ofstream fout("CVR.dat");
  
  // the header
  fout << "VARIABLES = \"t (s)\" \"T (K)\" \"P (Pa)\"" << endl;
  for(int k=0;k<nsp;k++) fout << " \"x_" << gas.speciesName(k) << "\"";
  fout << endl;
  fout << "ZONE" << endl;

  //
  // integrate to specified time
  //
  double tnow = 0.0;
  while (tnow < tfinal) {

    // take a step
    tnow = sim.step(tfinal);

    // output
    fout << tnow << ", "
	 << reactor.temperature() << ", "
	 << reactor.pressure() << ", ";
    ThermoPhase& c = reactor.contents();
    for (int k = 0; k < nsp; k++)
      fout << c.moleFraction(k) << ", ";
    fout << endl;

  } // endwile
  
  // close file
  fout.close();


  // set the state of the burned mixture
  reactor.contents().getMoleFractions(x);
  T = reactor.temperature();
  p = reactor.pressure();

}


/**
 * Integrate a Plug Flow Reactor from the given initial gas
 * state to steady state.  Return the gas contents.  The plug flow
 * reactor is represented by a linear chain of zero-dimensional 
 * reactors. The gas at the inlet to the first one has the specified 
 * inlet composition, and for all others the inlet
 * composition is fixed at the composition of the reactor immediately
 * upstream. Since in a PFR model there is no diffusion, the upstream
 * reactors are not affected by any downstream reactors, and therefore
 * the problem may be solved by simply marching from the first to last
 * reactors, integrating each one to steady state.
 *
 * \param gas An ideal gas mixture object.
 * \param tfinal The final time in [s].
 * \param mass_flow The mass flow rate in the reactor [kg/s].
 * \param length The total length of the reactor [m].
 * \param area The cross-sectional area of the reactor [m^2].
 * \param NReactors The number of WSR reactors to use.
 * \param isothermal True if isothermal reactor.
 * \return T The final temperature [K].
 * \return p The final pressure [Pa].
 * \return x Array of species mole fractions.
 */
void PFR( IdealGasMix &gas,
	  const double mass_flow,
	  const double length,
	  const double area,
	  const int NReactors,
	  const bool isothermal,
	  double &T,
	  double &p,
	  double *x )
{

  // iteration parameters
  const double TOL = 1.E-10;
  const int Niter = 20;

  // the number of species
  int nsp = gas.nSpecies();

  // allocate production rate arrays
  double* wdot = new double[nsp]; // net
  double* wold = new double[nsp]; // net

  // compute the reactor element properties
  double volume_n = length * area / NReactors;
  double length_n = length / NReactors;

  // the output file
  ofstream fout("PFR.dat");

  // the header
  fout << "VARIABLES = \"tres (s)\" \"len (m)\" \"T (K)\" \"P (Pa)\"" << endl;
  for(int k=0;k<nsp;k++) fout << " \"x_" << gas.speciesName(k) << "\"";
  fout << endl;
  fout << "ZONE" << endl;

  //
  // loop over each reactor in the chain
  //
  double tres = 0.0;
  for ( int n=0; n<NReactors; n++ ) {

    //------------------------------------------------
    // setup reactor network
    //------------------------------------------------

    // create a new PFR reactor
    Reactor reactor;
    reactor.insert(gas);
    reactor.setInitialVolume(volume_n);
    reactor.setEnergy(!isothermal); // isothermal

    // create a reservoir to represent the reactor immediately
    // upstream. Note that the gas object is set already to the
    // state of the upstream reactor
    Reservoir upstream;
    upstream.insert(gas);

    // create a reservoir for the reactor to exhaust into. The
    // composition of this reservoir is irrelevant.
    Reservoir downstream;
    downstream.insert(gas);

    // We need a valve between the reactor and the downstream reservoir.
    // This will determine the pressure in the reactor. Set Kv large
    // enough that the pressure difference is very small.
    Valve v;
    v.install(reactor, downstream);
    double Kv = 1000.0;
    v.setParameters(1, &Kv);
    
    // The mass flow rate into the reactor will be fixed by using a
    // MassFlowController obbject.
    MassFlowController m;
    m.install(upstream, reactor);
    m.setMassFlowRate(mass_flow);
      
    // the simulation
    ReactorNet sim;
    sim.addReactor(&reactor);

    //set relative and absolute tolerances on the simulation
    sim.setTolerances( /*rtol*/1.0e-6, /*atol*/1.0e-15);


    //------------------------------------------------
    // solve PFR
    //------------------------------------------------
    // guess a timestep size using an estimate for the residance time
    double dt = reactor.mass()/mass_flow;

    // initialize
    double tnow = 0.0;
    gas.getNetProductionRates(wold);

    //
    // solve this PFR element to steady state
    //
    while ( tnow < Niter*dt ) {

      // step in time
      tnow += dt;
      sim.advance(tnow);

      // compute max change and store old rates
      double max_change = 0.0;
      gas.getNetProductionRates(wdot);
      for (int k = 0; k < nsp; k++) {
	max_change = max( fabs(wdot[k]-wold[k]), max_change);
	wold[k] = wdot[k];
      }

      // check convergence
      if ( max_change < TOL ) break;

    } // endwhile - t

    // compute the residance time
    tres += reactor.mass()/v.massFlowRate();
            
    // the new gas state
    T = reactor.temperature();
    p = reactor.pressure();
    reactor.contents().getMoleFractions(x);
    gas.setState_TPX(reactor.temperature(), reactor.pressure(), x);

    //
    // write the gas mole fractions and surface coverages
    // vs. distance
    //
    fout << tres << ", "
 	 << length_n*(n+1) << ", "
	 << reactor.temperature() << ", "
	 << reactor.pressure() << ", ";
    ThermoPhase& c = reactor.contents();
    for (int k = 0; k < nsp; k++)
      fout << c.moleFraction(k) << ", ";
    fout << endl;

  } // endfor - reactors

  // close file
  fout.close();

  // clean up memory
  delete[] wdot;
  delete[] wold;

}


/**
 * Integrate a Plug Flow Reactor from the given initial gas
 * state to steady state.  Return the gas contents.  The plug flow
 * reactor is modeled using Cantera's 'FlowReactor' object.     
 *
 * \param gas An ideal gas mixture object.
 * \param tfinal The final time in [s].
 * \param mass_flow The mass flow rate in the reactor per unit area [kg/m^2/s].
 * \param isothermal True if isothermal reactor.
 * \return T The final temperature [K].
 * \return p The final pressure [Pa].
 * \return x Array of species mole fractions.
 */
void PFR( IdealGasMix &gas,
	  const double tfinal,
	  const double mass_flow,
	  const bool isothermal,
	  double &T,
	  double &p,
	  double *x )
{

  // the number of species
  int nsp = gas.nSpecies();

  //------------------------------------------------
  // setup reactor / simulation
  //------------------------------------------------
  // create a new PFR reactor
  FlowReactor reactor;
  reactor.insert(gas);
  reactor.setMassFlowRate(mass_flow);
  reactor.setEnergy(!isothermal); // isothermal

  // the simulation
  ReactorNet sim;
  sim.addReactor(&reactor);

  //------------------------------------------------
  // step through the PFR
  //------------------------------------------------

  // output file object
  ofstream fout("PFR.dat");
  
  // the header
  fout << "VARIABLES = \"tres (s)\" \"dist (m)\" \"vel (m/s)\" \"T (K)\" \"P (Pa)\"" << endl;
  for(int k=0;k<nsp;k++) fout << " \"x_" << gas.speciesName(k) << "\"";
  fout << endl;
  fout << "ZONE" << endl;

  //
  // integrate to specified time
  //
  double tres = 0.0;
  while (tres < tfinal) {

    // take a step
    tres = sim.step(tfinal);
    double dist = reactor.distance();
    double vel = reactor.speed();

    // output
    fout << tres << ", "
	 << dist << ", "
	 << vel << ", "
	 << reactor.temperature() << ", "
	 << reactor.pressure() << ", ";
    ThermoPhase& c = reactor.contents();
    for (int k = 0; k < nsp; k++)
      fout << c.moleFraction(k) << ", ";
    fout << endl;

  } // endwile
  
  // close file
  fout.close();


  // set the state of the burned mixture
  reactor.contents().getMoleFractions(x);
  T = reactor.temperature();
  p = reactor.pressure();
  
}


