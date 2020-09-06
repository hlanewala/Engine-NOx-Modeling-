/////////////////////////////////////////////////////////////////////
///
/// \file PlugFlow.cpp
/// 
/// \author H. Lanewala
/// 
/// \breif This file contains the main() function for the Plug Flow
///        reactor computations.
/// 
/// This file contains the main() function for the Plug Flow reactor
/// computations.
///
/////////////////////////////////////////////////////////////////////
#include "ReactorTools.h"

/**
 * This is the main function which calls the WSR and PFR utilities.
 */
int main() {
  
  //------------------------------------------------
  // Inputs
  //------------------------------------------------
  // inlet conditions
  double T = 1200.0;                    // Temperature [K]
  double p = 10.0*OneAtm;               // Pressure [Pa]
  string X = "H2:1.0, O2:0.5";          // initial molar composition

  // CVR parameters
  double vol = 1.0;      // volume [m^3]
  double time = 6.e-6;  // final time to integrate WSR to [s]


  //------------------------------------------------
  // Solve
  //------------------------------------------------
  try {

    //
    // create a gas object using grimech 3.0
    //
    IdealGasMix gas("h2o2.cti", "h2o2");
    int nsp = gas.nSpecies();
        
    // set the state of the uburned mixture
    gas.setState_TPX(T, p, X);


    //
    // Solve CVR
    //
    double* x = new double[nsp];
    CVR(gas, time, vol, T, p, x);
    gas.setState_TPX(T, p, x);


    // clean up memory
    delete[] x;

    // exit
    return 0;
  

  //------------------------------------------------
  // handle exceptions thrown by Cantera
  //------------------------------------------------
  } catch (CanteraError) {
    showErrors(cout);
    cout << " terminating... " << endl;
    appdelete();
    return 1;
  }

}
