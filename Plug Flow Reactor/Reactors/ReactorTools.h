/////////////////////////////////////////////////////////////////////
///
/// \file Reactor.h
/// 
/// \author H. Lanewala
/// 
/// \brief This header file contains the functin prototypes and 
///        class declarations for the plug flow reactor utilities.
/// 
/// This header file contains the functin prototypes and class 
/// declarations for the plug flow reactor utilities.
///
/////////////////////////////////////////////////////////////////////

/// Cantera includes
#include <cantera/Cantera.h>
#include <cantera/zerodim.h>
#include <cantera/IdealGasMix.h>

/**
 * Determine the mixture composition corresponding to the
 * specified equivalence ratio.
 *
 * \param gas An ideal gas mixture object.
 * \param fuel The fuel species name.
 * \param phi The equivalence ratio.
 * \return x Array of species mole fractions.
 */
void Composition( const IdealGasMix &gas, 
		  const string &fuel, 
		  const double phi, 
		  double *x );

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
	  double *x );


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
	  double *x );

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
	  double *x );


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
	  double *x );
