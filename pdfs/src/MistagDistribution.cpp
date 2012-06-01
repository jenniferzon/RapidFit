// $Id: MistagDistribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class MistagDistribution MistagDistribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "MistagDistribution.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( MistagDistribution );

//Constructor
MistagDistribution::MistagDistribution(PDFConfigurator* configurator) : BasePDF(),
	// Physics parameters
	GFgammaName	( configurator->getName("MistagGamma" ) ),
	GFmuName	( configurator->getName("MistagMu" ) ),
	GFbetaName	( configurator->getName("MistagBeta" ) ),
	GFshoulderName	( configurator->getName("MistagShoulder" ) ),
	// Observables
	GFxName( configurator->getName("mistag") )
	, gamma(0.), mu(0.), beta(0.)
{
	std::cout << "Constructing PDF: MistagDistribution " << std::endl ;

	MakePrototypes();
}


//Make the data point and parameter set
void MistagDistribution::MakePrototypes()
{
	// Observables
	allObservables.push_back( GFxName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( GFgammaName );
	parameterNames.push_back( GFmuName );
	parameterNames.push_back( GFbetaName );
	parameterNames.push_back( GFshoulderName );

	allParameters = ParameterSet(parameterNames);
}

//Destructor
MistagDistribution::~MistagDistribution()
{
}


bool MistagDistribution::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  	gamma = allParameters.GetPhysicsParameter( GFgammaName )->GetValue();
  	mu    = allParameters.GetPhysicsParameter( GFmuName )->GetValue();
  	beta  = allParameters.GetPhysicsParameter( GFbetaName )->GetValue();
  	shoulder  = allParameters.GetPhysicsParameter( GFshoulderName )->GetValue();
	
	return isOK;
}



//Calculate the function value
double MistagDistribution::Evaluate(DataPoint * measurement)
{
	// Get the observable
	double x = 0.5 - measurement->GetObservable( GFxName )->GetValue();
	
	double returnValue ;
	
	if( x > mu ) returnValue = TMath::GammaDist( x, gamma, mu, beta );
	else if( shoulder > 0.0) returnValue = shoulder ;
	else returnValue = 0.0000001 ;

  	return returnValue ;
}


// Normalisation
double MistagDistribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	return 1.023  +  shoulder*mu;
}
