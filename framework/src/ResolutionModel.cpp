/**
  @class ResolutionModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "ResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"

#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

RESMODEL_CREATOR( ResolutionModel );

//............................................
// Constructor 
ResolutionModel::ResolutionModel( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName("timeResolutionScale") ),
	eventResolutionName	( configurator->getName("eventResolution") ), isCacheValid(false),
    usePolynomialAcceptance(false)
{
        if( !quiet) cout << "ResolutionModel:: Instance created " << endl ;
        usePolynomialAcceptance = configurator->isTrue( "UsePolynomialAcceptance" ) ;
        if (usePolynomialAcceptance){
                if( configurator->getConfigurationValue( "PolynomialAcceptanceFile" ) != "" )
                {
			string fullFileName = StringProcessing::FindFileName( configurator->getConfigurationValue( "PolynomialAcceptanceFile" ), quiet );
                        if( !quiet ) cout << "ResolutionModel:: Constructing Polynomial Coefficients using file: " << fullFileName << endl ;
			int number_of_lines = 0;
			std::string line;
			std::ifstream myfile( fullFileName.c_str() );
			while (std::getline(myfile, line))	++number_of_lines;
			myfile.close();
			ifstream in;
			in.open(fullFileName.c_str());
			if( in.fail() )
			{
				cout << "PolynomialAcceptance::PolynomialAcceptance : failed to open acceptance file  '  " << fullFileName << "  '  " << endl;
				exit(1);
			}
			std::vector<double> knots;
			std::vector<std::vector<double>> coeffs;
			int ok = true;
			unsigned int this_line=0;
			while (ok)
			{
				if( this_line < number_of_lines - 7) knots.push_back(stream(in));
				else 
				{
					std::vector<double> coefflist;
					for(unsigned i=0 ; i<4 ; i++)	coefflist.push_back(stream(in));
					coeffs.push_back(coefflist);
				}
				++this_line;
				if( in.eof() == true || in.good() == false || this_line == number_of_lines ) ok = false;
			}
			in.close();
		        this->expPoly =  new Mathematics::ExpPolynomialAcceptance(knots, coeffs) ;
        		//expPoly->show() ;
                }
        }
}

//............................................
//Helpers	Only Used in Constructor
double ResolutionModel::stream(ifstream& thisStream)
{
	double tmpVal;
	thisStream >> tmpVal;
	return tmpVal;
}

//..........................
//This method allows the instance to add the parameters it needs to the list
void ResolutionModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( resScaleName );
	return ;
}

//..........................
//To take the current value of a parameter into the instance
void ResolutionModel::setParameters( ParameterSet & parameters )
{
	isCacheValid = ( resScale == parameters.GetPhysicsParameter( resScaleName )->GetValue() );
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	return ;
}

bool ResolutionModel::CacheValid() const
{
	return isCacheValid;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void ResolutionModel::addObservables( vector<string> & observableNames )
{
	observableNames.push_back( eventResolutionName );
	return ;
}

//..........................
//To take the current value of an obserable into the instance
void ResolutionModel::setObservables( DataPoint * measurement )
{
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return ;
}

//..........................
//To take the current value of an obserable into the instance
bool ResolutionModel::isPerEvent( ) {  return true ; }


//..............................
// Primitive Functions
double ResolutionModel::Exp( double time, double gamma ) {
    if(usePolynomialAcceptance) { return expPoly->Exp( time, gamma, eventResolution*resScale  ) ; }
	return Mathematics::Exp( time, gamma, eventResolution*resScale  ) ;
}
double ResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
    if(usePolynomialAcceptance) { return expPoly->ExpInt( tlow, thigh, gamma, eventResolution*resScale ) ;}
	return Mathematics::ExpInt( tlow, thigh, gamma, eventResolution*resScale ) ;
}
double ResolutionModel::ExpSin( double time, double gamma, double dms ) {
    if(usePolynomialAcceptance) {return expPoly->ExpSin( time, gamma, dms, eventResolution*resScale) ;}
	return Mathematics::ExpSin( time, gamma, dms, eventResolution*resScale) ;
}
double ResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
    if(usePolynomialAcceptance) {return expPoly->ExpSinInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;}
	return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}
double ResolutionModel::ExpCos( double time, double gamma, double dms ) {
    if(usePolynomialAcceptance) {return expPoly->ExpCos( time, gamma, dms, eventResolution*resScale) ;}
	return Mathematics::ExpCos( time, gamma, dms, eventResolution*resScale) ;
}
double ResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
    if(usePolynomialAcceptance) {return expPoly->ExpCosInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;}
	return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}

//..............................
// Wrappers of Primitive Functions for Robs clever stuff
double ResolutionModel::Exp_Wrapper( vector<double> input ) {
	return this->Exp( input[0], input[1]  ) ;
}

double ResolutionModel::ExpInt_Wrapper( vector<double> input ) {
	return this->ExpInt( input[0], input[1], input[2] ) ;
}

double ResolutionModel::ExpSin_Wrapper( vector<double> input ) {
	return this->ExpSin( input[0], input[1], input[2] ) ;
}
double ResolutionModel::ExpSinInt_Wrapper( vector<double> input ) {
	return this->ExpSinInt( input[0], input[1], input[2], input[3] ) ;
}

double ResolutionModel::ExpCos_Wrapper( vector<double> input ) {
	return this->ExpCos( input[0], input[1], input[2] ) ;
}
double ResolutionModel::ExpCosInt_Wrapper( vector<double> input ) {
	return this->ExpCosInt( input[0], input[1], input[2], input[3] )  ;
}



