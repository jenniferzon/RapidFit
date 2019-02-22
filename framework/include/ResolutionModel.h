/**
  @class ResolutionModel

  A class for implementing a decay time resolution model

  This one is a single Gaussian event-by-event model, with a single scale factor

  @author Pete Clarke
  @data 2013-06-03
  */

#pragma once
#ifndef Resolution_Model_H
#define Resolution_Model_H

//	RapidFir Headers
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "IResolutionModel.h"
#include "Observable.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

#include "Mathematics.h"

using namespace::std;

//=======================================

class ResolutionModel : public IResolutionModel
{
	public:

		ResolutionModel( PDFConfigurator* configurator, bool quiet=false ) ;

		void addParameters( vector<string> & parameterNames ) ;
		void setParameters( ParameterSet & parameters ) ;

		void addObservables( vector<string> & observableNames ) ;
		void setObservables( DataPoint * measurement ) ;

		double Exp( double time, double gamma ) ;
		double ExpInt( double tlow, double thigh, double gamma ) ;

		double ExpSin( double time, double gamma, double dms ) ;
		double ExpSinInt( double tlow, double thigh, double gamma, double dms ) ;

		double ExpCos( double time, double gamma, double dms ) ;
		double ExpCosInt( double tlow, double thigh, double gamma, double dms ) ;

		bool isPerEvent() ;

		//Wrappers
		double Exp_Wrapper( vector<double> input) ;
		double ExpInt_Wrapper( vector<double> input ) ;

		double ExpSin_Wrapper( vector<double> input ) ;
		double ExpSinInt_Wrapper( vector<double> input ) ;

		double ExpCos_Wrapper( vector<double> input ) ;
		double ExpCosInt_Wrapper( vector<double> input ) ;

		bool CacheValid() const;

	private:

		bool isCacheValid;
    
        bool usePolynomialAcceptance ;
        Mathematics::ExpPolynomialAcceptance * expPoly ;

		ObservableRef resScaleName;			// Scale to multiply e-by-e resolution
		double resScale ;

		ObservableRef eventResolutionName;  // Event-by-event resolution observable
		double eventResolution ;

		double stream(ifstream& stream);
		virtual unsigned int numComponents(){return 1;};
		virtual void requestComponent( unsigned int ){return;};

		virtual double GetFraction( unsigned int ){return 1.0;};

};

#endif

