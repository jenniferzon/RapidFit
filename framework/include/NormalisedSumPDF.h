/**
        @class NormalisedSumPDF

        An implementation of IPDF for adding the values of two other IPDFs, normalising them relative to each other.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-12
*/

#ifndef NORMALISED_SUM_PDF_H
#define NORMALISED_SUM_PDF_H

#include "IPDF.h"
#include "RapidFitIntegrator.h"

class NormalisedSumPDF : public IPDF
{
	public:
		NormalisedSumPDF();
		NormalisedSumPDF( IPDF*, IPDF*, PhaseSpaceBoundary* );
		NormalisedSumPDF( IPDF*, IPDF*, PhaseSpaceBoundary*, string );
		~NormalisedSumPDF();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters(ParameterSet*);

		//Return the integral of the function over the given boundary
		virtual double Integral(DataPoint*, PhaseSpaceBoundary*);

		//Return the function value at the given point
		virtual double Evaluate(DataPoint*);

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
                virtual vector<string> GetDoNotIntegrateList();

	private:
		void MakePrototypes( PhaseSpaceBoundary* );

		vector<string> prototypeDataPoint, prototypeParameterSet, doNotIntegrateList;
		IPDF * firstPDF;
	       	IPDF * secondPDF;
		RapidFitIntegrator * firstIntegrator;
		RapidFitIntegrator * secondIntegrator;
		double firstFraction, firstIntegralCorrection, secondIntegralCorrection;
		string fractionName;
		PhaseSpaceBoundary * integrationBoundary;
};

#endif