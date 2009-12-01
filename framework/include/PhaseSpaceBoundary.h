/**
        @class PhaseSpaceBoundary

        A collection of constraints on observables, defining the phase space in which a data point exists

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PHASE_SPACE_BOUNDARY_H
#define PHASE_SPACE_BOUNDARY_H

#include <vector>
#include <string>
#include "IConstraint.h"
#include "DataPoint.h"

using namespace std;

class PhaseSpaceBoundary
{
	public:
		PhaseSpaceBoundary();
		PhaseSpaceBoundary( vector<string> );
		~PhaseSpaceBoundary();

		vector<string> GetAllNames();
		bool SetConstraint( string, IConstraint* );
		bool SetConstraint( string, double, double, string );
		bool SetConstraint( string, vector<double>, string );
		IConstraint * GetConstraint(string);
		bool IsPointInBoundary( DataPoint* );
		//bool CheckBoundary( PhaseSpaceBoundary* );

	private:
		vector< IConstraint* > allConstraints;
		vector<string> allNames;
};

#endif