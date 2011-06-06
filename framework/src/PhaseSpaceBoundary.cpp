/**
  @class PhaseSpaceBoundary

  A collection of constraints on observables, defining the phase space in which a data point exists

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "StringProcessing.h"
#include "PhaseSpaceBoundary.h"
#include "ObservableContinuousConstraint.h"
#include "ObservableDiscreteConstraint.h"
#include "ObservableRef.h"
//	System Headers
#include <iostream>
#include <stdlib.h>

//Default constructor
PhaseSpaceBoundary::PhaseSpaceBoundary() : allConstraints(), allNames()
{
}

//Constructor with correct arguments
PhaseSpaceBoundary::PhaseSpaceBoundary( vector<string> NewNames ) : allConstraints(), allNames(NewNames)
{
	allConstraints.reserve(NewNames.size());
	//Populate the map
	for (unsigned int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allConstraints.push_back(NULL);
	}
}

//Destructor
PhaseSpaceBoundary::~PhaseSpaceBoundary()
{
}

//Return the names of all bounds stored
vector<string> PhaseSpaceBoundary::GetAllNames()
{
	return allNames;
}

//Retrieve a constraint by its name
IConstraint * PhaseSpaceBoundary::GetConstraint( pair<string,int>* wanted_constraint )
{
	if( wanted_constraint->second != -1 )
	{
		return allConstraints[unsigned(wanted_constraint->second)];
	} else {
		wanted_constraint->second = StringProcessing::VectorContains( &allNames, &(wanted_constraint->first) );
	}
	if( wanted_constraint->second == -1 ){
		cerr << "PhysicsParameter " << wanted_constraint->first << " not found" <<endl;
	}else{
		return allConstraints[unsigned(wanted_constraint->second)];}
	exit(-1);
}

IConstraint * PhaseSpaceBoundary::GetConstraint( ObservableRef& object )
{
        if( object.GetIndex() < 0 ) {
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		if( object.GetIndex() >= 0 ) return allConstraints[ (unsigned) object.GetIndex() ];
	} else {
		return allConstraints[ (unsigned) object.GetIndex() ];
	}
	cerr << "Observable name " << object.Name().c_str() << " not found" << endl;
	throw(-20);
}

//Retrieve a constraint by its name
IConstraint * PhaseSpaceBoundary::GetConstraint(string Name)
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Constraint on " << Name << " not found" << endl;
		exit(1);
		//return new ObservableContinuousConstraint( Name, 0.0, 0.0, "NameNotFoundError" );
	}
	else
	{
		return allConstraints[unsigned(nameIndex)];
	}
}

//Set a constraint by name
bool PhaseSpaceBoundary::SetConstraint( string Name, IConstraint * NewConstraint )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Constraint on " << Name << " not found" << endl;
		exit(1);
		//return false;
	}
	else
	{
		//Delete the old constraint before overwriting pointer
		delete allConstraints[unsigned(nameIndex)];
		allConstraints[unsigned(nameIndex)] = NewConstraint;
		return true;
	}
}

//Initialise bound
bool PhaseSpaceBoundary::SetConstraint( string Name, double Minimum, double Maximum, string Unit )
{
	ObservableContinuousConstraint * newConstraint = new ObservableContinuousConstraint( Name, Minimum, Maximum, Unit );
	bool returnValue = SetConstraint( Name, newConstraint );
	return returnValue;
}
bool PhaseSpaceBoundary::SetConstraint( string Name, vector<double> Values, string Unit )
{
	ObservableDiscreteConstraint * newConstraint = new ObservableDiscreteConstraint( Name, Values, Unit );
	bool returnValue = SetConstraint( Name, newConstraint );
	return returnValue;
}

//Returns whether a point is within the boundary
bool PhaseSpaceBoundary::IsPointInBoundary( DataPoint * TestDataPoint )
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		//Check if test Observable exists in the DataPoint
		Observable * testObservable = TestDataPoint->GetObservable( allNames[nameIndex] );
		if ( testObservable->GetUnit() == "NameNotFoundError" )
		{
			cerr << "Observable \"" << allNames[nameIndex] << "\" expected but not found" << endl;
			return false;
		}
		else
		{
			//Check if the Observable fits
			if ( !allConstraints[nameIndex]->CheckObservable(testObservable) )
			{
				//cerr << "Observable \"" << allNames[nameIndex] << "\" value (" << testObservable->GetValue() << ") is outside boundary" << endl;
				return false;
			}
		}
	}

	//Point is within the boundary
	return true;
}

/*
//Returns true if the argument PhaseSpaceBoundary fits within this PhaseSpaceBoundary
bool PhaseSpaceBoundary::CheckBoundary( PhaseSpaceBoundary * TestPhaseSpaceBoundary )
{
for ( int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
{
//Check if the test SingleBound exists in both PhaseSpaceBoundaries
IConstraint * testConstraint = TestPhaseSpaceBoundary->GetConstraint( allNames[nameIndex] );
if (testConstraint->GetUnit() == "NameNotFoundError")
{
cerr << "Constraint \"" << allNames[nameIndex] << "\" expected but not found" << endl;
return false;
}
else
{
//Check if the Observable fits
if ( !allConstraints[nameIndex]->CheckConstraint(testConstraint) )
{
cerr << "Observable \"" << allNames[nameIndex] << "\" value is outside boundary" << endl;
return false;
}

//Check if the SingleBound fits
if ( nameIterator->second.GetMinimum() > testBound.GetMinimum() || testBound.GetMinimum() > nameIterator->second.GetMaximum() )
{
cerr << "Bound \"" << nameIterator->first << "\" minimum is outside the boundary" << endl;
return false;
}

if ( nameIterator->second.GetMaximum() < testBound.GetMaximum() || testBound.GetMaximum() < nameIterator->second.GetMinimum())
{
cerr << "Bound \"" << nameIterator->first << "\" maximum is outside the boundary" << endl;
return false;
}
}
}

//Test boundary fits
return true;
}*/
