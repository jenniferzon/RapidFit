
//	RapidFit Headers
#include "FoamIntegrator.h"
#include "StatisticsFunctions.h"
//	System Headers
#include <math.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
FoamIntegrator::FoamIntegrator() : allIntegrators(), discreteNames(), discreteValues()
{
}

//Constructor with correct arguments
FoamIntegrator::FoamIntegrator( IPDF * InputPDF, IDataSet * InputData ) : allIntegrators(), discreteNames(), discreteValues()
{
	//Calculate all possible combinations of discrete observables
	vector<string> continuousNames;
	vector<string> allNames = InputPDF->GetPrototypeDataPoint();
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, InputData->GetBoundary(), discreteNames, continuousNames, discreteValues );

	//Create data points for each combination (using data averaged values for continuous, non-integrable functions)
	vector<string> dataPointDescriptions;
	vector<double> dataPointWeights;
	vector<DataPoint> combinationPoints =  StatisticsFunctions::DataAverage( InputData, discreteCombinations, discreteValues, discreteNames, continuousNames, dataPointDescriptions, dataPointWeights );

	//Make a foam for each discrete combination
	for (unsigned int combinationIndex = 0; combinationIndex < combinationPoints.size(); ++combinationIndex )
	{
		MakeFoam* combinationFoam = new MakeFoam( InputPDF, InputData->GetBoundary(), &( combinationPoints[combinationIndex] ) );
		allIntegrators.push_back(combinationFoam);
	}
}

//Destructor
FoamIntegrator::~FoamIntegrator()
{
}

//Select and run the correct integrator
double FoamIntegrator::Integral( DataPoint * InputPoint, PhaseSpaceBoundary * InputBoundary )
{
	//	Stupid gcc
	(void)InputBoundary;
	//The integral won't work if the boundary has changed, but you might want a check that it's the same

	//Use the data point to find the index of the correct foam
	int combinationIndex = 0;
	int incrementValue = 1;
	for ( int discreteIndex = int(discreteNames.size()) - 1; discreteIndex >= 0; --discreteIndex )
	{
		//Retrieve the observable value
		Observable * temporaryObservable = InputPoint->GetObservable( discreteNames[unsigned(discreteIndex)] );
		double currentValue = temporaryObservable->GetValue();

		//Calculate the index
		for (unsigned int valueIndex = 0; valueIndex < discreteValues[unsigned(discreteIndex)].size(); ++valueIndex )
		{
			if ( fabs(discreteValues[unsigned(discreteIndex)][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
			{
				combinationIndex += ( incrementValue * int(valueIndex) );
				incrementValue *= int(discreteValues[unsigned(discreteIndex)].size());
				break;
			}
		}
	}

	//Use the foam to integrate
	return allIntegrators[unsigned(combinationIndex)]->Integral();
}
