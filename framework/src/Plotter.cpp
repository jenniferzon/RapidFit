/**
  @class Plotter

  A class for plotting PDF projections onto histograms

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-13
  */

#include "Plotter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include <iostream>
#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "StringProcessing.h"

using namespace std;

//Default constructor
Plotter::Plotter()
{
}

//Constructor with correct arguments
Plotter::Plotter( IPDF * NewPDF, IDataSet * NewDataSet ) : plotPDF(NewPDF), plotData(NewDataSet), pdfIntegrator( new RapidFitIntegrator(NewPDF) )
{
}

//Destructor
Plotter::~Plotter()
{
}

//Create a root file containing a projection plot over each observable
void Plotter::PlotAllObservables( string FileName )
{
	//Work out what to plot
	vector<double> combinationWeights;
	vector<string> combinationDescriptions;
	vector<DataPoint> allCombinations = GetDiscreteCombinations( combinationWeights, combinationDescriptions );
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	//Make the file
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );

	//Loop over all (integrable) observables
	vector<string> observableNames = plotPDF->GetPrototypeDataPoint();
	vector<string>::iterator observableIterator;
	for ( observableIterator = observableNames.begin(); observableIterator != observableNames.end(); observableIterator++ )
	{
		//Check the observable can be plotted
		bool continuous = !( plotData->GetBoundary()->GetConstraint( *observableIterator )->IsDiscrete() );
		bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &(*observableIterator) ) == -1 );
		if ( continuous && doIntegrate )
		{
			cout << "Projecting " << *observableIterator << endl;
			MakeObservablePlots( *observableIterator, allCombinations, combinationWeights, combinationDescriptions, rootFile );
		}
	}

	rootFile->Close();
}

//Create a root file containing a projection plot over one observable
void Plotter::PlotSingleObservable( string FileName, string ObservableName )
{
	//Work out what to plot
	vector<double> combinationWeights;
	vector<string> combinationDescriptions;
	vector<DataPoint> allCombinations = GetDiscreteCombinations( combinationWeights, combinationDescriptions );
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	//Check the observable can be plotted
	bool continuous = !( plotData->GetBoundary()->GetConstraint( ObservableName )->IsDiscrete() );
	bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &ObservableName ) == -1 );
	if ( continuous && doIntegrate )
	{
		//Make the file and do the plotting
		TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
		cout << "Projecting " << ObservableName << endl;
		MakeObservablePlots( ObservableName, allCombinations, combinationWeights, combinationDescriptions, rootFile );
		rootFile->Close();
	}
}

//Make a directory containing all the plots for a particular observable
//NOTE: this function assumes there is an open TFile already
void Plotter::MakeObservablePlots( string ObservableName, vector<DataPoint> AllCombinations, vector<double> CombinationWeights, vector<string> CombinationDescriptions, TFile * PlotFile )
{
	//////////////////////////////////////////////////////////////////
	//Make the histogram of this observable

	//Get the data needed to make the histogram
	double maximum, minimum;
	int binNumber;
	vector<double> observableValues = GetStatistics( ObservableName, minimum, maximum, binNumber );
	double binInterval = ( maximum - minimum ) / (double)binNumber;

	//Make the histogram
	string histogramName = ObservableName + "ProjectionPlot";
	string histogramTitle = ObservableName + " projection plot";
	TH1F * dataHistogram = new TH1F( histogramName.c_str(), histogramTitle.c_str(), binNumber, minimum, maximum );

	//Loop over all data points and add them to the histogram
	for ( int dataIndex = 0; dataIndex < observableValues.size(); dataIndex++)
	{
		dataHistogram->Fill( observableValues[dataIndex] );
	}

	dataHistogram->Sumw2();

	//Histogram complete
	////////////////////////////////////////////////////////////////////

	//Make and enter the directory
	string directoryName = ObservableName + "Directory";
	string directoryTitle = ObservableName + " directory";
	TDirectory * combinationDirectory = PlotFile->mkdir( directoryName.c_str(), directoryTitle.c_str() );
	combinationDirectory->cd();

	////////////////////////////////////////////////////////////////////
	//Calculate the projection for each combination, and the data averaged projection

	vector<double> dataAverageProjection;
	double averageIntegral = 0.0;
	double plotInterval;
	
	int plotNumber;
	if (ObservableName == "time" ) plotNumber = 128;
	else plotNumber = 32;

	//Initialise the data averaged projection
	for ( int pointIndex = 0; pointIndex < plotNumber; pointIndex++ )
	{
		dataAverageProjection.push_back(0.0);
	}

	//Loop over all discrete combinations
	for ( int combinationIndex = 0; combinationIndex < AllCombinations.size(); combinationIndex++ )
	{
		//Calculate the projection for this combination
		vector<double> projectionValueVector = ProjectObservable( AllCombinations[combinationIndex], ObservableName, minimum, maximum, plotNumber, plotInterval );
		double projectionIntegral = pdfIntegrator->Integral( &( AllCombinations[combinationIndex] ), plotData->GetBoundary() );

		//Update the data average values, and make the projection graph arrays
		double projectionValueArray[plotNumber];
		double observableValueArray[plotNumber];
		for ( int pointIndex = 0; pointIndex < plotNumber; pointIndex++ )
		{
			//Projection graph
			projectionValueArray[pointIndex] = projectionValueVector[pointIndex] * binInterval * observableValues.size() / projectionIntegral;
			observableValueArray[pointIndex] = minimum + ( plotInterval * pointIndex );

			//Data average
			dataAverageProjection[pointIndex] += projectionValueVector[pointIndex] * CombinationWeights[combinationIndex];
		}
		averageIntegral += projectionIntegral * CombinationWeights[combinationIndex];

		//Plot the projection
		MakePlotCanvas( ObservableName, CombinationDescriptions[combinationIndex], dataHistogram, observableValueArray, projectionValueArray, plotNumber );
	}

	//Projecting complete
	/////////////////////////////////////////////////////////////////////

	//Leave the directory
	combinationDirectory->GetMotherDir()->cd();

	//Plot the data averaged projection
	double projectionValueArray[plotNumber];
	double observableValueArray[plotNumber];
	for ( int pointIndex = 0; pointIndex < plotNumber; pointIndex++ )
	{
		projectionValueArray[pointIndex] = dataAverageProjection[pointIndex] * binInterval * observableValues.size() / averageIntegral;
		observableValueArray[pointIndex] = minimum + ( plotInterval * pointIndex );
	}
	MakePlotCanvas( ObservableName, "", dataHistogram, observableValueArray, projectionValueArray, plotNumber );
}

//Create a projection plot for a given observable
//NOTE: this function assumes there is an open TFile already
void Plotter::MakePlotCanvas( string ObservableName, string Description, TH1F * Histogram, double * ProjectionXValues, double * ProjectionYValues, int PlotNumber )
{
	//Make the graph
	TMultiGraph * graph = new TMultiGraph();
	TGraphErrors * projectionGraph = new TGraphErrors( PlotNumber, ProjectionXValues, ProjectionYValues );
	graph->Add(projectionGraph);
	
	//Formatting
	EdStyle * greigFormat = new EdStyle();	
	projectionGraph->SetLineColor(2);
	projectionGraph->SetLineWidth(4);
	string xTitle = ObservableName + " (" + plotData->GetDataPoint(0)->GetObservable(ObservableName)->GetUnit() + ")";
	Histogram->GetXaxis()->SetTitle( xTitle.c_str() );
	
	string drawOptions = "C";
	//Note "same" option is not required if graph is drawn second. "A" option - requesting axes - will overwrite whatever was there before

	//Stick both objects on one canvas
	string canvasName = ObservableName + "Projection" + Description;
	string canvasTitle = ObservableName + " projection " + Description;
	TCanvas * bothPlots = new TCanvas( canvasName.c_str(), canvasTitle.c_str() );
 	gStyle->SetMarkerStyle(0);
	Histogram->SetStats(0);

	double ymin = Histogram->GetBinContent(Histogram->GetMinimumBin());
	double ymax = Histogram->GetBinContent(Histogram->GetMaximumBin());
	
	if (ObservableName == "time" ){
		// Attempted hack to get round the fact that ROOT is rubbish!
		//TH1F * tmpHist = new TH1F("tmp", "tmp", 1, -2., -1.);
		//tmpHist->Fill(-2.);
		//tmpHist->Draw();
		Histogram->SetStats(1);
		bothPlots->SetLogy();
		Histogram->GetXaxis()->SetRangeUser(-2., 15.); // This does not work because of STUPID ROOT!!
		Histogram->GetYaxis()->SetRangeUser(0.1, 2*ymax);
		//projectionGraph->SetMinimum(-2.); // Don't think this works either	
	}
	else
	{
		Histogram->GetYaxis()->SetRangeUser(0., ymax + ymax/10.);
	}
	
	Histogram->Draw();
	graph->Draw( drawOptions.c_str() );
	
	bothPlots->Update();
	bothPlots->Write();
}

//Return the values tracing the PDF projection
vector<double> Plotter::ProjectObservable( DataPoint InputPoint, string ObservableName )
{
	double maximum, minimum, plotInterval;
	int plotNumber;
	GetStatistics( ObservableName, minimum, maximum, plotNumber );
	return ProjectObservable( InputPoint, ObservableName, minimum, maximum, 16, plotInterval );
}

//Return the values tracing the PDF projection
vector<double> Plotter::ProjectObservable( DataPoint InputPoint, string ObservableName, double Minimum, double Maximum, int PlotNumber, double & PlotInterval )
{
	PlotInterval = ( Maximum - Minimum ) / (double)( PlotNumber - 1 );

	//Find the value of the observable projection at each data point
	vector<double> pointValues;
	for ( int pointIndex = 0; pointIndex < PlotNumber; pointIndex++ )
	{
		//cout << "Integration " << pointIndex + 1 << " of " << PlotNumber << endl;
		double observableValue = Minimum + ( PlotInterval * pointIndex );
		InputPoint.GetObservable(ObservableName)->SetValue(observableValue);
		pointValues.push_back( pdfIntegrator->ProjectObservable( &InputPoint, plotData->GetBoundary(), ObservableName ) );
	}
	return pointValues;
}

//Run the statistics functions
vector<double> Plotter::GetStatistics( string ObservableName, double & Minimum, double & Maximum, int & BinNumber )
{
	//Make a vector of the values found for the observable
	vector<double> observableValues;
	for ( int dataIndex = 0; dataIndex < plotData->GetDataNumber(); dataIndex++ )
	{
		observableValues.push_back( plotData->GetDataPoint(dataIndex)->GetObservable(ObservableName)->GetValue() );
	}

	//Find out number and range of data points needed
	Maximum = StatisticsFunctions::Maximum(observableValues);
	Minimum = StatisticsFunctions::Minimum(observableValues);
	BinNumber = StatisticsFunctions::OptimumBinNumber(observableValues);
	return observableValues;
}

//Return a list of data points
//Each should take the data average value of each continuous observable
//Each should represent one combination of possible discrete values
vector<DataPoint> Plotter::GetDiscreteCombinations( vector<double> & DataPointWeights, vector<string> & DataPointDescriptions )
{
	//Calculate all possiblle combinations of discrete observables
	vector< vector<double> > discreteValues;
	vector<string> discreteNames, continuousNames;
	vector<string> allNames = plotPDF->GetPrototypeDataPoint();
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, plotData->GetBoundary(), discreteNames, continuousNames, discreteValues );

	//Initialise the data averaging
	vector<double> continuousSums;
	vector<long> combinationCounts;
	for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
	{
		continuousSums.push_back(0.0);
	}
	for ( int combinationIndex = 0; combinationIndex < discreteCombinations.size(); combinationIndex++ )
	{
		combinationCounts.push_back(0);
	}

	//Examine the data set. Find the average value for each continuous observable, and the weight for each discrete combination
	for ( int dataIndex = 0; dataIndex < plotData->GetDataNumber(); dataIndex++ )
	{
		DataPoint * readDataPoint = plotData->GetDataPoint(dataIndex);

		//Sum the continuous values, in preparation for taking the average
		for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
		{
			continuousSums[continuousIndex] += readDataPoint->GetObservable( continuousNames[continuousIndex] )->GetValue();
		}

		//Calculate the index for the discrete combination, and increment the corresponding count
		int combinationIndex = 0;
		int incrementValue = 1;
		for ( int discreteIndex = discreteNames.size() - 1; discreteIndex >= 0; discreteIndex-- )
		{
			double currentValue = readDataPoint->GetObservable( discreteNames[discreteIndex] )->GetValue();

			for ( int valueIndex = 0; valueIndex < discreteValues[discreteIndex].size(); valueIndex++ )
			{
				if ( discreteValues[discreteIndex][valueIndex] == currentValue )
				{
					combinationIndex += ( incrementValue * valueIndex );
					incrementValue *= discreteValues[discreteIndex].size();
					break;
				}
			}
		}
		combinationCounts[combinationIndex]++;
	}

	//Calculate averages and weights
	vector<double> combinationWeights;
	double dataNumber = (double)plotData->GetDataNumber();
	for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
	{
		continuousSums[continuousIndex] /= dataNumber;
	}

	for ( int combinationIndex = 0; combinationIndex < discreteCombinations.size(); combinationIndex++ )
	{
		combinationWeights.push_back( (double)combinationCounts[combinationIndex] / dataNumber );
	}

	//Create the data points to return
	vector<DataPoint> newDataPoints;
	vector<string> allDescriptions;
	DataPoint templateDataPoint = *( plotData->GetDataPoint(0) );
	for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
	{
		Observable * newValue = templateDataPoint.GetObservable( continuousNames[continuousIndex] );
		newValue->SetValue( continuousSums[continuousIndex] );
		templateDataPoint.SetObservable( continuousNames[continuousIndex], newValue );
	}
	for ( int combinationIndex = 0; combinationIndex < discreteCombinations.size(); combinationIndex++ )
	{
		string description = "(";

		//Output the discrete values for this combination
		for ( int discreteIndex = 0; discreteIndex < discreteNames.size(); discreteIndex++ )
		{
			//Set the data point
			Observable * newValue = templateDataPoint.GetObservable( discreteNames[discreteIndex] );
			newValue->SetValue( discreteCombinations[combinationIndex][discreteIndex] );
			templateDataPoint.SetObservable( discreteNames[discreteIndex], newValue );

			//Make the description
			char value[100];
			sprintf( value, "%f", discreteCombinations[combinationIndex][discreteIndex] );
			string addToDescription;
			if ( discreteIndex == discreteNames.size() - 1 )
			{
				addToDescription = discreteNames[discreteIndex] + "=" + value;
			}
			else
			{
				addToDescription = discreteNames[discreteIndex] + "=" + value + "_";
			}
			description.append(addToDescription);
		}

		description.append(")");
		allDescriptions.push_back(description);
		newDataPoints.push_back(templateDataPoint);
	}

	//Output the results
	DataPointDescriptions = allDescriptions;
	DataPointWeights = combinationWeights;
	return newDataPoints;
}