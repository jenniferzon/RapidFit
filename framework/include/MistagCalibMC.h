
#ifndef _MISTAGCALIB_MC_H_
#define _MISTAGCALIB_MC_H_

#include "IMistagCalib.h"
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DataPoint.h"

#include <string>

class IMistagCalib;

using namespace::std;

class MistagCalibMC : public IMistagCalib
{
	public:
		MistagCalibMC( PDFConfigurator* configurator );
		~MistagCalibMC();

		void addParameters( vector<string>& parameterNames ) const;
		void setParameters( const ParameterSet& parameters );

		void addObservables( vector<string>& observableNames ) const;
		void setObservables( DataPoint* measurement );

		double D1() const;

		double D2() const;

		void Print() const;

		double q() const;

		double mistagBbar() const;

		double mistagB() const;

		bool eventIsTagged() const;

	protected:

		MistagCalibMC( const MistagCalibMC& );
		MistagCalibMC& operator=(const MistagCalibMC );

		double mistagOS() const;
		double mistagSS() const;
		double mistagOSOther() const;
		double mistagSSOther() const;

		double mistagSSB() const;
		double mistagOSB() const;
		double mistagSSBbar() const;
		double mistagOSBbar() const;

		bool OSTagged() const;

		bool SSTagged() const;


		double RealD2() const;
		double RealD1() const;

		int GetCombinedTag() const;

		bool _OSTagged, _SSTagged;

		int _tagOS, _tagSS, _combinedtag;
		double _mistagOS, _mistagSS;
		double _mistagP0_OS, _mistagP1_OS, _mistagSetPoint_OS, _mistagDeltaP1_OS, _mistagDeltaP0_OS, _mistagDeltaSetPoint_OS;
		double _mistagP0_SS, _mistagP1_SS, _mistagSetPoint_SS, _mistagDeltaP1_SS, _mistagDeltaP0_SS, _mistagDeltaSetPoint_SS;
        double _mistagP2_SS, _mistagDeltaP2_SS;

		double _storedD1, _storedD2;

		ObservableRef tagOSName, tagSSName, mistagOSName, mistagSSName;
		ObservableRef mistagP1Name_OS, mistagP0Name_OS, mistagSetPointName_OS, mistagDeltaP1Name_OS, mistagDeltaP0Name_OS, mistagDeltaSetPointName_OS;
		ObservableRef mistagP1Name_SS, mistagP0Name_SS, mistagSetPointName_SS, mistagDeltaP1Name_SS, mistagDeltaP0Name_SS, mistagDeltaSetPointName_SS;
        ObservableRef mistagP2Name_SS, mistagDeltaP2Name_SS;

		bool _debugMistag;

		bool _onTuple;
		bool _useFixedEta;

		bool _floatCalib;

		bool _untagged;
};

#endif

