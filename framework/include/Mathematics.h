/*!
 * @namespace Mathematics
 *
 * Namespace holding some common mathematical functions that are used
 * by many PDFs. For example, gaussian's convolved with trigonometric
 * functions.
 *
 * @author Greig A Cowan, Pete Clalrke greig.cowan@cern.ch
 */

#pragma once
#ifndef RAPIDFIT_MATHEMATICS
#define RAPIDFIT_MATHEMATICS

///	ROOT Headers
#include "TMath.h"
#include <cmath>
#include "RooMath.h"
#include "RooSentinel.h"
///	RapidFit Headers
#include "IDataSet.h"
#include "IPDF.h"
#include "RapidFitIntegrator.h"
#include "PhaseSpaceBoundary.h"
///	System Headers
#include <vector>

using namespace::std;

namespace Mathematics
{
	//	This comes up A LOT, save some CPU cycles and cache it!
	static const double _pi = atan(1.0)*4.;;
	static const double sqrt_2 = sqrt(2.);
	static const double _over_sqrt_2 = 1./sqrt_2;
	static const double _over_sqrt_2pi = 1./sqrt(2.*_pi);
	static const double invpi=TMath::InvPi();
	static const double global_frac = 0.28125 * invpi;
	static const double third= 1./3.;
	static const double twothird = 2.*third;
	static const double fourthird = 1.+third;
	static const double eightthird = 2.*fourthird;
	static const double root_6 = sqrt(6.);
	static const double root_3 = sqrt(3.);
	static const double root_2 = sqrt(2.);
	static const double rootpi= sqrt(_pi) ; //sqrt(atan2(0.,-1.));
    static const complex<double> rootminus1(0.,1.);

	inline double SQRT_2(){ return sqrt_2; }//=sqrt(2.);
	inline double _Over_SQRT_2(){ return _over_sqrt_2; }// = 1./SQRT_2;
	inline double _Over_SQRT_2PI(){ return _over_sqrt_2pi; }
	inline double _Over_PI(){ return invpi; }// = TMath::InvPi();//1./TMath::Pi();
	inline double Global_Frac(){ return global_frac; }// = 0.28125 * _Over_PI;
	inline double Root_6(){ return root_6; }
	inline double Root_3(){ return root_3; }
	inline double Root_2(){ return root_2; }
	inline double Third(){ return third; }
	inline double TwoThird(){ return twothird; }
	inline double FourThird(){ return fourthird; }
	inline double EightThird(){ return eightthird; }
	inline double Rootpi(){ return rootpi; }
	inline double Pi(){ return _pi; }

	complex<double> evalCerfApprox( double swt, double u, double c );
	complex<double> evalCerf( double swt, double u, double c );

	double evalCerfRe( double swt, double u, double c );

	double evalCerfIm( double swt, double u, double c );

	//--------------------------- exp and exp*sin and exp*cos time functions -------------------------
	// time functions for use in PDFs with resolution

	double Exp( double t, double gamma, double resolution );
	double Exp_Wrapper( vector<double> input );

	double ExpInt( double tlow, double thigh, double gamma, double resolution  );
	double ExpInt_Wrapper( vector<double> input );

	//Added these to include an upper time acceptance of form ( 1 - beta*t)
	double Exp_betaAcceptance( double t, double gamma, double resolution, double betaParameter );
	double Exp_betaAcceptance_Wrapper( vector<double> input );

	double ExpInt_betaAcceptance( double tlow, double thigh, double gamma, double resolution, double betaParameter );
	double ExpInt_betaAcceptance_Wrapper( vector<double> input );

	double ExpCosh( double t, double gamma, double deltaGamma, double resolution );
	double ExpCosh_Wrapper( vector<double> input );

	double ExpCoshInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	double ExpCoshInt_Wrapper( vector<double> input );

	double ExpSinh( double t, double gamma, double deltaGamma, double resolution );
	double ExpSinh_Wrapper( vector<double> input );

	double ExpSinhInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	double ExpSinhInt_Wrapper( vector<double> input );

	double ExpCos( double t, double gamma, double deltaM, double resolution );
	double ExpCos_Wrapper( vector<double> input );

	double ExpCosInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	double ExpCosInt_Wrapper( vector<double> input );

	double ExpSin( double t, double gamma, double deltaM, double resolution );
	double ExpSin_Wrapper( vector<double> input );

	double ExpSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	double ExpSinInt_Wrapper( vector<double> input );

	pair<double,double> ExpCosSin( double t, double gamma, double deltaM, double resolution );
	pair<double,double> ExpCosSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution );

	double expErfInt( double tlimit, double tau, double sigma);
	double expErfInt_Wrapper( vector<double> input );

	void getBs2JpsiPhiAngularFunctions( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6, const double cosTheta, const double phi, const double cosPsi);
	void getBs2JpsiPhiAngularFunctionsWithSwave( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6, double & f7, double & f8, double & f9, double & f10, const double cosTheta, const double phi, const double cosPsi);

	vector<double> calculateAcceptanceWeights( IDataSet * dataSet, IPDF * PDF );
    int calculateAcceptanceCoefficients( IDataSet * dataSet, IPDF * PDF );
    void calculateAcceptanceWeightsWithSwave( IDataSet * dataSet, IPDF * PDF );
    
    
    //New functions for spline code
    complex<double> splmom_kn( int n, complex<double> z ) ;
    double splmom_fact(int n) ;
    complex<double> splmom_erfc( complex<double> z ) ;
    complex<double> splmom_mn( int n, double x1, double x2, complex<double> z,
                              double theErf1, double theErf2,
                              complex<double> theErfc1, complex<double> theErfc2,
                              complex<double> theExp1, complex<double>theExp2 ) ;
    complex<double> splmom_integral_n( int n,  double res, double x1, double x2, complex<double> z ,
                                      double theErf1, double theErf2,
                                      complex<double> theErfc1, complex<double> theErfc2,
                                      complex<double> theExp1, complex<double>theExp2 ) ;

    //................................................................................
    // This was invented to deal with the introduction of 3rd order splines for the acceptance
    // The result is the acceptance being divided into a set od sectors
    // For each sector there are a set of 3rd order polynomial coefficients
    class ExpPolynomialAcceptance {
        
    private:
        const int nseg ;
        const int ncoeff  ;
        std::vector<double> knots ;
        std::vector<std::vector<double>> coeffs ;
        double value( double t ) ;
        complex<double> integral( double tlo, double thi, complex<double> q, double res) ;
        
    public:
        ExpPolynomialAcceptance() ;
        ExpPolynomialAcceptance( std::vector<double> knots, std::vector<std::vector<double>> coeffs ) ;
        void show() ;
        double Exp(double t, double gamma, double resolution) ;
        double ExpInt(double tlow, double thigh, double gamma, double resolution);
        double ExpCos(double t, double gamma, double deltaM, double resolution);
        double ExpCosInt(double tlow, double thigh, double gamma, double deltaM, double resolution );
        double ExpSin(double t, double gamma, double deltaM, double resolution);
        double ExpSinInt(double tlow, double thigh, double gamma, double deltaM, double resolution );
    };
    
    

}

#endif

