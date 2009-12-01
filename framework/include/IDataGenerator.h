/**
        @class IDataGenerator

        Interface for any class used for generating toy data from a PDF

 	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-10
*/


#ifndef I_DATA_GENERATOR_H
#define I_DATA_GENERATOR_H

#include "IDataSet.h"

class IDataGenerator
{
	public:
		virtual int GenerateData(int) = 0;
		virtual IDataSet * GetDataSet() = 0;
};

#endif