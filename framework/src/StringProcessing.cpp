/**
  @class StringProcessing

  Collection of static functions for string processing

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#include "StringProcessing.h"
#include <iostream>

//Split a string every time you find a given character
vector<string> StringProcessing::SplitString( string Input, char SplitCharacter )
{
	vector<string> splitParts;

	while(true)
	{
		//Search for the character
		int position = CharacterPosition( Input, SplitCharacter );

		if ( position == -1 )
		{
			//Ignore empty strings
			if ( Input != "" )
			{
				splitParts.push_back( Input );
			}

			//If it's not found, you've reached the end
			break;
		}
		else
		{
			//Split the string at the character postion
			string tempString( Input, 0, position );
			string newInput( Input, position + 1 );

			//Ignore empty strings
			if ( tempString != "" )
			{
				splitParts.push_back( tempString );
			}
			Input = newInput;
		}
	}

	return splitParts;
}

//Return the position of the first instance of a character in a string
int StringProcessing::CharacterPosition( string Input, char SearchCharacter )
{
	for (int characterIndex = 0; characterIndex < Input.size(); characterIndex++)
	{
		//Look for the character
		if ( Input[characterIndex] == SearchCharacter )
		{
			return characterIndex;
		}

		//If you get to the end, character not found
		if ( characterIndex == Input.size() - 1 )
		{
			return -1;
		}
	}
}

//Return the position of all instances of a string in another string
vector<int> StringProcessing::StringPositions( string Input, string SearchString )
{
	vector<int> positions;
	int numberDiscarded = 0;

	while (true)
	{
		int firstCharacterPosition = CharacterPosition( Input, SearchString[0] );
		if ( firstCharacterPosition == -1 )
		{
			//If you can't find the first character of the string, you won't find the string!
			return positions;
		}
		else
		{
			bool found = true;
			for ( int characterIndex = 0; characterIndex < SearchString.size(); characterIndex++ )
			{
				//Compare each subsequent character
				if ( Input[ characterIndex + firstCharacterPosition ] != SearchString[characterIndex] )
				{
					found = false;
					break;
				}
			}

			if (found)
			{
				positions.push_back( firstCharacterPosition + numberDiscarded );
			}

			//Search the rest of the string
			string tempString( Input, firstCharacterPosition + 1 );
			numberDiscarded += firstCharacterPosition + 1;
			Input = tempString;
		}
	}

	return positions;
}

//Remove any instances of a particular character in a string
void StringProcessing::RemoveCharacter( string & Input, char SearchCharacter )
{
	char passedCharacters[ Input.size() + 1 ];
	int addedCharacters = 0;

	for ( int characterIndex = 0; characterIndex < Input.size(); characterIndex++ )
	{
		if ( Input[characterIndex] != SearchCharacter )
		{
			passedCharacters[addedCharacters] = Input[characterIndex];
			addedCharacters++;
		}
	}

	passedCharacters[addedCharacters] = '\0';
	Input = passedCharacters;
}

//Remove white space from passed lines
void StringProcessing::RemoveWhiteSpace( vector<string> & newContent )
{
	vector<string> output;

	for ( int lineIndex = 0; lineIndex < newContent.size(); lineIndex++ )
	{
		//Remove tabs
		RemoveCharacter( newContent[lineIndex], '\t' );

		//Remove empty lines
		if ( newContent[lineIndex] != "" )
		{
			output.push_back( newContent[lineIndex] );
		}
	}

	newContent = output;
}

//Return a vector containing all the unique strings from the two input vectors
vector<string> StringProcessing::CombineUniques( vector<string> VectorOne, vector<string> VectorTwo )
{
	vector<string> result;
	vector<string>::iterator stringIterator;

	//Don't assume VectorOne is unique
	for ( stringIterator = VectorOne.begin(); stringIterator != VectorOne.end(); stringIterator++ )
	{
		if ( VectorContains( &result, &(*stringIterator) ) == -1 )
		{
			result.push_back( *stringIterator );
		}
	}

	//Now add in VectorTwo
	for ( stringIterator = VectorTwo.begin(); stringIterator != VectorTwo.end(); stringIterator++ )
	{
		if ( VectorContains( &result, &(*stringIterator) ) == -1 )
		{
			result.push_back( *stringIterator );
		}
	}

	return result;
}

//Return the position of a search string within a vector of strings, or -1 if not found
int StringProcessing::VectorContains( vector<string> * InputVector, string * SearchString )
{
	for ( int searchIndex = 0; searchIndex < InputVector->size(); searchIndex++ )
	{
		if ( (*InputVector)[searchIndex] == (*SearchString) )
		{
			//Found the search string
			return searchIndex;
		}
	}

	//If you've got this far, it wasn't found
	return -1;
}