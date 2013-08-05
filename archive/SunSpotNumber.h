//
//  SunSpotNumber.h
//  SunSpotNumberViaMarkovNetwork
//
//  Created by Masoud Mirmomeni on 10/19/12.
//  Copyright (c) 2012 Masoud Mirmomeni. All rights reserved.
//

#ifndef SunSpotNumberViaMarkovNetwork_SunSpotNumber_h
#define SunSpotNumberViaMarkovNetwork_SunSpotNumber_h


#include <vector>

using namespace std;

class SunSpotNumber {
    
    
public:
    
    SunSpotNumber(int vec[] , int output = 0 , unsigned PredHor = 1 , unsigned InDim = 1)
    {
        for (int i = 0; i < InDim;i++)
            Input.push_back(vec[i]);
        Output = output;
        PredictionHorizon = PredHor;
        InputDimension = InDim;
    }
    SunSpotNumber(SunSpotNumber const& );
    SunSpotNumber& operator= (SunSpotNumber const&);
    
    vector<int> GetInput() const { return Input; }
    
    int GetOutput() const { return Output; }
    
    unsigned GetPredictionHorizon() const { return PredictionHorizon; }

    unsigned GetInputDimension() const { return InputDimension; }

    friend ostream& operator<<( ostream& Out, const SunSpotNumber& ssn );
    
private:
    vector<int> Input;
    int Output;
    unsigned PredictionHorizon;
    unsigned InputDimension;
};

int* IntToBinary(int );

#endif
