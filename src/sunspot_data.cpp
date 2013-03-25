//
//  main.cpp
//  SunSpotNumberViaMarkovNetwork
//
//  Created by Masoud Mirmomeni on 10/19/12.
//  Copyright (c) 2012 Masoud Mirmomeni. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <vector>
#include "SunSpotNumber.h"
#include <string>
#include <cstring>
#include <sstream>
#include <array>

using namespace std;

int main(int argc, const char * argv[])
{
    
    cout.setf(ios::boolalpha);
    // insert code here...
    vector<SunSpotNumber> SSN;
    int Number , TotalNumber = 0 , Counter = 0;
    double * SunSpotTotal;
    
    vector <string> VectorLine;
    string Line;
    int Year , InitialYear , FinalYear;
    cout << "Enter the years you want to read (Beginning Ending): ";
    cin>> InitialYear >> FinalYear;
    
    for (Year = InitialYear ; Year <=FinalYear ; Year++)
    {
        if (Year % 4 == 0)
            TotalNumber += 366;
        else
            TotalNumber += 365;
    }
    
    SunSpotTotal = new (nothrow) double [TotalNumber];
    
    for (Year = InitialYear ; Year <=FinalYear ; Year++)
    {
        bool LeapYear;
        
        
        if (Year % 4 == 0)
        {
            LeapYear = true;
            Number = 366;
        }
        
        else
        {
            LeapYear = false;
            Number = 365;
        }
        
        string FileAddress = "/Users/mirmomeny/Desktop/SunSpotNumberViaMarkovNetwork/SunSpotNumber";
        
        string YearString;       // string which will contain the year
        
        ostringstream Convert;   // stream used for the conversion
        
        Convert << Year;
        
        YearString = Convert.str();
        
        FileAddress = FileAddress + YearString;
        
        ifstream MyFile (FileAddress);
        string SubString;
        double * SunSpot;
        
        if (MyFile.is_open())
        {
            for (int i = 1; i <= 4; i++)
            {
                getline (MyFile,Line);
            }
            
            if (LeapYear)
            {
                
                int MonthsDay[12] = {0,31,60,91,121,152,182,213,244,274,305,335};
                SunSpot = new (nothrow) double [366];
                
                for (int i = 0 ; i < 29 ; i++)
                {
                    getline (MyFile,Line);
                    istringstream LineStream(Line);
                    
                    for (int j = 0 ; j < 12 ; j++)
                    {
                        LineStream >> SubString;
                        stringstream StringSubStream (SubString);
                        int Temp;
                        StringSubStream >> Temp;
                        SunSpot[i + MonthsDay[j]] = Temp;
                        
                    }
                }
                
                int MonthsDayTemp1[11] = {29,89,120,150,181,211,242,273,303,334,364};
                
                getline (MyFile,Line);
                istringstream LineStream(Line);
                
                for (int j = 0 ; j < 11 ; j++)
                {
                    LineStream >> SubString;
                    stringstream StringSubStream (SubString);
                    int Temp;
                    StringSubStream >> Temp;
                    SunSpot[MonthsDayTemp1[j]] = Temp;
                    
                }
                
                int MonthsDayTemp2[7] = {30,90,151,212,243,304,365};
                
                getline (MyFile,Line);
                istringstream LineStream2(Line);
                
                for (int j = 0 ; j < 7 ; j++)
                {
                    LineStream2 >> SubString;
                    stringstream StringSubStream (SubString);
                    int Temp;
                    StringSubStream >> Temp;
                    SunSpot[MonthsDayTemp2[j]] = Temp;
                    
                }
                
            }
            else
            {
                int MonthsDay[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
                
                SunSpot = new (nothrow) double [365];
                
                for (int i = 0 ; i < 28 ; i++)
                {
                    getline (MyFile,Line);
                    istringstream LineStream(Line);
                    
                    for (int j = 0 ; j < 12 ; j++)
                    {
                        LineStream >> SubString;
                        stringstream StringSubStream (SubString);
                        int Temp;
                        StringSubStream >> Temp;
                        SunSpot[i + MonthsDay[j]] = Temp;
                        
                    }
                }
                
                int MonthsDayTemp1[11] = {28,87,118,148,179,209,240,271,301,331,362};
                
                for (int i = 0 ; i < 2 ; i++)
                {
                    getline (MyFile,Line);
                    istringstream LineStream(Line);
                    
                    for (int j = 0 ; j < 11 ; j++)
                    {
                        LineStream >> SubString;
                        stringstream StringSubStream (SubString);
                        int Temp;
                        StringSubStream >> Temp;
                        SunSpot[i + MonthsDayTemp1[j]] = Temp;
                        
                    }
                }
                
                int MonthsDayTemp2[7] = {30,89,150,211,242,303,364};
                
                getline (MyFile,Line);
                istringstream LineStream2(Line);
                
                for (int j = 0 ; j < 7 ; j++)
                {
                    LineStream2 >> SubString;
                    stringstream StringSubStream (SubString);
                    int Temp;
                    StringSubStream >> Temp;
                    SunSpot[MonthsDayTemp2[j]] = Temp;
                    
                }
            }
            
            for (int i = 0; i < Number ; i++)
                SunSpotTotal[Counter + i] = SunSpot[i];
            
            Counter += Number;
            
            MyFile.close();
            delete [] SunSpot;
            
        }
        else
        {
            cout << "Unable to open file \n";
            break;
        }
    }
    
    cout << "What is the dimension of the input vector (e.g. 5)? ";
    int InputDimension;
    cin >> InputDimension;
    
    cout << "What is the prediction horizon (e.g. 1)? ";
    int PredictionHorizon;
    cin >> PredictionHorizon;
    
    
    int TempArray[InputDimension];
    
    for (int i = 0; i < TotalNumber - InputDimension - PredictionHorizon + 1; i++)
    {
        
        for (int j = 0 ; j < InputDimension ; j++)
            TempArray[j] = SunSpotTotal[i + j];
        int TempOutput = SunSpotTotal[i + InputDimension + PredictionHorizon - 1];
        SSN.push_back(SunSpotNumber(TempArray , TempOutput , PredictionHorizon , InputDimension));
    }
    
    delete [] SunSpotTotal;
    
    ofstream MyFile;
    MyFile.open ("MarkovNetworkBinarykData.txt");
    MyFile <<"This file contains the integer representation of sunspot numbers."<<endl;
    MyFile<<"The first (n-1)th integers of each row represent the input SSNs for Markov network, and the last integer represents the corresponding output to those inputs (prediction target)."<<endl;
    MyFile<<"This file contains sunspot number from "<<InitialYear<<" to " <<FinalYear<<"."<<endl;
    
    
    MyFile <<"Number of data in this file is:"<<endl;
    MyFile<<SSN.size()<<endl;
    
    
    SunSpotNumber TempSSN(TempArray);
    vector<int> TempInput;
    int TempOutput;
    
    for (int i = 0 ; i <SSN.size() ; i++)
    {
        TempSSN=SSN[i];
        TempInput = TempSSN.GetInput();
        TempOutput = TempSSN.GetOutput();
    
        for (int j = 0 ; j < TempSSN.GetInputDimension() ; j++)
        {
            
            /* int* TempBinary;
            TempBinary = IntToBinary(TempInput[j]);
            
            for (int k = 0; k<9; k++) {
                MyFile <<TempBinary[k]<<" ";
            } */
            
            MyFile << TempInput[j] << " ";
            
            
        }
        
        /* int* TempBinary = IntToBinary(TempOutput);
                
        for (int k = 0; k<8; k++) {
            MyFile <<TempBinary[k]<<" ";
        }*/
        
        MyFile <<TempOutput<<endl;
    }
    
    MyFile.close();
    cout<<"Preprocessing is done :-)!" << endl;
    return 0;
}

