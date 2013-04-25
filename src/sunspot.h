/* sunspot.h
 *
 * This file is part of sunspot.
 *
 * Copyright 2012 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _SUNSPOT_H_
#define _SUNSPOT_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <ea/fitness_function.h>
#include <ea/meta_data.h>
#include <mkv/parse.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace boost::numeric::ublas;
using namespace ea;
using namespace std;
// Common meta-data needed for sunspot prediction.
LIBEA_MD_DECL(SUNSPOT_INPUT, "sunspot.input", std::string);

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;
    
	// create a working copy of the input
	matrix<T> A(input);
    
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
    
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;
    
	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));
    
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
    
	return true;
}

/*! Fitness function for sunspot number prediction.
 */
struct sunspot_fitness : fitness_function<unary_fitness<double>, constantS, absoluteS, stochasticS> {
    typedef boost::numeric::ublas::matrix<int> matrix_type; //!< Type for matrix that will store raw sunspot numbers.
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::vector<int> vector_type; //!< Type for a vector of sunspot numbers.

    mkv::markov_network::desc_type _desc; //!< Description of the markov network we'll be using.
    matrix_type _input; //!< All historical sunspot number data used during fitness evaluation (inputs to MKV network).
    matrix_type _observed; //!< Observed (real) historical sunspot numbers.
    vector_type _IntegerObserved;
    vector_type _IntegerInput;
    vector_type _IntegerObservedED;
    
    
    
    //! Estimating the embedding dimension.
	template <typename Embedding, typename Nonlinearity, typename EA>
	unsigned embedding_dimension(Embedding& d , Nonlinearity& n , EA& ea) {
        namespace bnu=boost::numeric::ublas;
        // input data can be found here (defined in config file or command line):
        std::string filename=get<SUNSPOT_INPUT>(ea);
        std::ifstream MyFile (filename.c_str());
        
        
        std::string Line;
        int MatrixSize=0;
        
        if (MyFile.is_open())
        {
            for (int i = 1; i <= 4; i++)
            {
                getline (MyFile,Line);
            }
            MyFile>>MatrixSize;
        }
        else
        {
            std::cerr<<"Sorry, this file cannot be read."<<std::endl;
            return 0;
        }
        
        _IntegerInput.resize(MatrixSize - 1);
        int TempSSN = 0;
        MyFile >>TempSSN;
        
        for (int i = 0; i < MatrixSize - 1; i++)
        {
            _IntegerInput(i) = TempSSN;
            MyFile >>TempSSN;
            _IntegerObservedED(i) = TempSSN;
            
        }
        
        const int MAX_ED = 7;
        const int MAX_NONLINEARITY = 4;
        matrix_type _Parameters;
        matrix_type _Training;
        matrix_type _TrainingEstimationMatrix;
        matrix_type _TrainVector;
        matrix_type _IntegerEstimatedED;
        vector_type _TrainError;
        
        
        _TrainError.resize(MAX_ED * MAX_NONLINEARITY);
        int ParameterOrder [MAX_NONLINEARITY][MAX_ED] = {{1,2,3,4,5,6,7},{8,10,13,17,22,28,35},{36,39,45,55,70,91,119},{120,124,134,154,189,245,329}};
        int NumParameters = 0;
        
        NumParameters = boost::math::factorial<int>(MAX_ED + MAX_NONLINEARITY) / (boost::math::factorial<int>(MAX_ED) * boost::math::factorial<int>(MAX_NONLINEARITY));
        _Training.resize(MatrixSize - MAX_ED - 1 , NumParameters);
        _TrainVector.resize(MatrixSize - MAX_ED - 1,1);
        _IntegerEstimatedED.resize(MatrixSize - MAX_ED - 1,1);;
        
        for (int i = 0; i < MatrixSize - MAX_ED - 1; i++)
        {
            _Training(i,0) = 1;
            
            for (int j = 1 ; j <= MAX_ED ; j++)
                _Training(i,j) = _IntegerObservedED(i + MAX_ED - j);
           
            _Training(i,8)   = _IntegerObservedED(i + MAX_ED - 1)^2;
            
            _Training(i,9)   = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,10)  = _IntegerObservedED(i + MAX_ED - 2)^2;
            
            _Training(i,11)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,12)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,13)  = _IntegerObservedED(i + MAX_ED - 3)^2;
            
            _Training(i,14)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,15)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,16)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,17)  = _IntegerObservedED(i + MAX_ED - 4)^2;
            
            _Training(i,18)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,19)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,20)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,21)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,22)  = _IntegerObservedED(i + MAX_ED - 5)^2;
            
            _Training(i,23)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,24)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,25)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,26)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,27)  = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,28)  = _IntegerObservedED(i + MAX_ED - 6)^2;
            
            _Training(i,29)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,30)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,31)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,32)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,33)  = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,34)  = _IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,35)  = _IntegerObservedED(i + MAX_ED - 7)^2;
            
            
            
            
            
            _Training(i,36)  = _IntegerObservedED(i + MAX_ED - 1)^3;
            
            _Training(i,37)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,38)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2;
            _Training(i,39)  = _IntegerObservedED(i + MAX_ED - 2)^3;
            
            _Training(i,40)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,41)  = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,42)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,43)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^2;
            _Training(i,44)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2;
            _Training(i,45)  = _IntegerObservedED(i + MAX_ED - 3)^3;
            
            _Training(i,46)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,47)  = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,48)  = _IntegerObservedED(i + MAX_ED - 3)^3*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,49)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,50)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,51)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,52)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,53)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,54)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,55)  = _IntegerObservedED(i + MAX_ED - 4)^3;
            
            _Training(i,56)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,57)  = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,58)  = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,59)  = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,60)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,61)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,62)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,63)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,64)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,65)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,66)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,67)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,68)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,69)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,70)  = _IntegerObservedED(i + MAX_ED - 5)^3;
            
            _Training(i,71)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,72)  = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,73)  = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,74)  = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,75)  = _IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,76)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,77)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,78)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,79)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,80)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,81)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,82)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,83)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,84)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,85)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,86)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,87)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,88)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,89)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,90)  = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,91)  = _IntegerObservedED(i + MAX_ED - 6)^3;
            
            _Training(i,92)  = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,93)  = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,94)  = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,95)  = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,96)  = _IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,97)  = _IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,98)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,99)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,100) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,101) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,102) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,103) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,104) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,105) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,106) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,107) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,108) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,109) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,110) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,111) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,112) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,113) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,114) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,115) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,116) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,117) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,118) = _IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,119) = _IntegerObservedED(i + MAX_ED - 7)^3;
            
            
            
            
            
            
            _Training(i,120) = _IntegerObservedED(i + MAX_ED - 1)^4;
            
            _Training(i,121) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,122) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)^2;
            _Training(i,123) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^3;
            _Training(i,124) = _IntegerObservedED(i + MAX_ED - 2)^4;
            
            _Training(i,125) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,126) = _IntegerObservedED(i + MAX_ED - 2)^3*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,127) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,128) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,129) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3)^2;
            _Training(i,130) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2;
            _Training(i,131) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3)^2;
            _Training(i,132) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^3;
            _Training(i,133) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^3;
            _Training(i,134) = _IntegerObservedED(i + MAX_ED - 3)^4;
            
            _Training(i,135) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,136) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,137) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,138) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,139) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,140) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,141) = _IntegerObservedED(i + MAX_ED - 2)^3*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,142) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,143) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,144) = _IntegerObservedED(i + MAX_ED - 3)^3*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,145) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,146) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,147) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,148) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,149) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,150) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4)^2;
            _Training(i,151) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)^3;
            _Training(i,152) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^3;
            _Training(i,153) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^3;
            _Training(i,154) = _IntegerObservedED(i + MAX_ED - 4)^4;
            
            _Training(i,155) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,156) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,157) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,158) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,159) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,160) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,161) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,162) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,163) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,164) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,165) = _IntegerObservedED(i + MAX_ED - 2)^3*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,166) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,167) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,168) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,169) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,170) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,171) = _IntegerObservedED(i + MAX_ED - 3)^3*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,172) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,173) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,174) = _IntegerObservedED(i + MAX_ED - 4)^3*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,175) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,176) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,177) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,178) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,179) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,180) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,181) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,182) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,183) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,184) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5)^2;
            _Training(i,185) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)^3;
            _Training(i,186) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)^3;
            _Training(i,187) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 3)^3;
            _Training(i,188) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^3;
            _Training(i,189) = _IntegerObservedED(i + MAX_ED - 5)^4;
            
            _Training(i,190) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,191) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,192) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,193) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,194) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,195) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,196) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,197) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,198) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,199) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,200) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,201) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,202) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,203) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,204) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,205) = _IntegerObservedED(i + MAX_ED - 2)^3*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,206) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,207) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,208) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,209) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,210) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,211) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,212) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,213) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,214) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,215) = _IntegerObservedED(i + MAX_ED - 3)^3*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,216) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,217) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,218) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,219) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,220) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,221) = _IntegerObservedED(i + MAX_ED - 4)^3*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,222) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,223) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,224) = _IntegerObservedED(i + MAX_ED - 5)^3*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,225) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,226) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,227) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,228) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,229) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,230) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,231) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,232) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,233) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,234) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,235) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,236) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,237) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,238) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,239) = _IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6)^2;
            _Training(i,240) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)^3;
            _Training(i,241) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)^3;
            _Training(i,242) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)^3;
            _Training(i,243) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^3;
            _Training(i,244) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^3;
            _Training(i,245) = _IntegerObservedED(i + MAX_ED - 6)^4;
            
            _Training(i,246) = _IntegerObservedED(i + MAX_ED - 1)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,247) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,248) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,249) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,250) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,251) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,252) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,253) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,254) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,255) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,256) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,257) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,258) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,259) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,260) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,261) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,262) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,263) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,264) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,265) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,266) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,267) = _IntegerObservedED(i + MAX_ED - 2)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,268) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,269) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,270) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,271) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,272) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,273) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,274) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,275) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,276) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,277) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,278) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,279) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,280) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,281) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,282) = _IntegerObservedED(i + MAX_ED - 3)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,283) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,284) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,285) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,286) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,287) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,288) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,289) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,290) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,291) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,292) = _IntegerObservedED(i + MAX_ED - 4)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,293) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,294) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,295) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,296) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,297) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,298) = _IntegerObservedED(i + MAX_ED - 5)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,299) = _IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,300) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,301) = _IntegerObservedED(i + MAX_ED - 6)^3*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,302) = _IntegerObservedED(i + MAX_ED - 1)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,303) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,304) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,305) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,306) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,307) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,308) = _IntegerObservedED(i + MAX_ED - 2)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,309) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,310) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,311) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,312) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,313) = _IntegerObservedED(i + MAX_ED - 3)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,314) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,315) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,316) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,317) = _IntegerObservedED(i + MAX_ED - 4)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,318) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,319) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,320) = _IntegerObservedED(i + MAX_ED - 5)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,321) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,322) = _IntegerObservedED(i + MAX_ED - 6)^2*_IntegerObservedED(i + MAX_ED - 7)^2;
            _Training(i,323) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,324) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,325) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,326) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,327) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,328) = _IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7)^3;
            _Training(i,329) = _IntegerObservedED(i + MAX_ED - 7)^4;
            _TrainVector(i,1)  = _IntegerObservedED(i + MAX_ED);
        }
        
        matrix_type TempMatrix(MatrixSize - MAX_ED - 1 , 1 , 1);
        int ColumnCounter;
    
        for (int j = 1; j <= MAX_NONLINEARITY; j++)
        {
            for (int i = 1; i <= MAX_ED; i++)
            {
                NumParameters = boost::math::factorial<int>(i + j) / (boost::math::factorial<int>(i) * boost::math::factorial<int>(j));
                _Parameters.resize(NumParameters,1);
                _TrainingEstimationMatrix.resize(MatrixSize - MAX_ED - 1 , NumParameters);
                ColumnCounter = 0;
                column(_TrainingEstimationMatrix, ColumnCounter) = column(TempMatrix, 1);
                
                for (int p = 0; p < j ; p++)
                {
                    for (int k = ParameterOrder[p][0];k <= ParameterOrder[p][i-1];k++)
                    {
                        ColumnCounter++;
                        column(_TrainingEstimationMatrix, ColumnCounter) = column(_Training, k);
                    }
                }
                                
                matrix_type _TrainingTranspose = boost::numeric::ublas::trans(_TrainingEstimationMatrix);
                matrix_type _TrainingSquare    = boost::numeric::ublas::prod(_TrainingTranspose, _TrainingEstimationMatrix);
                matrix_type _TrainingInverse;
                _TrainingInverse.resize(NumParameters , NumParameters);
                
                InvertMatrix(_TrainingSquare, _TrainingInverse);
                matrix_type _LeftMatrix = boost::numeric::ublas::prod(_TrainingTranspose, _TrainingEstimationMatrix);
                _Parameters = boost::numeric::ublas::prod(_LeftMatrix, _TrainVector);
                
                
                /*****************
                 Error Calculation
                 *****************/
                
                _IntegerEstimatedED = boost::numeric::ublas::prod (_TrainingEstimationMatrix , _Parameters);
                bnu::vector<double> err = column(_IntegerEstimatedED,1) - column(_TrainVector , 1);
                _TrainError((j - 1) * MAX_ED + (i - 1)) = sqrt(1/static_cast<double>(err.size()) * bnu::inner_prod(err,err));
            }
        }
    
    unsigned f = 0;
    double MinError = 100000000000000000000000000;
    
    for (unsigned i = 0 ; i < _TrainError.size() ; i++)
    {
        if (_TrainError(i) < MinError)
        {
            MinError = _TrainError(i);
            f = i % MAX_ED + 1;
            d = i % MAX_ED + 1;
            n = i / MAX_ED + 1;

        }
    }
        
    return f;
    }
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    //! Initialize this fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        // this simply parses the geometry of the MKV network from the configuration
        // file:
        mkv::parse_desc(get<MKV_DESC>(ea), _desc);
        
        // input data can be found here (defined in config file or command line):
        std::string filename=get<SUNSPOT_INPUT>(ea);
        std::ifstream MyFile (filename.c_str());
        
        
        std::string Line;
        int MatrixSize=0;
        
        if (MyFile.is_open())
        {
            for (int i = 1; i <= 4; i++)
            {
                getline (MyFile,Line);
            }
            MyFile>>MatrixSize;
        }
        else
        {
            std::cerr<<"Sorry, this file cannot be read."<<std::endl;
            return;
        }
        
        const int NumberOfDigits = 9;
        
        // read in the historical data for sunspot numbers, and split it into:

        // _input: a matrix where each row vector i is assumed to be the complete
        // **binary input vector** to the MKV network at time i.
        _input.resize(MatrixSize,NumberOfDigits); // dummy initializer; replace with real size and data
        
        // _observed: a vector where element i corresponds to the observed (real) sunspot
        // number corresponding to row i in _inputs.
        _observed.resize(MatrixSize,NumberOfDigits); // dummy initializer; replace with real size and data
        _IntegerObserved.resize(MatrixSize);
        int TempSSN = 0;
        
        for (int i = 0; i < MatrixSize; i++)
        {
            MyFile >>TempSSN;
            _IntegerObserved(i) = TempSSN;
            std::bitset<NumberOfDigits> TempBinarySSN = ~std::bitset<NumberOfDigits>(_IntegerObserved(i));
            
            for (int j = 0; j < NumberOfDigits; j++)
            {
                _input(i,j)=TempBinarySSN[NumberOfDigits - j - 1];
            }
            
            MyFile >>TempSSN;
            _IntegerObserved(i) = TempSSN;
            TempBinarySSN = ~std::bitset<NumberOfDigits>(_IntegerObserved(i));
            
            for (int j = 0; j < NumberOfDigits; j++)
            {
                _observed(i,j)=TempBinarySSN[NumberOfDigits - j - 1];
            }
        }
        
        
    }
    
    //! Calculate fitness of ind.
	template <typename Individual, typename RNG, typename EA>
	double operator()(Individual& ind, RNG& rng, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        mkv::markov_network net = mkv::make_markov_network(_desc, ind.repr().begin(), ind.repr().end(), rng, ea);

        // vector of outputs from the MKV network, initialized to the same size
        // as R:
        vector_type output(_observed.size1());
        
        // run each row of _inputs through the MKV network for a single update,
        // place the result in the output vector:
        for(std::size_t i=0; i<_input.size1(); ++i) {
            
            row_type r(_input,i);
            
            mkv::update(net, 1, r.begin());
            // convert the binary output from the MKV network to an integer;
            // this uses the +/- encoding that has been shown to be effective with
            // MKV networks:
            output(i) = static_cast<double>(algorithm::range_pair2int(net.begin_output(), net.end_output()));
        }
    
        // fitness is 1.0/(1.0+sqrt((observed-output)^2)) -- RMSE:
        bnu::vector<int> err = _IntegerObserved - output;
        double f = 100.0/(1.0+sqrt(1.0/static_cast<double>(err.size()) * bnu::inner_prod(err,err)));
        return f;
    }
};
        
        
#endif
