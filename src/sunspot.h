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
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <cmath>
#include <ea/fitness_function.h>
#include <ea/meta_data.h>
#include <mkv/parse.h>
#include <fstream>
#include <string>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>

/*#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
*/

using namespace boost::numeric::ublas;
using namespace ea;
using namespace std;
// using namespace Eigen;
// Common meta-data needed for sunspot prediction.
LIBEA_MD_DECL(SUNSPOT_INPUT, "sunspot.input", std::string);

int MatrixSize;
const int MAX_ED = 7;
const int MAX_NONLINEARITY = 4;
int ParameterOrder [MAX_NONLINEARITY][MAX_ED] = {{1,2,3,4,5,6,7},{8,10,13,17,22,28,35},{36,39,45,55,70,91,119},{120,124,134,154,189,245,329}};


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
    typedef boost::numeric::ublas::matrix<double> matrix_type_estimated;
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::vector<int> vector_type; //!< Type for a vector of sunspot numbers.
    typedef boost::numeric::ublas::vector<double> vector_type_distance; //!< Type for a vector of sunspot numbers.
    
    mkv::markov_network::desc_type _desc; //!< Description of the markov network we'll be using.
    matrix_type _input; //!< All historical sunspot number data used during fitness evaluation (inputs to MKV network).
    matrix_type           _observed; //!< Observed (real) historical sunspot numbers.
    vector_type           _IntegerObserved;
    vector_type           _IntegerInput;
    vector_type_distance  _IntegerObservedED;
    matrix_type_estimated _Training;
    
    
    
    // Estimating the embedding dimension.
	template <typename Embedding, typename Nonlinearity, typename EA>
	unsigned embedding_dimension(Embedding& d , Nonlinearity& n , EA& ea) {
        namespace bnu=boost::numeric::ublas;
        // input data can be found here (defined in config file or command line):
        std::string filename=get<SUNSPOT_INPUT>(ea);
        std::ifstream MyFile (filename.c_str());
        
        
        std::string Line;
        MatrixSize=0;
        
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
        
        matrix_type_estimated _Parameters;
        matrix_type           _TrainingEstimationMatrix;
        matrix_type_estimated _IntegerEstimatedED;
        vector_type_distance  _TrainError;
        matrix_type           _TrainVector;
        
        
        _TrainError.resize(MAX_ED * MAX_NONLINEARITY);
        int NumParameters = 0;
        
        NumParameters = boost::math::factorial<int>(MAX_ED + MAX_NONLINEARITY) / (boost::math::factorial<int>(MAX_ED) * boost::math::factorial<int>(MAX_NONLINEARITY));
        _Training.resize(MatrixSize - MAX_ED - 1 , NumParameters);
        _TrainVector.resize(MatrixSize - MAX_ED - 1,1);
        _IntegerEstimatedED.resize(MatrixSize - MAX_ED - 1,1);;
        
        for (int i = 0; i < MatrixSize - MAX_ED - 1; i++)
        {
            _Training(i,0) = 1.0;
            
            for (int j = 1 ; j <= MAX_ED ; j++)
                _Training(i,j) = _IntegerObservedED(i + MAX_ED - j);
           
            _Training(i,8)   = pow(_IntegerObservedED(i + MAX_ED - 1),2);
            
            _Training(i,9)   = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,10)  = pow(_IntegerObservedED(i + MAX_ED - 2),2);
            
            _Training(i,11)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,12)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,13)  = pow(_IntegerObservedED(i + MAX_ED - 3),2);
            
            _Training(i,14)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,15)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,16)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,17)  = pow(_IntegerObservedED(i + MAX_ED - 4),2);
            
            _Training(i,18)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,19)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,20)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,21)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,22)  = pow(_IntegerObservedED(i + MAX_ED - 5),2);
            
            _Training(i,23)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,24)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,25)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,26)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,27)  = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,28)  = pow(_IntegerObservedED(i + MAX_ED - 6),2);
            
            _Training(i,29)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,30)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,31)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,32)  = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,33)  = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,34)  = _IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,35)  = pow(_IntegerObservedED(i + MAX_ED - 7),2);
            
            
            
            
            
            _Training(i,36)  = pow(_IntegerObservedED(i + MAX_ED - 1),3);
            
            _Training(i,37)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,38)  = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2);
            _Training(i,39)  = pow(_IntegerObservedED(i + MAX_ED - 2),3);
            
            _Training(i,40)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,41)  = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,42)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,43)  = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),2);
            _Training(i,44)  = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2);
            _Training(i,45)  = pow(_IntegerObservedED(i + MAX_ED - 3),3);
            
            _Training(i,46)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,47)  = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,48)  = pow(_IntegerObservedED(i + MAX_ED - 3),3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,49)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,50)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,51)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,52)  = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,53)  = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,54)  = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,55)  = pow(_IntegerObservedED(i + MAX_ED - 4),3);
            
            _Training(i,56)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,57)  = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,58)  = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,59)  = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,60)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,61)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,62)  = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,63)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,64)  = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,65)  = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,66)  = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,67)  = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,68)  = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,69)  = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,70)  = pow(_IntegerObservedED(i + MAX_ED - 5),3);
            
            _Training(i,71)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,72)  = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,73)  = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,74)  = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,75)  = pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6);
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
            _Training(i,86)  = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,87)  = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,88)  = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,89)  = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,90)  = _IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,91)  = pow(_IntegerObservedED(i + MAX_ED - 6),3);
            
            _Training(i,92)  = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,93)  = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,94)  = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,95)  = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,96)  = pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,97)  = pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
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
            _Training(i,113) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,114) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,115) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,116) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,117) = _IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,118) = _IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,119) = pow(_IntegerObservedED(i + MAX_ED - 7),3);
            
            
            
            
            
            
            _Training(i,120) = pow(_IntegerObservedED(i + MAX_ED - 1),4);
            
            _Training(i,121) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 2);
            _Training(i,122) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 2),2);
            _Training(i,123) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),3);
            _Training(i,124) = pow(_IntegerObservedED(i + MAX_ED - 2),4);
            
            _Training(i,125) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,126) = pow(_IntegerObservedED(i + MAX_ED - 2),3)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,127) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,128) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3);
            _Training(i,129) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 3),2);
            _Training(i,130) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2);
            _Training(i,131) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*pow(_IntegerObservedED(i + MAX_ED - 3),2);
            _Training(i,132) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),3);
            _Training(i,133) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),3);
            _Training(i,134) = pow(_IntegerObservedED(i + MAX_ED - 3),4);
            
            _Training(i,135) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,136) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,137) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,138) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,139) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,140) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,141) = pow(_IntegerObservedED(i + MAX_ED - 2),3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,142) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,143) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,144) = pow(_IntegerObservedED(i + MAX_ED - 3),3)*_IntegerObservedED(i + MAX_ED - 4);
            _Training(i,145) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,146) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,147) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,148) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,149) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,150) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*pow(_IntegerObservedED(i + MAX_ED - 4),2);
            _Training(i,151) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 4),3);
            _Training(i,152) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),3);
            _Training(i,153) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),3);
            _Training(i,154) = pow(_IntegerObservedED(i + MAX_ED - 4),4);
            
            _Training(i,155) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,156) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,157) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,158) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,159) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,160) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,161) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,162) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,163) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,164) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,165) = pow(_IntegerObservedED(i + MAX_ED - 2),3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,166) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,167) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,168) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,169) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,170) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,171) = pow(_IntegerObservedED(i + MAX_ED - 3),3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,172) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,173) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,174) = pow(_IntegerObservedED(i + MAX_ED - 4),3)*_IntegerObservedED(i + MAX_ED - 5);
            _Training(i,175) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,176) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,177) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,178) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,179) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,180) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,181) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,182) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,183) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,184) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*pow(_IntegerObservedED(i + MAX_ED - 5),2);
            _Training(i,185) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 5),3);
            _Training(i,186) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 5),3);
            _Training(i,187) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 3),3);
            _Training(i,188) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),3);
            _Training(i,189) = pow(_IntegerObservedED(i + MAX_ED - 5),4);
            
            _Training(i,190) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,191) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,192) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,193) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,194) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,195) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,196) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,197) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,198) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,199) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,200) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,201) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,202) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,203) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,204) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,205) = pow(_IntegerObservedED(i + MAX_ED - 2),3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,206) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,207) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,208) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,209) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,210) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,211) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,212) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,213) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,214) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,215) = pow(_IntegerObservedED(i + MAX_ED - 3),3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,216) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,217) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,218) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,219) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,220) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,221) = pow(_IntegerObservedED(i + MAX_ED - 4),3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,222) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,223) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,224) = pow(_IntegerObservedED(i + MAX_ED - 5),3)*_IntegerObservedED(i + MAX_ED - 6);
            _Training(i,225) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,226) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,227) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,228) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,229) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,230) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,231) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,232) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,233) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,234) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,235) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,236) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,237) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,238) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,239) = pow(_IntegerObservedED(i + MAX_ED - 5),2)*pow(_IntegerObservedED(i + MAX_ED - 6),2);
            _Training(i,240) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 6),3);
            _Training(i,241) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 6),3);
            _Training(i,242) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 6),3);
            _Training(i,243) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),3);
            _Training(i,244) = _IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),3);
            _Training(i,245) = pow(_IntegerObservedED(i + MAX_ED - 6),4);
            
            _Training(i,246) = pow(_IntegerObservedED(i + MAX_ED - 1),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,247) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,248) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,249) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,250) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,251) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
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
            _Training(i,262) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,263) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,264) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,265) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,266) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,267) = pow(_IntegerObservedED(i + MAX_ED - 2),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,268) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,269) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,270) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,271) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,272) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,273) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,274) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,275) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,276) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,277) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,278) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,279) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,280) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,281) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,282) = pow(_IntegerObservedED(i + MAX_ED - 3),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,283) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,284) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,285) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,286) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,287) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,288) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,289) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,290) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,291) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,292) = pow(_IntegerObservedED(i + MAX_ED - 4),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,293) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,294) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,295) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,296) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,297) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,298) = pow(_IntegerObservedED(i + MAX_ED - 5),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,299) = pow(_IntegerObservedED(i + MAX_ED - 5),2)*_IntegerObservedED(i + MAX_ED - 6)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,300) = _IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 6),2)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,301) = pow(_IntegerObservedED(i + MAX_ED - 6),3)*_IntegerObservedED(i + MAX_ED - 7);
            _Training(i,302) = pow(_IntegerObservedED(i + MAX_ED - 1),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,303) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,304) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,305) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,306) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,307) = _IntegerObservedED(i + MAX_ED - 1)*_IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,308) = pow(_IntegerObservedED(i + MAX_ED - 2),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,309) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,310) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,311) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,312) = _IntegerObservedED(i + MAX_ED - 2)*_IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,313) = pow(_IntegerObservedED(i + MAX_ED - 3),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,314) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,315) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,316) = _IntegerObservedED(i + MAX_ED - 3)*_IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,317) = pow(_IntegerObservedED(i + MAX_ED - 4),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,318) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,319) = _IntegerObservedED(i + MAX_ED - 4)*_IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,320) = pow(_IntegerObservedED(i + MAX_ED - 5),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,321) = _IntegerObservedED(i + MAX_ED - 5)*_IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,322) = pow(_IntegerObservedED(i + MAX_ED - 6),2)*pow(_IntegerObservedED(i + MAX_ED - 7),2);
            _Training(i,323) = _IntegerObservedED(i + MAX_ED - 1)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,324) = _IntegerObservedED(i + MAX_ED - 2)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,325) = _IntegerObservedED(i + MAX_ED - 3)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,326) = _IntegerObservedED(i + MAX_ED - 4)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,327) = _IntegerObservedED(i + MAX_ED - 5)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,328) = _IntegerObservedED(i + MAX_ED - 6)*pow(_IntegerObservedED(i + MAX_ED - 7),3);
            _Training(i,329) = pow(_IntegerObservedED(i + MAX_ED - 7),4);
            
            
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
                matrix_type_estimated _TrainingSquare    = boost::numeric::ublas::prod(_TrainingTranspose, _TrainingEstimationMatrix);
                matrix_type_estimated _TrainingInverse;
                _TrainingInverse.resize(NumParameters , NumParameters);
                
                InvertMatrix(_TrainingSquare, _TrainingInverse);
                matrix_type_estimated _LeftMatrix = boost::numeric::ublas::prod(_TrainingTranspose, _TrainingEstimationMatrix);
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
    double MinError = std::numeric_limits<double>::max();
    
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
    
    
    
    // Estimating the Lyapunov exponent first approach.
	template <typename Embedding, typename Nonlinearity, typename Lyapunov, typename EA>
	double lyapunov_estimation(Embedding& d , Nonlinearity& n , Lyapunov& l, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        // input data can be found here (defined in config file or command line):
        
        /**************
         Initialization
         **************/
        
        int NumParameters = boost::math::factorial<int>(d + n) / (boost::math::factorial<int>(d) * boost::math::factorial<int>(n));
        
        matrix_type_estimated _Parameters;
        matrix_type_estimated _Regressor;
        matrix_type_estimated _RegressorTranspose;
        matrix_type_estimated _Variance;
        matrix_type_estimated _TempOne;
        matrix_type_estimated _TempTwo;
        matrix_type_estimated _TempThree;
        matrix_type_estimated _TempFour;
        matrix_type_estimated _Denominator;
        matrix_type_estimated _Identity;
        matrix_type_estimated _InverseDenominator;
        matrix_type_estimated _EvolutionTerm;
        matrix_type_estimated _Jacobian;
        
        _Jacobian.resize(d,d);
        
        
        for (int i = 0; i < d - 1; i++)
        {
            for (int j = 0; j < d; j++)
            {
                _Jacobian(i,j)=0;
            }
            _Jacobian(i,i+1)=1;
        }
        
        _Parameters.resize(NumParameters,1);
        _Identity.resize(1,1);
        _Identity(0,0)=1;
        
        matrix_type_estimated _ActualInput;
        matrix_type_estimated _EstimatedInput;
        matrix_type_estimated _Error;
        
        _ActualInput.resize(1,1);
        _EstimatedInput.resize(1,1);
        _Error.resize(1,1);
        
        
        _Regressor.resize(NumParameters,1);
        _Variance.resize(NumParameters,NumParameters);
        
        srand((unsigned)time(0));
        
        for(int i=0; i<NumParameters; i++)
            _Parameters(i,0) = (rand()%10)+1;
            
        for(int i=0; i<NumParameters; i++){
            
            for(int j=0; j<NumParameters; j++)
                _Variance(i,j)=0;
            _Variance(i,i)=1000000;
        }
        
        /**********
         Estimation
         **********/
        
        for (int i = 0; i < MatrixSize - MAX_ED - 1; i++)
        {
            
            /********************
             Parameter Estimation
             ********************/
            
            int ColumnCounter = 0;
            _Regressor(0,0)   = 1;
            
            for (int p = 0; p < n ; p++)
            {
                for (int k = ParameterOrder[p][0];k <= ParameterOrder[p][d-1];k++)
                {
                    ColumnCounter++;
                    _Regressor(ColumnCounter , 0) = _Training(i , k);
                }
            }
            
            
            
            _ActualInput(0,0) = _IntegerObservedED(i + MAX_ED);
            
            _RegressorTranspose = boost::numeric::ublas::trans(_Regressor);
            _TempOne            = boost::numeric::ublas::prod(_RegressorTranspose, _Variance);
            _TempTwo            = boost::numeric::ublas::prod(_Variance , _Regressor);
            _TempThree          = boost::numeric::ublas::prod(_TempTwo, _TempOne);
            _TempFour           = boost::numeric::ublas::prod(_RegressorTranspose, _TempTwo);
            _Denominator        = _Identity + _TempFour;
            InvertMatrix(_Denominator, _InverseDenominator);
            _EvolutionTerm = boost::numeric::ublas::prod(_InverseDenominator, _TempThree);
            _Variance -= _EvolutionTerm;
            matrix_type_estimated _TransParameters = boost::numeric::ublas::trans(_Parameters);
            _EstimatedInput = boost::numeric::ublas::prod(_TransParameters , _Regressor);
            
            _Error = _ActualInput - _EstimatedInput;
            
            matrix_type_estimated _EvolutionParameterOne = boost::numeric::ublas::prod(_Variance , _Regressor);
            matrix_type_estimated _EvolutionParameter    = boost::numeric::ublas::prod(_Error , _EvolutionParameterOne);
            _Parameters += _EvolutionParameter;
            
            /*******************
             Jacobian Estimation
             *******************/
            
            switch (n)
            {
                case 1:
                    
                    switch (d)
                    {
                        case 1:
                            _Jacobian(0,0) = _Parameters(1,0);
                            break;
                        case 2:
                            _Jacobian(1,0) = _Parameters(1,0);
                            _Jacobian(1,1) = _Parameters(2,0);
                            break;
                        case 3:
                            _Jacobian(2,0) = _Parameters(1,0);
                            _Jacobian(2,1) = _Parameters(2,0);
                            _Jacobian(2,2) = _Parameters(3,0);
                            break;
                        case 4:
                            _Jacobian(3,0) = _Parameters(1,0);
                            _Jacobian(3,1) = _Parameters(2,0);
                            _Jacobian(3,2) = _Parameters(3,0);
                            _Jacobian(3,3) = _Parameters(4,0);
                            break;
                        case 5:
                            _Jacobian(4,0) = _Parameters(1,0);
                            _Jacobian(4,1) = _Parameters(2,0);
                            _Jacobian(4,2) = _Parameters(3,0);
                            _Jacobian(4,3) = _Parameters(4,0);
                            _Jacobian(4,4) = _Parameters(5,0);
                            break;
                        case 6:
                            _Jacobian(5,0) = _Parameters(1,0);
                            _Jacobian(5,1) = _Parameters(2,0);
                            _Jacobian(5,2) = _Parameters(3,0);
                            _Jacobian(5,3) = _Parameters(4,0);
                            _Jacobian(5,4) = _Parameters(5,0);
                            _Jacobian(5,5) = _Parameters(6,0);
                            break;
                        case 7:
                            _Jacobian(6,0) = _Parameters(1,0);
                            _Jacobian(6,1) = _Parameters(2,0);
                            _Jacobian(6,2) = _Parameters(3,0);
                            _Jacobian(6,3) = _Parameters(4,0);
                            _Jacobian(6,4) = _Parameters(5,0);
                            _Jacobian(6,5) = _Parameters(6,0);
                            _Jacobian(6,6) = _Parameters(7,0);
                            break;
                    }
                    break;
                    
                case 2:
                    switch (d)
                    {
                    case 1:
                        _Jacobian(0,0) = _Parameters(1,0) + 2 * _Parameters(2,0) * _Training(i,1);
                        break;
                    case 2:
                        _Jacobian(1,0) = _Parameters(1,0) + 2 * _Parameters(3,0) * _Training(i,1) + _Parameters(4,0) * _Training(i,2);
                        _Jacobian(1,1) = _Parameters(2,0) + _Parameters(4,0) * _Training(i,1) + 2 * _Parameters(5,0) * _Training(i,2);
                        break;
                    case 3:
                        _Jacobian(2,0) = _Parameters(1,0)+2*_Parameters(4,0)*_Training(i,1)+_Parameters(5,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,3);
                        _Jacobian(2,1) = _Parameters(2,0)+2*_Parameters(6,0)*_Training(i,2)+_Parameters(5,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,3);
                        _Jacobian(2,2) = _Parameters(3,0)+2*_Parameters(9,0)*_Training(i,3)+_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2);
                        break;
                    case 4:
                        _Jacobian(3,0) = _Parameters(1,0)+2*_Parameters(5,0)*_Training(i,1)+_Parameters(6,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,3)+_Parameters(11,0)*_Training(i,4);
                        _Jacobian(3,1) = _Parameters(2,0)+2*_Parameters(7,0)*_Training(i,2)+_Parameters(6,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4);
                        _Jacobian(3,2) = _Parameters(3,0)+2*_Parameters(10,0)*_Training(i,3)+_Parameters(8,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,2)+_Parameters(13,0)*_Training(i,4);
                        _Jacobian(3,3) = _Parameters(4,0)+2*_Parameters(14,0)*_Training(i,4)+_Parameters(11,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,2)+_Parameters(13,0)*_Training(i,4);
                        break;
                    case 5:
                        _Jacobian(4,0) = _Parameters(1,0)+2*_Parameters(6,0)*_Training(i,1)+_Parameters(7,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4)+_Parameters(16,0)*_Training(i,5);
                        _Jacobian(4,1) = _Parameters(2,0)+2*_Parameters(8,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5);
                        _Jacobian(4,2) = _Parameters(3,0)+2*_Parameters(11,0)*_Training(i,3)+_Parameters(9,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5);
                        _Jacobian(4,3) = _Parameters(4,0)+2*_Parameters(15,0)*_Training(i,4)+_Parameters(12,0)*_Training(i,1)+_Parameters(13,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5);
                        _Jacobian(4,4) = _Parameters(5,0)+2*_Parameters(20,0)*_Training(i,5)+_Parameters(16,0)*_Training(i,1)+_Parameters(17,0)*_Training(i,2)+_Parameters(18,0)*_Training(i,3)+_Parameters(19,0)*_Training(i,4);
                        break;
                    case 6:
                        _Jacobian(5,0) = _Parameters(1,0)+2*_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5)+_Parameters(22,0)*_Training(i,6);
                        _Jacobian(5,1) = _Parameters(2,0)+2*_Parameters(9,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6);
                        _Jacobian(5,2) = _Parameters(3,0)+2*_Parameters(12,0)*_Training(i,3)+_Parameters(10,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,2)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6);
                        _Jacobian(5,3) = _Parameters(4,0)+2*_Parameters(16,0)*_Training(i,4)+_Parameters(13,0)*_Training(i,1)+_Parameters(14,0)*_Training(i,2)+_Parameters(15,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6);
                        _Jacobian(5,4) = _Parameters(5,0)+2*_Parameters(21,0)*_Training(i,5)+_Parameters(17,0)*_Training(i,1)+_Parameters(18,0)*_Training(i,2)+_Parameters(19,0)*_Training(i,3)+_Parameters(20,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,6);
                        _Jacobian(5,5) = _Parameters(6,0)+2*_Parameters(27,0)*_Training(i,6)+_Parameters(22,0)*_Training(i,1)+_Parameters(23,0)*_Training(i,2)+_Parameters(24,0)*_Training(i,3)+_Parameters(25,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,5);
                        break;
                    case 7:
                        _Jacobian(5,0) = _Parameters(1,0)+2*_Parameters(8,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,2)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6)+_Parameters(29,0)*_Training(i,7);
                        _Jacobian(5,1) = _Parameters(2,0)+2*_Parameters(10,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,3)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6)+_Parameters(30,0)*_Training(i,7);
                        _Jacobian(5,2) = _Parameters(3,0)+2*_Parameters(13,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,2)+_Parameters(16,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6)+_Parameters(31,0)*_Training(i,7);
                        _Jacobian(5,3) = _Parameters(4,0)+2*_Parameters(17,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,1)+_Parameters(15,0)*_Training(i,2)+_Parameters(16,0)*_Training(i,4)+_Parameters(21,0)*_Training(i,5)+_Parameters(26,0)*_Training(i,6)+_Parameters(32,0)*_Training(i,7);
                        _Jacobian(5,4) = _Parameters(5,0)+2*_Parameters(22,0)*_Training(i,5)+_Parameters(18,0)*_Training(i,1)+_Parameters(19,0)*_Training(i,2)+_Parameters(20,0)*_Training(i,3)+_Parameters(21,0)*_Training(i,4)+_Parameters(27,0)*_Training(i,6)+_Parameters(33,0)*_Training(i,7);
                        _Jacobian(5,5) = _Parameters(6,0)+2*_Parameters(28,0)*_Training(i,6)+_Parameters(23,0)*_Training(i,1)+_Parameters(24,0)*_Training(i,2)+_Parameters(25,0)*_Training(i,3)+_Parameters(26,0)*_Training(i,4)+_Parameters(27,0)*_Training(i,5)+_Parameters(34,0)*_Training(i,7);
                        _Jacobian(6,6) = _Parameters(7,0)+2*_Parameters(35,0)*_Training(i,7)+_Parameters(29,0)*_Training(i,1)+_Parameters(30,0)*_Training(i,2)+_Parameters(31,0)*_Training(i,3)+_Parameters(32,0)*_Training(i,4)+_Parameters(33,0)*_Training(i,5)+_Parameters(34,0)*_Training(i,6);
                        break;
                    }
                    break;
                    
                case 3:
                    switch (d)
                    {
                    case 1:
                            _Jacobian(0,0) = _Parameters(1,0)+2*_Parameters(2,0)*_Training(i,1)+3*_Parameters(3,0)*pow(_Training(i,1),2);
                            break;
                    case 2:
                            _Jacobian(1,0) = _Parameters(1,0)+2*_Parameters(3,0)*_Training(i,1)+_Parameters(4,0)*_Training(i,2)+3*_Parameters(6,0)*pow(_Training(i,1),2)+2*_Parameters(7,0)*_Training(i,1)*_Training(i,2)+_Parameters(8,0)*pow(_Training(i,2),2);
                            
                            _Jacobian(1,1) = _Parameters(2,0)+2*_Parameters(5,0)*_Training(i,5)+_Parameters(4,0)*_Training(i,1)+3*_Parameters(9,0)*pow(_Training(i,2),2)+2*_Parameters(8,0)*_Training(i,1)*_Training(i,2)+_Parameters(7,0)*pow(_Training(i,1),2);
                            break;
                    case 3:
                        
                        _Jacobian(2,0) = _Parameters(1,0)+2*_Parameters(4,0)*_Training(i,1)+_Parameters(5,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,3)+3*_Parameters(10,0)*pow(_Training(i,1),2)+2*_Parameters(11,0)*_Training(i,1)*_Training(i,2)+_Parameters(12,0)*pow(_Training(i,2),2)+2*_Parameters(14,0)*_Training(i,1)*_Training(i,3)+_Parameters(16,0)*_Training(i,2)*_Training(i,3)+_Parameters(17,0)*pow(_Training(i,3),2);
                        
                        _Jacobian(2,1) = _Parameters(2,0)+2*_Parameters(6,0)*_Training(i,2)+_Parameters(5,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,3)+3*_Parameters(13,0)*pow(_Training(i,2),2)+2*_Parameters(12,0)*_Training(i,1)*_Training(i,2)+_Parameters(11,0)*pow(_Training(i,1),2)+2*_Parameters(15,0)*_Training(i,2)*_Training(i,3)+_Parameters(16,0)*_Training(i,1)*_Training(i,3)+_Parameters(18,0)*pow(_Training(i,3),2);
                        
                        _Jacobian(2,2) = _Parameters(3,0)+2*_Parameters(9,0)*_Training(i,3)+_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2)+3*_Parameters(19,0)*pow(_Training(i,3),2)+2*_Parameters(17,0)*_Training(i,1)*_Training(i,3)+_Parameters(14,0)*pow(_Training(i,1),2)+2*_Parameters(18,0)*_Training(i,2)*_Training(i,3)+_Parameters(16,0)*_Training(i,1)*_Training(i,2)+_Parameters(15,0)*pow(_Training(i,2),2);
                        break;
                    case 4:
                        _Jacobian(3,0) = _Parameters(1,0)+2*_Parameters(5,0)*_Training(i,1)+_Parameters(6,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,3)+_Parameters(11,0)*_Training(i,4)+3*_Parameters(15,0)*pow(_Training(i,1),2)+2*_Parameters(16,0)*_Training(i,1)*_Training(i,2)+_Parameters(17,0)*pow(_Training(i,2),2)+2*_Parameters(19,0)*_Training(i,1)*_Training(i,3)+_Parameters(21,0)*_Training(i,2)*_Training(i,3)+_Parameters(22,0)*pow(_Training(i,3),2)+2*_Parameters(25,0)*_Training(i,1)*_Training(i,4)+_Parameters(28,0)*_Training(i,2)*_Training(i,4)+_Parameters(29,0)*_Training(i,3)*_Training(i,4)+_Parameters(31,0)*pow(_Training(i,4),2);
                        
                        _Jacobian(3,1) = _Parameters(2,0)+2*_Parameters(7,0)*_Training(i,2)+_Parameters(6,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4)+3*_Parameters(18,0)*pow(_Training(i,2),2)+2*_Parameters(17,0)*_Training(i,1)*_Training(i,2)+_Parameters(16,0)*pow(_Training(i,1),2)+2*_Parameters(20,0)*_Training(i,2)*_Training(i,3)+_Parameters(21,0)*_Training(i,1)*_Training(i,3)+_Parameters(23,0)*pow(_Training(i,3),2)+2*_Parameters(26,0)*_Training(i,2)*_Training(i,4)+_Parameters(28,0)*_Training(i,1)*_Training(i,4)+_Parameters(30,0)*_Training(i,3)*_Training(i,4)+_Parameters(32,0)*pow(_Training(i,4),2);
                        
                        _Jacobian(3,2) = _Parameters(3,0)+2*_Parameters(10,0)*_Training(i,3)+_Parameters(9,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,1)+_Parameters(13,0)*_Training(i,4)+3*_Parameters(24,0)*pow(_Training(i,3),2)+2*_Parameters(22,0)*_Training(i,1)*_Training(i,3)+_Parameters(19,0)*pow(_Training(i,1),2)+2*_Parameters(23,0)*_Training(i,2)*_Training(i,3)+_Parameters(21,0)*_Training(i,1)*_Training(i,2)+_Parameters(20,0)*pow(_Training(i,2),2)+2*_Parameters(27,0)*_Training(i,3)*_Training(i,4)+_Parameters(29,0)*_Training(i,1)*_Training(i,4)+_Parameters(30,0)*_Training(i,2)*_Training(i,4)+_Parameters(33,0)*pow(_Training(i,4),2);
                        
                        _Jacobian(3,3) = _Parameters(4,0)+_Parameters(11,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,2)+_Parameters(13,0)*_Training(i,3)+2*_Parameters(14,0)*_Training(i,4)+_Parameters(25,0)*pow(_Training(i,1),2)+_Parameters(26,0)*pow(_Training(i,2),2)+_Parameters(27,0)*pow(_Training(i,3),2)+_Parameters(28,0)*_Training(i,1)*_Training(i,2)+_Parameters(29,0)*_Training(i,1)*_Training(i,3)+_Parameters(30,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(31,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(32,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(33,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(34,0)*pow(_Training(i,4),2);
                        break;
                    case 5:
                        _Jacobian(4,0) = _Parameters(1,0)+2*_Parameters(6,0)*_Training(i,1)+_Parameters(7,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4)+_Parameters(16,0)*_Training(i,5)+3*_Parameters(21,0)*pow(_Training(i,1),2)+2*_Parameters(22,0)*_Training(i,1)*_Training(i,2)+_Parameters(23,0)*pow(_Training(i,2),2)+2*_Parameters(25,0)*_Training(i,1)*_Training(i,3)+_Parameters(27,0)*_Training(i,2)*_Training(i,3)+_Parameters(28,0)*pow(_Training(i,3),2)+2*_Parameters(31,0)*_Training(i,1)*_Training(i,4)+_Parameters(34,0)*_Training(i,2)*_Training(i,4)+_Parameters(35,0)*_Training(i,3)*_Training(i,4)+_Parameters(37,0)*pow(_Training(i,4),2)+2*_Parameters(41,0)*_Training(i,1)*_Training(i,5)+_Parameters(45,0)*_Training(i,2)*_Training(i,5)+_Parameters(46,0)*_Training(i,3)*_Training(i,5)+_Parameters(47,0)*_Training(i,4)*_Training(i,5)+_Parameters(51,0)*pow(_Training(i,5),2);
                        
                        _Jacobian(4,1) = _Parameters(2,0)+2*_Parameters(8,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5)+3*_Parameters(24,0)*pow(_Training(i,2),2)+2*_Parameters(23,0)*_Training(i,1)*_Training(i,2)+_Parameters(22,0)*pow(_Training(i,1),2)+2*_Parameters(26,0)*_Training(i,2)*_Training(i,3)+_Parameters(27,0)*_Training(i,1)*_Training(i,3)+_Parameters(29,0)*pow(_Training(i,3),2)+2*_Parameters(32,0)*_Training(i,2)*_Training(i,4)+_Parameters(34,0)*_Training(i,1)*_Training(i,4)+_Parameters(36,0)*_Training(i,3)*_Training(i,4)+_Parameters(38,0)*pow(_Training(i,4),2)+2*_Parameters(42,0)*_Training(i,2)*_Training(i,5)+_Parameters(45,0)*_Training(i,1)*_Training(i,5)+_Parameters(48,0)*_Training(i,3)*_Training(i,5)+_Parameters(49,0)*_Training(i,4)*_Training(i,5)+_Parameters(52,0)*pow(_Training(i,5),2);
                        
                        _Jacobian(4,2) = _Parameters(3,0)+2*_Parameters(11,0)*_Training(i,3)+_Parameters(9,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+3*_Parameters(30,0)*pow(_Training(i,3),2)+2*_Parameters(28,0)*_Training(i,1)*_Training(i,3)+_Parameters(25,0)*pow(_Training(i,1),2)+2*_Parameters(29,0)*_Training(i,2)*_Training(i,3)+_Parameters(27,0)*_Training(i,1)*_Training(i,2)+_Parameters(26,0)*pow(_Training(i,2),2)+2*_Parameters(33,0)*_Training(i,3)*_Training(i,4)+_Parameters(35,0)*_Training(i,1)*_Training(i,4)+_Parameters(36,0)*_Training(i,2)*_Training(i,4)+_Parameters(39,0)*pow(_Training(i,4),2)+2*_Parameters(43,0)*_Training(i,3)*_Training(i,5)+_Parameters(46,0)*_Training(i,1)*_Training(i,5)+_Parameters(48,0)*_Training(i,2)*_Training(i,5)+_Parameters(50,0)*_Training(i,4)*_Training(i,5)+_Parameters(53,0)*pow(_Training(i,5),2);
                        
                        _Jacobian(4,3) = _Parameters(4,0)+2*_Parameters(15,0)*_Training(i,4)+_Parameters(12,0)*_Training(i,1)+_Parameters(13,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,3)+_Parameters(19,0)*_Training(i,5)+3*_Parameters(40,0)*pow(_Training(i,4),2)+2*_Parameters(37,0)*_Training(i,1)*_Training(i,4)+_Parameters(31,0)*pow(_Training(i,1),2)+2*_Parameters(38,0)*_Training(i,2)*_Training(i,4)+_Parameters(34,0)*_Training(i,1)*_Training(i,2)+_Parameters(32,0)*pow(_Training(i,2),2)+2*_Parameters(39,0)*_Training(i,3)*_Training(i,4)+_Parameters(35,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*_Training(i,2)*_Training(i,3)+_Parameters(33,0)*pow(_Training(i,3),2)+2*_Parameters(44,0)*_Training(i,4)*_Training(i,5)+_Parameters(47,0)*_Training(i,1)*_Training(i,5)+_Parameters(49,0)*_Training(i,2)*_Training(i,5)+_Parameters(50,0)*_Training(i,3)*_Training(i,5)+_Parameters(54,0)*pow(_Training(i,5),2);
                            
                        _Jacobian(4,4) = _Parameters(5,0)+2*_Parameters(20,0)*_Training(i,5)+_Parameters(16,0)*_Training(i,1)+_Parameters(17,0)*_Training(i,2)+_Parameters(18,0)*_Training(i,3)+_Parameters(19,0)*_Training(i,4)+3*_Parameters(55,0)*pow(_Training(i,5),2)+2*_Parameters(51,0)*_Training(i,1)*_Training(i,5)+_Parameters(41,0)*pow(_Training(i,1),2)+2*_Parameters(52,0)*_Training(i,2)*_Training(i,5)+_Parameters(45,0)*_Training(i,1)*_Training(i,2)+_Parameters(42,0)*pow(_Training(i,2),2)+2*_Parameters(53,0)*_Training(i,3)*_Training(i,5)+_Parameters(46,0)*_Training(i,1)*_Training(i,3)+_Parameters(48,0)*_Training(i,2)*_Training(i,3)+_Parameters(43,0)*pow(_Training(i,3),2)+2*_Parameters(54,0)*_Training(i,4)*_Training(i,5)+_Parameters(47,0)*_Training(i,1)*_Training(i,4)+_Parameters(49,0)*_Training(i,2)*_Training(i,4)+_Parameters(50,0)*_Training(i,3)*_Training(i,4)+_Parameters(44,0)*pow(_Training(i,4),2);
                            
                        break;
                    case 6:
                        _Jacobian(5,0) = _Parameters(1,0)+2*_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5)+_Parameters(22,0)*_Training(i,6)+3*_Parameters(28,0)*pow(_Training(i,1),2)+2*_Parameters(29,0)*_Training(i,1)*_Training(i,2)+_Parameters(30,0)*pow(_Training(i,2),2)+2*_Parameters(32,0)*_Training(i,1)*_Training(i,2)+_Parameters(34,0)*_Training(i,2)*_Training(i,3)+_Parameters(35,0)*pow(_Training(i,3),2)+2*_Parameters(38,0)*_Training(i,1)*_Training(i,4)+_Parameters(41,0)*_Training(i,2)*_Training(i,4)+_Parameters(42,0)*_Training(i,3)*_Training(i,4)+_Parameters(44,0)*pow(_Training(i,4),2)+2*_Parameters(48,0)*_Training(i,1)*_Training(i,5)+_Parameters(52,0)*_Training(i,2)*_Training(i,5)+_Parameters(53,0)*_Training(i,3)*_Training(i,5)+_Parameters(54,0)*_Training(i,4)*_Training(i,5)+_Parameters(58,0)*pow(_Training(i,5),2)+2*_Parameters(63,0)*_Training(i,1)*_Training(i,6)+_Parameters(68,0)*_Training(i,2)*_Training(i,6)+_Parameters(69,0)*_Training(i,3)*_Training(i,6)+_Parameters(70,0)*_Training(i,4)*_Training(i,6)+_Parameters(71,0)*_Training(i,5)*_Training(i,6)+_Parameters(78,0)*pow(_Training(i,6),2);
                            
                        _Jacobian(5,1) = _Parameters(2,0)+2*_Parameters(9,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6)+3*_Parameters(31,0)*pow(_Training(i,2),2)+2*_Parameters(30,0)*_Training(i,1)*_Training(i,2)+_Parameters(29,0)*pow(_Training(i,1),2)+2*_Parameters(33,0)*_Training(i,2)*_Training(i,3)+_Parameters(34,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*pow(_Training(i,3),2)+2*_Parameters(39,0)*_Training(i,2)*_Training(i,4)+_Parameters(41,0)*_Training(i,1)*_Training(i,4)+_Parameters(43,0)*_Training(i,3)*_Training(i,4)+_Parameters(45,0)*pow(_Training(i,4),2)+2*_Parameters(49,0)*_Training(i,2)*_Training(i,5)+_Parameters(52,0)*_Training(i,1)*_Training(i,5)+_Parameters(55,0)*_Training(i,3)*_Training(i,5)+_Parameters(56,0)*_Training(i,4)*_Training(i,5)+_Parameters(59,0)*pow(_Training(i,5),2)+2*_Parameters(64,0)*_Training(i,2)*_Training(i,6)+_Parameters(68,0)*_Training(i,1)*_Training(i,6)+_Parameters(72,0)*_Training(i,3)*_Training(i,6)+_Parameters(73,0)*_Training(i,4)*_Training(i,6)+_Parameters(74,0)*_Training(i,5)*_Training(i,6)+_Parameters(79,0)*pow(_Training(i,6),2);
                            
                        _Jacobian(5,2) = _Parameters(3,0)+_Parameters(10,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,2)+2*_Parameters(12,0)*_Training(i,3)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6)+_Parameters(32,0)*pow(_Training(i,1),2)+_Parameters(33,0)*pow(_Training(i,2),2)+_Parameters(34,0)*_Training(i,1)*_Training(i,2)+2*_Parameters(35,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*_Training(i,2)*_Training(i,3)+3*_Parameters(37,0)*pow(_Training(i,3),2)+2*_Parameters(40,0)*_Training(i,3)*_Training(i,4)+_Parameters(42,0)*_Training(i,1)*_Training(i,4)+_Parameters(43,0)*_Training(i,2)*_Training(i,4)+_Parameters(46,0)*pow(_Training(i,4),2)+2*_Parameters(50,0)*_Training(i,3)*_Training(i,5)+_Parameters(53,0)*_Training(i,1)*_Training(i,5)+_Parameters(55,0)*_Training(i,2)*_Training(i,5)+_Parameters(57,0)*_Training(i,4)*_Training(i,5)+_Parameters(60,0)*pow(_Training(i,5),2)+2*_Parameters(65,0)*_Training(i,3)*_Training(i,6)+_Parameters(69,0)*_Training(i,1)*_Training(i,6)+_Parameters(72,0)*_Training(i,2)*_Training(i,6)+_Parameters(75,0)*_Training(i,4)*_Training(i,6)+_Parameters(76,0)*_Training(i,5)*_Training(i,6)+_Parameters(80,0)*pow(_Training(i,6),2);
                            
                        _Jacobian(5,3) = _Parameters(4,0)+_Parameters(13,0)*_Training(i,1)+_Parameters(14,0)*_Training(i,2)+_Parameters(15,0)*_Training(i,3)+2*_Parameters(16,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6)+_Parameters(38,0)*pow(_Training(i,1),2)+_Parameters(39,0)*pow(_Training(i,2),2)+_Parameters(40,0)*pow(_Training(i,3),2)+_Parameters(41,0)*_Training(i,1)*_Training(i,2)+_Parameters(42,0)*_Training(i,1)*_Training(i,3)+_Parameters(43,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(44,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(45,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(46,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(47,0)*pow(_Training(i,4),2)+2*_Parameters(51,0)*_Training(i,4)*_Training(i,5)+_Parameters(54,0)*_Training(i,1)*_Training(i,5)+_Parameters(56,0)*_Training(i,2)*_Training(i,5)+_Parameters(57,0)*_Training(i,3)*_Training(i,5)+_Parameters(61,0)*pow(_Training(i,5),2)+2*_Parameters(66,0)*_Training(i,4)*_Training(i,6)+_Parameters(70,0)*_Training(i,1)*_Training(i,6)+_Parameters(73,0)*_Training(i,2)*_Training(i,6)+_Parameters(75,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,5)*_Training(i,6)+_Parameters(81,0)*pow(_Training(i,6),2);
                            
                        _Jacobian(5,4) = _Parameters(5,0)+_Parameters(17,0)*_Training(i,1)+_Parameters(18,0)*_Training(i,2)+_Parameters(19,0)*_Training(i,3)+2*_Parameters(21,0)*_Training(i,5)+_Parameters(20,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,6)+_Parameters(48,0)*pow(_Training(i,1),2)+_Parameters(49,0)*pow(_Training(i,2),2)+_Parameters(50,0)*pow(_Training(i,3),2)+_Parameters(51,0)*pow(_Training(i,4),2)+_Parameters(52,0)*_Training(i,1)*_Training(i,2)+_Parameters(53,0)*_Training(i,1)*_Training(i,3)+_Parameters(54,0)*_Training(i,1)*_Training(i,4)+_Parameters(55,0)*_Training(i,2)*_Training(i,3)+_Parameters(56,0)*_Training(i,2)*_Training(i,4)+_Parameters(57,0)*_Training(i,3)*_Training(i,4)+2*_Parameters(58,0)*_Training(i,1)*_Training(i,5)+2*_Parameters(59,0)*_Training(i,2)*_Training(i,5)+2*_Parameters(60,0)*_Training(i,3)*_Training(i,5)+2*_Parameters(61,0)*_Training(i,4)*_Training(i,5)+3*_Parameters(62,0)*pow(_Training(i,5),2)+2*_Parameters(67,0)*_Training(i,5)*_Training(i,6)+_Parameters(71,0)*_Training(i,1)*_Training(i,6)+_Parameters(74,0)*_Training(i,2)*_Training(i,6)+_Parameters(76,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,4)*_Training(i,6)+_Parameters(82,0)*pow(_Training(i,6),2);
                            
                        _Jacobian(5,5) = _Parameters(6,0)+_Parameters(22,0)*_Training(i,1)+_Parameters(23,0)*_Training(i,2)+_Parameters(24,0)*_Training(i,3)+2*_Parameters(27,0)*_Training(i,6)+_Parameters(25,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,5)+_Parameters(63,0)*pow(_Training(i,1),2)+_Parameters(64,0)*pow(_Training(i,2),2)+_Parameters(65,0)*pow(_Training(i,3),2)+_Parameters(66,0)*pow(_Training(i,4),2)+_Parameters(68,0)*_Training(i,1)*_Training(i,2)+_Parameters(69,0)*_Training(i,1)*_Training(i,3)+_Parameters(70,0)*_Training(i,1)*_Training(i,4)+_Parameters(71,0)*_Training(i,1)*_Training(i,5)+_Parameters(72,0)*_Training(i,2)*_Training(i,3)+_Parameters(73,0)*_Training(i,2)*_Training(i,4)+_Parameters(74,0)*_Training(i,2)*_Training(i,5)+_Parameters(75,0)*_Training(i,3)*_Training(i,4)+_Parameters(76,0)*_Training(i,3)*_Training(i,5)+_Parameters(77,0)*_Training(i,4)*_Training(i,5)+2*_Parameters(78,0)*_Training(i,1)*_Training(i,6)+2*_Parameters(79,0)*_Training(i,2)*_Training(i,6)+2*_Parameters(80,0)*_Training(i,3)*_Training(i,6)+2*_Parameters(81,0)*_Training(i,4)*_Training(i,6)+2*_Parameters(85,0)*_Training(i,5)*_Training(i,6)+2*_Parameters(67,0)*pow(_Training(i,5),2)+3*_Parameters(83,0)*pow(_Training(i,6),2);
                            
                        break;
                    case 7:
                        _Jacobian(6,0) = _Parameters(1,0)+2*_Parameters(8,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,2)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6)+_Parameters(29,0)*_Training(i,7)+3*_Parameters(36,0)*pow(_Training(i,1),2)+2*_Parameters(37,0)*_Training(i,1)*_Training(i,2)+_Parameters(38,0)*pow(_Training(i,2),2)+2*_Parameters(40,0)*_Training(i,1)*_Training(i,3)+_Parameters(42,0)*_Training(i,2)*_Training(i,3)+_Parameters(43,0)*pow(_Training(i,3),2)+2*_Parameters(46,0)*_Training(i,1)*_Training(i,4)+_Parameters(49,0)*_Training(i,2)*_Training(i,4)+_Parameters(50,0)*_Training(i,3)*_Training(i,4)+_Parameters(52,0)*pow(_Training(i,4),2)+2*_Parameters(56,0)*_Training(i,1)*_Training(i,5)+_Parameters(60,0)*_Training(i,2)*_Training(i,5)+_Parameters(61,0)*_Training(i,3)*_Training(i,5)+_Parameters(62,0)*_Training(i,4)*_Training(i,5)+_Parameters(66,0)*pow(_Training(i,5),2)+2*_Parameters(71,0)*_Training(i,1)*_Training(i,6)+_Parameters(76,0)*_Training(i,2)*_Training(i,6)+_Parameters(77,0)*_Training(i,3)*_Training(i,6)+_Parameters(78,0)*_Training(i,4)*_Training(i,6)+_Parameters(79,0)*_Training(i,5)*_Training(i,6)+_Parameters(86,0)*pow(_Training(i,6),2)+2*_Parameters(92,0)*_Training(i,1)*_Training(i,7)+_Parameters(98,0)*_Training(i,2)*_Training(i,7)+_Parameters(99,0)*_Training(i,3)*_Training(i,7)+_Parameters(100,0)*_Training(i,4)*_Training(i,7)+_Parameters(101,0)*_Training(i,5)*_Training(i,7)+_Parameters(102,0)*_Training(i,6)*_Training(i,7)+_Parameters(113,0)*pow(_Training(i,7),2);
                            
                        _Jacobian(6,1) = _Parameters(2,0)+2*_Parameters(10,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,3)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6)+_Parameters(30,0)*_Training(i,7)+3*_Parameters(39,0)*pow(_Training(i,2),2)+2*_Parameters(38,0)*_Training(i,1)*_Training(i,2)+_Parameters(37,0)*pow(_Training(i,1),2)+2*_Parameters(41,0)*_Training(i,2)*_Training(i,3)+_Parameters(42,0)*_Training(i,1)*_Training(i,3)+_Parameters(44,0)*pow(_Training(i,3),2)+2*_Parameters(47,0)*_Training(i,2)*_Training(i,4)+_Parameters(49,0)*_Training(i,1)*_Training(i,4)+_Parameters(51,0)*_Training(i,3)*_Training(i,4)+_Parameters(53,0)*pow(_Training(i,4),2)+2*_Parameters(57,0)*_Training(i,2)*_Training(i,5)+_Parameters(60,0)*_Training(i,1)*_Training(i,5)+_Parameters(63,0)*_Training(i,3)*_Training(i,5)+_Parameters(64,0)*_Training(i,4)*_Training(i,5)+_Parameters(67,0)*pow(_Training(i,5),2)+2*_Parameters(72,0)*_Training(i,2)*_Training(i,6)+_Parameters(76,0)*_Training(i,1)*_Training(i,6)+_Parameters(80,0)*_Training(i,3)*_Training(i,6)+_Parameters(81,0)*_Training(i,4)*_Training(i,6)+_Parameters(82,0)*_Training(i,5)*_Training(i,6)+_Parameters(87,0)*pow(_Training(i,6),2)+2*_Parameters(93,0)*_Training(i,2)*_Training(i,7)+_Parameters(98,0)*_Training(i,1)*_Training(i,7)+_Parameters(103,0)*_Training(i,3)*_Training(i,7)+_Parameters(104,0)*_Training(i,4)*_Training(i,7)+_Parameters(105,0)*_Training(i,5)*_Training(i,7)+_Parameters(106,0)*_Training(i,6)*_Training(i,7)+_Parameters(114,0)*pow(_Training(i,7),2);
                            
                        _Jacobian(6,2) = _Parameters(3,0)+2*_Parameters(13,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,2)+_Parameters(11,0)*_Training(i,1)+_Parameters(16,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6)+_Parameters(31,0)*_Training(i,7)+3*_Parameters(45,0)*pow(_Training(i,3),2)+2*_Parameters(44,0)*_Training(i,3)*_Training(i,2)+_Parameters(41,0)*pow(_Training(i,2),2)+2*_Parameters(43,0)*_Training(i,1)*_Training(i,3)+_Parameters(42,0)*_Training(i,1)*_Training(i,2)+_Parameters(40,0)*pow(_Training(i,1),2)+2*_Parameters(48,0)*_Training(i,3)*_Training(i,4)+_Parameters(50,0)*_Training(i,1)*_Training(i,4)+_Parameters(51,0)*_Training(i,2)*_Training(i,4)+_Parameters(54,0)*pow(_Training(i,4),2)+2*_Parameters(58,0)*_Training(i,3)*_Training(i,5)+_Parameters(61,0)*_Training(i,1)*_Training(i,5)+_Parameters(63,0)*_Training(i,2)*_Training(i,5)+_Parameters(65,0)*_Training(i,4)*_Training(i,5)+_Parameters(68,0)*pow(_Training(i,5),2)+2*_Parameters(73,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,1)*_Training(i,6)+_Parameters(80,0)*_Training(i,2)*_Training(i,6)+_Parameters(83,0)*_Training(i,4)*_Training(i,6)+_Parameters(84,0)*_Training(i,5)*_Training(i,6)+_Parameters(88,0)*pow(_Training(i,6),2)+2*_Parameters(94,0)*_Training(i,3)*_Training(i,7)+_Parameters(99,0)*_Training(i,1)*_Training(i,7)+_Parameters(103,0)*_Training(i,2)*_Training(i,7)+_Parameters(107,0)*_Training(i,4)*_Training(i,7)+_Parameters(108,0)*_Training(i,5)*_Training(i,7)+_Parameters(109,0)*_Training(i,6)*_Training(i,7)+_Parameters(115,0)*pow(_Training(i,7),2);
                            

                        _Jacobian(6,3) = _Parameters(4,0)+_Parameters(14,0)*_Training(i,1)+_Parameters(15,0)*_Training(i,2)+_Parameters(16,0)*_Training(i,3)+2*_Parameters(17,0)*_Training(i,4)+_Parameters(21,0)*_Training(i,5)+_Parameters(26,0)*_Training(i,6)+_Parameters(32,0)*_Training(i,7)+_Parameters(46,0)*pow(_Training(i,1),2)+_Parameters(47,0)*pow(_Training(i,2),2)+_Parameters(48,0)*pow(_Training(i,3),2)+_Parameters(49,0)*_Training(i,1)*_Training(i,2)+_Parameters(50,0)*_Training(i,1)*_Training(i,3)+_Parameters(51,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(52,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(53,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(53,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(55,0)*pow(_Training(i,4),2)+2*_Parameters(59,0)*_Training(i,4)*_Training(i,5)+_Parameters(62,0)*_Training(i,1)*_Training(i,5)+_Parameters(64,0)*_Training(i,2)*_Training(i,5)+_Parameters(65,0)*_Training(i,3)*_Training(i,5)+_Parameters(69,0)*pow(_Training(i,5),2)+2*_Parameters(74,0)*_Training(i,4)*_Training(i,6)+_Parameters(78,0)*_Training(i,1)*_Training(i,6)+_Parameters(81,0)*_Training(i,2)*_Training(i,6)+_Parameters(83,0)*_Training(i,3)*_Training(i,6)+_Parameters(85,0)*_Training(i,5)*_Training(i,6)+_Parameters(89,0)*pow(_Training(i,6),2)+2*_Parameters(95,0)*_Training(i,1)*_Training(i,7)+_Parameters(100,0)*_Training(i,1)*_Training(i,7)+_Parameters(104,0)*_Training(i,2)*_Training(i,7)+_Parameters(107,0)*_Training(i,3)*_Training(i,7)+_Parameters(110,0)*_Training(i,5)*_Training(i,7)+_Parameters(111,0)*_Training(i,6)*_Training(i,7)+_Parameters(116,0)*pow(_Training(i,7),2);
                            
                        _Jacobian(6,4) = _Parameters(5,0)+_Parameters(18,0)*_Training(i,1)+_Parameters(19,0)*_Training(i,2)+_Parameters(20,0)*_Training(i,3)+_Parameters(21,0)*_Training(i,4)+2*_Parameters(22,0)*_Training(i,5)+_Parameters(27,0)*_Training(i,6)+_Parameters(33,0)*_Training(i,7)+_Parameters(56,0)*pow(_Training(i,1),2)+_Parameters(57,0)*pow(_Training(i,2),2)+_Parameters(58,0)*pow(_Training(i,3),2)+_Parameters(59,0)*pow(_Training(i,4),2)+_Parameters(60,0)*_Training(i,1)*_Training(i,2)+_Parameters(61,0)*_Training(i,1)*_Training(i,3)+_Parameters(62,0)*_Training(i,1)*_Training(i,4)+_Parameters(63,0)*_Training(i,2)*_Training(i,3)+_Parameters(64,0)*_Training(i,2)*_Training(i,4)+_Parameters(65,0)*_Training(i,3)*_Training(i,4)+2*_Parameters(66,0)*_Training(i,1)*_Training(i,5)+2*_Parameters(67,0)*_Training(i,2)*_Training(i,5)+2*_Parameters(68,0)*_Training(i,3)*_Training(i,5)+2*_Parameters(69,0)*_Training(i,4)*_Training(i,5)+3*_Parameters(70,0)*pow(_Training(i,5),2)+2*_Parameters(75,0)*_Training(i,5)*_Training(i,6)+_Parameters(79,0)*_Training(i,1)*_Training(i,6)+_Parameters(82,0)*_Training(i,2)*_Training(i,6)+_Parameters(84,0)*_Training(i,3)*_Training(i,6)+_Parameters(85,0)*_Training(i,4)*_Training(i,6)+_Parameters(90,0)*pow(_Training(i,6),2)+2*_Parameters(96,0)*_Training(i,5)*_Training(i,7)+_Parameters(101,0)*_Training(i,1)*_Training(i,7)+_Parameters(105,0)*_Training(i,2)*_Training(i,7)+_Parameters(108,0)*_Training(i,3)*_Training(i,7)+_Parameters(110,0)*_Training(i,4)*_Training(i,7)+_Parameters(112,0)*_Training(i,6)*_Training(i,7)+_Parameters(117,0)*pow(_Training(i,7),2);
                            
                        _Jacobian(6,5) = _Parameters(6,0)+_Parameters(23,0)*_Training(i,1)+_Parameters(24,0)*_Training(i,2)+_Parameters(25,0)*_Training(i,3)+_Parameters(26,0)*_Training(i,4)+_Parameters(27,0)*_Training(i,5)+2*_Parameters(28,0)*_Training(i,6)+_Parameters(34,0)*_Training(i,7)+_Parameters(71,0)*pow(_Training(i,1),2)+_Parameters(72,0)*pow(_Training(i,2),2)+_Parameters(73,0)*pow(_Training(i,3),2)+_Parameters(74,0)*pow(_Training(i,4),2)+_Parameters(75,0)*pow(_Training(i,5),2)+_Parameters(76,0)*_Training(i,1)*_Training(i,2)+_Parameters(77,0)*_Training(i,1)*_Training(i,3)+_Parameters(78,0)*_Training(i,1)*_Training(i,4)+_Parameters(79,0)*_Training(i,1)*_Training(i,5)+_Parameters(80,0)*_Training(i,2)*_Training(i,3)+_Parameters(81,0)*_Training(i,2)*_Training(i,4)+_Parameters(82,0)*_Training(i,2)*_Training(i,5)+_Parameters(83,0)*_Training(i,3)*_Training(i,4)+_Parameters(84,0)*_Training(i,3)*_Training(i,5)+_Parameters(85,0)*_Training(i,4)*_Training(i,5)+2*_Parameters(86,0)*_Training(i,1)*_Training(i,6)+2*_Parameters(87,0)*_Training(i,2)*_Training(i,6)+2*_Parameters(88,0)*_Training(i,3)*_Training(i,6)+2*_Parameters(89,0)*_Training(i,4)*_Training(i,6)+2*_Parameters(90,0)*_Training(i,5)*_Training(i,6)+3*_Parameters(91,0)*pow(_Training(i,6),2)+2*_Parameters(97,0)*_Training(i,6)*_Training(i,7)+_Parameters(102,0)*_Training(i,1)*_Training(i,7)+_Parameters(106,0)*_Training(i,2)*_Training(i,7)+_Parameters(109,0)*_Training(i,3)*_Training(i,7)+_Parameters(111,0)*_Training(i,4)*_Training(i,7)+_Parameters(112,0)*_Training(i,5)*_Training(i,7)+_Parameters(118,0)*pow(_Training(i,7),2);
                            
                        _Jacobian(6,6) = _Parameters(7,0)+_Parameters(29,0)*_Training(i,1)+_Parameters(30,0)*_Training(i,2)+_Parameters(31,0)*_Training(i,3)+_Parameters(32,0)*_Training(i,4)+_Parameters(33,0)*_Training(i,5)+_Parameters(34,0)*_Training(i,6)+2*_Parameters(35,0)*_Training(i,7)+_Parameters(92,0)*pow(_Training(i,1),2)+_Parameters(93,0)*pow(_Training(i,2),2)+_Parameters(94,0)*pow(_Training(i,3),2)+_Parameters(95,0)*pow(_Training(i,4),2)+_Parameters(96,0)*pow(_Training(i,5),2)+_Parameters(97,0)*pow(_Training(i,6),2)+_Parameters(98,0)*_Training(i,1)*_Training(i,2)+_Parameters(99,0)*_Training(i,1)*_Training(i,3)+_Parameters(100,0)*_Training(i,1)*_Training(i,4)+_Parameters(101,0)*_Training(i,1)*_Training(i,5)+_Parameters(102,0)*_Training(i,1)*_Training(i,6)+_Parameters(103,0)*_Training(i,2)*_Training(i,3)+_Parameters(104,0)*_Training(i,2)*_Training(i,4)+_Parameters(105,0)*_Training(i,2)*_Training(i,5)+_Parameters(106,0)*_Training(i,2)*_Training(i,6)+_Parameters(107,0)*_Training(i,3)*_Training(i,4)+_Parameters(108,0)*_Training(i,3)*_Training(i,5)+_Parameters(109,0)*_Training(i,3)*_Training(i,6)+_Parameters(110,0)*_Training(i,4)*_Training(i,5)+_Parameters(111,0)*_Training(i,4)*_Training(i,6)+_Parameters(112,0)*_Training(i,5)*_Training(i,6)+2*_Parameters(113,0)*_Training(i,1)*_Training(i,7)+2*_Parameters(114,0)*_Training(i,2)*_Training(i,7)+2*_Parameters(115,0)*_Training(i,3)*_Training(i,7)+2*_Parameters(116,0)*_Training(i,4)*_Training(i,7)+2*_Parameters(117,0)*_Training(i,5)*_Training(i,7)+2*_Parameters(118,0)*_Training(i,6)*_Training(i,7)+3*_Parameters(119,0)*pow(_Training(i,7),2);
                            
                        break;
                    }
                    break;
                    
                case 4:
                    switch (d)
                    {
                        case 1:
                            _Jacobian(0,0) = _Parameters(1,0)+2*_Parameters(2,0)*_Training(i,1)+3*_Parameters(3,0)*pow(_Training(i,1),2)+4*_Parameters(4,0)*pow(_Training(i,1),3);
                            break;
                        case 2:
                            _Jacobian(1,0) = _Parameters(1,0)+2*_Parameters(3,0)*_Training(i,1)+_Parameters(4,0)*_Training(i,2)+3*_Parameters(6,0)*pow(_Training(i,1),2)+2*_Parameters(7,0)*_Training(i,1)*_Training(i,2)+_Parameters(8,0)*pow(_Training(i,2),2)+4*_Parameters(10,0)*pow(_Training(i,1),3)+3*_Parameters(11,0)*pow(_Training(i,1),2)*_Training(i,2)+2*_Parameters(12,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(13,0)*pow(_Training(i,2),3);
                            
                            _Jacobian(1,1) = _Parameters(2,0)+2*_Parameters(5,0)*_Training(i,5)+_Parameters(4,0)*_Training(i,1)+3*_Parameters(9,0)*pow(_Training(i,2),2)+2*_Parameters(8,0)*_Training(i,1)*_Training(i,2)+_Parameters(7,0)*pow(_Training(i,1),2)+4*_Parameters(14,0)*pow(_Training(i,2),3)+3*_Parameters(13,0)*pow(_Training(i,2),2)*_Training(i,1)+2*_Parameters(12,0)*_Training(i,2)*pow(_Training(i,1),2)+_Parameters(11,0)*pow(_Training(i,1),3);
                            break;
                    case 3:
                        
                            _Jacobian(2,0) = _Parameters(1,0)+2*_Parameters(4,0)*_Training(i,1)+_Parameters(5,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,3)+3*_Parameters(10,0)*pow(_Training(i,1),2)+2*_Parameters(11,0)*_Training(i,1)*_Training(i,2)+_Parameters(12,0)*pow(_Training(i,2),2)+2*_Parameters(14,0)*_Training(i,1)*_Training(i,3)+_Parameters(16,0)*_Training(i,2)*_Training(i,3)+_Parameters(17,0)*pow(_Training(i,3),2)+4*_Parameters(20,0)*pow(_Training(i,1),3)+3*_Parameters(21,0)*pow(_Training(i,1),2)*_Training(i,2)+2*_Parameters(22,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(23,0)*pow(_Training(i,2),3)+3*_Parameters(25,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(27,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(28,0)*pow(_Training(i,2),2)*_Training(i,3)+2*_Parameters(29,0)*_Training(i,1)*pow(_Training(i,3),2)+_Parameters(30,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(32,0)*pow(_Training(i,3),3);

                        
                        _Jacobian(2,1) = _Parameters(2,0)+2*_Parameters(6,0)*_Training(i,2)+_Parameters(5,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,3)+3*_Parameters(13,0)*pow(_Training(i,2),2)+2*_Parameters(12,0)*_Training(i,1)*_Training(i,2)+_Parameters(11,0)*pow(_Training(i,1),2)+2*_Parameters(15,0)*_Training(i,2)*_Training(i,3)+_Parameters(16,0)*_Training(i,1)*_Training(i,3)+_Parameters(18,0)*pow(_Training(i,3),2)+3*_Parameters(21,0)*pow(_Training(i,1),3)+2*_Parameters(22,0)*pow(_Training(i,1),2)*_Training(i,2)+3*_Parameters(23,0)*_Training(i,1)*pow(_Training(i,2),2)+4*_Parameters(24,0)*pow(_Training(i,2),3)+3*_Parameters(26,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(27,0)*_Training(i,1)*_Training(i,3)+2*_Parameters(28,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(30,0)*_Training(i,1)*pow(_Training(i,3),2)+2*_Parameters(31,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(33,0)*pow(_Training(i,3),3);
                            

                        _Jacobian(2,2) = _Parameters(3,0)+2*_Parameters(9,0)*_Training(i,3)+_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2)+3*_Parameters(19,0)*pow(_Training(i,3),2)+2*_Parameters(17,0)*_Training(i,1)*_Training(i,3)+_Parameters(14,0)*pow(_Training(i,1),2)+2*_Parameters(18,0)*_Training(i,2)*_Training(i,3)+_Parameters(16,0)*_Training(i,1)*_Training(i,2)+_Parameters(15,0)*pow(_Training(i,2),2)+_Parameters(25,0)*pow(_Training(i,1),3)+_Parameters(26,0)*pow(_Training(i,2),3)+_Parameters(27,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(28,0)*_Training(i,1)*pow(_Training(i,2),2)+2*_Parameters(29,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(30,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+2*_Parameters(31,0)*pow(_Training(i,2),2)*_Training(i,3)+3*_Parameters(32,0)*_Training(i,1)*pow(_Training(i,3),2)+3*_Parameters(34,0)*_Training(i,2)*pow(_Training(i,3),2)+4*_Parameters(33,0)*pow(_Training(i,3),3);
                            
                        break;
                    
                        case 4:
                        
                        _Jacobian(3,0) = _Parameters(1,0)+2*_Parameters(5,0)*_Training(i,1)+_Parameters(6,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,3)+_Parameters(11,0)*_Training(i,4)+3*_Parameters(15,0)*pow(_Training(i,1),2)+2*_Parameters(16,0)*_Training(i,1)*_Training(i,2)+_Parameters(17,0)*pow(_Training(i,2),2)+2*_Parameters(19,0)*_Training(i,1)*_Training(i,3)+_Parameters(21,0)*_Training(i,2)*_Training(i,3)+_Parameters(22,0)*pow(_Training(i,3),2)+2*_Parameters(25,0)*_Training(i,1)*_Training(i,4)+_Parameters(28,0)*_Training(i,2)*_Training(i,4)+_Parameters(29,0)*_Training(i,3)*_Training(i,4)+_Parameters(31,0)*pow(_Training(i,4),2)+4*_Parameters(35,0)*pow(_Training(i,1),3)+3*_Parameters(36,0)*pow(_Training(i,1),2)*_Training(i,2)+2*_Parameters(37,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(38,0)*pow(_Training(i,2),3)+3*_Parameters(40,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(42,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(43,0)*pow(_Training(i,2),2)*_Training(i,3)+2*_Parameters(44,0)*_Training(i,1)*pow(_Training(i,3),2)+_Parameters(45,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(47,0)*pow(_Training(i,3),3)+3*_Parameters(50,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(51,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(52,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+2*_Parameters(53,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(54,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(55,0)*pow(_Training(i,3),2)*_Training(i,4)+2*_Parameters(60,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(61,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(62,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(66,0)*pow(_Training(i,4),3);
                        
                        _Jacobian(3,1) = _Parameters(2,0)+2*_Parameters(7,0)*_Training(i,2)+_Parameters(6,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4)+3*_Parameters(18,0)*pow(_Training(i,2),2)+2*_Parameters(17,0)*_Training(i,1)*_Training(i,2)+_Parameters(16,0)*pow(_Training(i,1),2)+2*_Parameters(20,0)*_Training(i,2)*_Training(i,3)+_Parameters(21,0)*_Training(i,1)*_Training(i,3)+_Parameters(23,0)*pow(_Training(i,3),2)+2*_Parameters(26,0)*_Training(i,2)*_Training(i,4)+_Parameters(28,0)*_Training(i,1)*_Training(i,4)+_Parameters(30,0)*_Training(i,3)*_Training(i,4)+_Parameters(32,0)*pow(_Training(i,4),2)+3*_Parameters(36,0)*pow(_Training(i,1),3)+2*_Parameters(37,0)*pow(_Training(i,1),2)*_Training(i,2)+3*_Parameters(38,0)*_Training(i,1)*pow(_Training(i,2),2)+4*_Parameters(39,0)*pow(_Training(i,2),3)+3*_Parameters(41,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(42,0)*_Training(i,1)*_Training(i,3)+2*_Parameters(43,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(45,0)*_Training(i,1)*pow(_Training(i,3),2)+2*_Parameters(46,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(48,0)*pow(_Training(i,3),3)+_Parameters(51,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(53,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+_Parameters(54,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+3*_Parameters(56,0)*pow(_Training(i,2),2)*_Training(i,4)+2*_Parameters(57,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(58,0)*pow(_Training(i,3),2)*_Training(i,4)+_Parameters(61,0)*_Training(i,1)*pow(_Training(i,4),2)+2*_Parameters(63,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(64,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(67,0)*pow(_Training(i,4),3);


                        
                        _Jacobian(3,2) = _Parameters(3,0)+2*_Parameters(10,0)*_Training(i,3)+_Parameters(9,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,1)+_Parameters(13,0)*_Training(i,4)+3*_Parameters(24,0)*pow(_Training(i,3),2)+2*_Parameters(22,0)*_Training(i,1)*_Training(i,3)+_Parameters(19,0)*pow(_Training(i,1),2)+2*_Parameters(23,0)*_Training(i,2)*_Training(i,3)+_Parameters(21,0)*_Training(i,1)*_Training(i,2)+_Parameters(20,0)*pow(_Training(i,2),2)+2*_Parameters(27,0)*_Training(i,3)*_Training(i,4)+_Parameters(29,0)*_Training(i,1)*_Training(i,4)+_Parameters(30,0)*_Training(i,2)*_Training(i,4)+_Parameters(33,0)*pow(_Training(i,4),2)+_Parameters(40,0)*pow(_Training(i,1),3)+_Parameters(41,0)*pow(_Training(i,2),3)+_Parameters(42,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(43,0)*_Training(i,1)*pow(_Training(i,2),2)+2*_Parameters(44,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(45,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+2*_Parameters(46,0)*pow(_Training(i,2),2)*_Training(i,3)+3*_Parameters(47,0)*_Training(i,1)*pow(_Training(i,3),2)+3*_Parameters(48,0)*_Training(i,2)*pow(_Training(i,3),2)+4*_Parameters(49,0)*pow(_Training(i,3),3)+_Parameters(52,0)*pow(_Training(i,1),2)*_Training(i,4)+_Parameters(54,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(55,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+_Parameters(57,0)*pow(_Training(i,2),2)*_Training(i,4)+2*_Parameters(58,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+3*_Parameters(59,0)*pow(_Training(i,3),2)*_Training(i,4)+_Parameters(62,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(64,0)*_Training(i,2)*pow(_Training(i,4),2)+2*_Parameters(65,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(68,0)*pow(_Training(i,4),3);
                            

                        
                        _Jacobian(3,3) = _Parameters(4,0)+_Parameters(11,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,2)+_Parameters(13,0)*_Training(i,3)+2*_Parameters(14,0)*_Training(i,4)+_Parameters(25,0)*pow(_Training(i,1),2)+_Parameters(26,0)*pow(_Training(i,2),2)+_Parameters(27,0)*pow(_Training(i,3),2)+_Parameters(28,0)*_Training(i,1)*_Training(i,2)+_Parameters(29,0)*_Training(i,1)*_Training(i,3)+_Parameters(30,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(31,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(32,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(33,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(34,0)*pow(_Training(i,4),2)+_Parameters(50,0)*pow(_Training(i,1),3)+_Parameters(51,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(52,0)*pow(_Training(i,1),2)*_Training(i,3)+_Parameters(53,0)*pow(_Training(i,2),2)*_Training(i,1)+_Parameters(54,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(55,0)*pow(_Training(i,3),2)*_Training(i,1)+_Parameters(56,0)*pow(_Training(i,3),3)+_Parameters(57,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(58,0)*pow(_Training(i,3),2)*_Training(i,2)+_Parameters(59,0)*pow(_Training(i,3),3)+2*_Parameters(60,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(61,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(62,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+_Parameters(63,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(64,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+2*_Parameters(65,0)*pow(_Training(i,3),2)*_Training(i,4)+3*_Parameters(66,0)*pow(_Training(i,4),2)*_Training(i,1)+3*_Parameters(67,0)*pow(_Training(i,4),2)*_Training(i,2)+3*_Parameters(68,0)*pow(_Training(i,4),2)*_Training(i,3)+4*_Parameters(69,0)*pow(_Training(i,4),3);
                        
                        break;
                    
                    case 5:
                        
                        _Jacobian(4,0) = _Parameters(1,0)+2*_Parameters(6,0)*_Training(i,1)+_Parameters(7,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,4)+_Parameters(16,0)*_Training(i,5)+3*_Parameters(21,0)*pow(_Training(i,1),2)+2*_Parameters(22,0)*_Training(i,1)*_Training(i,2)+_Parameters(23,0)*pow(_Training(i,2),2)+2*_Parameters(25,0)*_Training(i,1)*_Training(i,3)+_Parameters(27,0)*_Training(i,2)*_Training(i,3)+_Parameters(28,0)*pow(_Training(i,3),2)+2*_Parameters(31,0)*_Training(i,1)*_Training(i,4)+_Parameters(34,0)*_Training(i,2)*_Training(i,4)+_Parameters(35,0)*_Training(i,3)*_Training(i,4)+_Parameters(37,0)*pow(_Training(i,4),2)+2*_Parameters(41,0)*_Training(i,1)*_Training(i,5)+_Parameters(45,0)*_Training(i,2)*_Training(i,5)+_Parameters(46,0)*_Training(i,3)*_Training(i,5)+_Parameters(47,0)*_Training(i,4)*_Training(i,5)+_Parameters(51,0)*pow(_Training(i,5),2)+4*_Parameters(56,0)*pow(_Training(i,1),3)+3*_Parameters(57,0)*pow(_Training(i,1),2)*_Training(i,2)+2*_Parameters(58,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(59,0)*pow(_Training(i,2),3)+3*_Parameters(61,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(63,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(64,0)*pow(_Training(i,2),2)*_Training(i,3)+2*_Parameters(65,0)*_Training(i,1)*pow(_Training(i,3),2)+_Parameters(66,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(68,0)*pow(_Training(i,3),3)+3*_Parameters(71,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(72,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(73,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+2*_Parameters(74,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(75,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(76,0)*pow(_Training(i,3),2)*_Training(i,4)+2*_Parameters(81,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(82,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(83,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(87,0)*pow(_Training(i,4),3)+3*_Parameters(91,0)*pow(_Training(i,1),2)*_Training(i,5)+2*_Parameters(92,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+2*_Parameters(93,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+2*_Parameters(94,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+_Parameters(95,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+_Parameters(96,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+_Parameters(97,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+_Parameters(98,0)*pow(_Training(i,2),2)*_Training(i,5)+_Parameters(99,0)*pow(_Training(i,3),2)*_Training(i,5)+_Parameters(100,0)*pow(_Training(i,4),2)*_Training(i,5)+2*_Parameters(111,0)*_Training(i,1)*pow(_Training(i,5),2)+_Parameters(112,0)*_Training(i,2)*pow(_Training(i,5),2)+_Parameters(113,0)*_Training(i,3)*pow(_Training(i,5),2)+_Parameters(114,0)*_Training(i,4)*pow(_Training(i,5),2)+_Parameters(121,0)*pow(_Training(i,5),3);
                        
                        _Jacobian(4,1) = _Parameters(2,0)+2*_Parameters(8,0)*_Training(i,2)+_Parameters(7,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5)+3*_Parameters(24,0)*pow(_Training(i,2),2)+2*_Parameters(23,0)*_Training(i,1)*_Training(i,2)+_Parameters(22,0)*pow(_Training(i,1),2)+2*_Parameters(26,0)*_Training(i,2)*_Training(i,3)+_Parameters(27,0)*_Training(i,1)*_Training(i,3)+_Parameters(29,0)*pow(_Training(i,3),2)+2*_Parameters(32,0)*_Training(i,2)*_Training(i,4)+_Parameters(34,0)*_Training(i,1)*_Training(i,4)+_Parameters(36,0)*_Training(i,3)*_Training(i,4)+_Parameters(38,0)*pow(_Training(i,4),2)+2*_Parameters(42,0)*_Training(i,2)*_Training(i,5)+_Parameters(45,0)*_Training(i,1)*_Training(i,5)+_Parameters(48,0)*_Training(i,3)*_Training(i,5)+_Parameters(49,0)*_Training(i,4)*_Training(i,5)+_Parameters(52,0)*pow(_Training(i,5),2)+3*_Parameters(57,0)*pow(_Training(i,1),3)+2*_Parameters(58,0)*pow(_Training(i,1),2)*_Training(i,2)+3*_Parameters(59,0)*_Training(i,1)*pow(_Training(i,2),2)+4*_Parameters(60,0)*pow(_Training(i,2),3)+3*_Parameters(62,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(63,0)*_Training(i,1)*_Training(i,3)+2*_Parameters(64,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(66,0)*_Training(i,1)*pow(_Training(i,3),2)+2*_Parameters(67,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(69,0)*pow(_Training(i,3),3)+_Parameters(72,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(74,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+_Parameters(75,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+3*_Parameters(77,0)*pow(_Training(i,2),2)*_Training(i,4)+2*_Parameters(78,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(79,0)*pow(_Training(i,3),2)*_Training(i,4)+_Parameters(82,0)*_Training(i,1)*pow(_Training(i,4),2)+2*_Parameters(84,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(85,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(88,0)*pow(_Training(i,4),3)+_Parameters(92,0)*pow(_Training(i,1),2)*_Training(i,5)+2*_Parameters(99,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+_Parameters(95,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+_Parameters(96,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+2*_Parameters(102,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+2*_Parameters(103,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+_Parameters(104,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+3*_Parameters(101,0)*pow(_Training(i,2),2)*_Training(i,5)+_Parameters(105,0)*pow(_Training(i,3),2)*_Training(i,5)+_Parameters(106,0)*pow(_Training(i,4),2)*_Training(i,5)+_Parameters(114,0)*_Training(i,1)*pow(_Training(i,5),2)+2*_Parameters(115,0)*_Training(i,2)*pow(_Training(i,5),2)+_Parameters(116,0)*_Training(i,3)*pow(_Training(i,5),2)+_Parameters(117,0)*_Training(i,4)*pow(_Training(i,5),2)+_Parameters(122,0)*pow(_Training(i,5),3);

                            
                        _Jacobian(4,2) = _Parameters(3,0)+2*_Parameters(11,0)*_Training(i,3)+_Parameters(9,0)*_Training(i,1)+_Parameters(10,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+3*_Parameters(30,0)*pow(_Training(i,3),2)+2*_Parameters(28,0)*_Training(i,1)*_Training(i,3)+_Parameters(25,0)*pow(_Training(i,1),2)+2*_Parameters(29,0)*_Training(i,2)*_Training(i,3)+_Parameters(27,0)*_Training(i,1)*_Training(i,2)+_Parameters(26,0)*pow(_Training(i,2),2)+2*_Parameters(33,0)*_Training(i,3)*_Training(i,4)+_Parameters(35,0)*_Training(i,1)*_Training(i,4)+_Parameters(36,0)*_Training(i,2)*_Training(i,4)+_Parameters(39,0)*pow(_Training(i,4),2)+2*_Parameters(43,0)*_Training(i,3)*_Training(i,5)+_Parameters(46,0)*_Training(i,1)*_Training(i,5)+_Parameters(48,0)*_Training(i,2)*_Training(i,5)+_Parameters(50,0)*_Training(i,4)*_Training(i,5)+_Parameters(53,0)*pow(_Training(i,5),2)+_Parameters(61,0)*pow(_Training(i,1),3)+_Parameters(62,0)*pow(_Training(i,2),3)+_Parameters(63,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(64,0)*_Training(i,1)*pow(_Training(i,2),2)+2*_Parameters(65,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(66,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+2*_Parameters(67,0)*pow(_Training(i,2),2)*_Training(i,3)+3*_Parameters(68,0)*_Training(i,1)*pow(_Training(i,3),2)+3*_Parameters(69,0)*_Training(i,2)*pow(_Training(i,3),2)+4*_Parameters(70,0)*pow(_Training(i,3),3)+_Parameters(73,0)*pow(_Training(i,1),2)*_Training(i,4)+_Parameters(75,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(76,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+_Parameters(78,0)*pow(_Training(i,2),2)*_Training(i,4)+2*_Parameters(79,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+3*_Parameters(80,0)*pow(_Training(i,3),2)*_Training(i,4)+_Parameters(83,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(85,0)*_Training(i,2)*pow(_Training(i,4),2)+2*_Parameters(86,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(89,0)*pow(_Training(i,4),3)+_Parameters(93,0)*pow(_Training(i,1),2)*_Training(i,5)+_Parameters(95,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+2*_Parameters(99,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+_Parameters(97,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+2*_Parameters(105,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+_Parameters(104,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+2*_Parameters(108,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+_Parameters(102,0)*pow(_Training(i,2),2)*_Training(i,5)+_Parameters(107,0)*pow(_Training(i,3),2)*_Training(i,5)+_Parameters(109,0)*pow(_Training(i,4),2)*_Training(i,5)+_Parameters(113,0)*_Training(i,1)*pow(_Training(i,5),2)+_Parameters(116,0)*_Training(i,2)*pow(_Training(i,5),2)+2*_Parameters(118,0)*_Training(i,3)*pow(_Training(i,5),2)+_Parameters(119,0)*_Training(i,4)*pow(_Training(i,5),2)+_Parameters(123,0)*pow(_Training(i,5),3);
                            
                        _Jacobian(4,3) = _Parameters(4,0)+2*_Parameters(15,0)*_Training(i,4)+_Parameters(12,0)*_Training(i,1)+_Parameters(13,0)*_Training(i,2)+_Parameters(14,0)*_Training(i,3)+_Parameters(19,0)*_Training(i,5)+3*_Parameters(40,0)*pow(_Training(i,4),2)+2*_Parameters(37,0)*_Training(i,1)*_Training(i,4)+_Parameters(31,0)*pow(_Training(i,1),2)+2*_Parameters(38,0)*_Training(i,2)*_Training(i,4)+_Parameters(34,0)*_Training(i,1)*_Training(i,2)+_Parameters(32,0)*pow(_Training(i,2),2)+2*_Parameters(39,0)*_Training(i,3)*_Training(i,4)+_Parameters(35,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*_Training(i,2)*_Training(i,3)+_Parameters(33,0)*pow(_Training(i,3),2)+2*_Parameters(44,0)*_Training(i,4)*_Training(i,5)+_Parameters(47,0)*_Training(i,1)*_Training(i,5)+_Parameters(49,0)*_Training(i,2)*_Training(i,5)+_Parameters(50,0)*_Training(i,3)*_Training(i,5)+_Parameters(54,0)*pow(_Training(i,5),2)+_Parameters(71,0)*pow(_Training(i,1),3)+_Parameters(72,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(73,0)*pow(_Training(i,1),2)*_Training(i,3)+_Parameters(74,0)*pow(_Training(i,2),2)*_Training(i,1)+_Parameters(75,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(76,0)*pow(_Training(i,3),2)*_Training(i,1)+_Parameters(77,0)*pow(_Training(i,3),3)+_Parameters(78,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(79,0)*pow(_Training(i,3),2)*_Training(i,2)+_Parameters(80,0)*pow(_Training(i,3),3)+2*_Parameters(81,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(82,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(83,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+_Parameters(84,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(85,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+2*_Parameters(86,0)*pow(_Training(i,3),2)*_Training(i,4)+3*_Parameters(87,0)*pow(_Training(i,4),2)*_Training(i,1)+3*_Parameters(88,0)*pow(_Training(i,4),2)*_Training(i,2)+3*_Parameters(89,0)*pow(_Training(i,4),2)*_Training(i,3)+4*_Parameters(90,0)*pow(_Training(i,4),3)+_Parameters(94,0)*pow(_Training(i,1),2)*_Training(i,5)+_Parameters(96,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+_Parameters(97,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+2*_Parameters(100,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+_Parameters(104,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+2*_Parameters(106,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+2*_Parameters(109,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+_Parameters(103,0)*pow(_Training(i,2),2)*_Training(i,5)+_Parameters(108,0)*pow(_Training(i,3),2)*_Training(i,5)+3*_Parameters(110,0)*pow(_Training(i,4),2)*_Training(i,5)+_Parameters(114,0)*_Training(i,1)*pow(_Training(i,5),2)+_Parameters(117,0)*_Training(i,2)*pow(_Training(i,5),2)+_Parameters(119,0)*_Training(i,3)*pow(_Training(i,5),2)+2*_Parameters(120,0)*_Training(i,4)*pow(_Training(i,5),2)+_Parameters(124,0)*pow(_Training(i,5),3);
                            
                            
                        _Jacobian(4,4) = _Parameters(5,0)+2*_Parameters(20,0)*_Training(i,5)+_Parameters(16,0)*_Training(i,1)+_Parameters(17,0)*_Training(i,2)+_Parameters(18,0)*_Training(i,3)+_Parameters(19,0)*_Training(i,4)+3*_Parameters(55,0)*pow(_Training(i,5),2)+2*_Parameters(51,0)*_Training(i,1)*_Training(i,5)+_Parameters(41,0)*pow(_Training(i,1),2)+2*_Parameters(52,0)*_Training(i,2)*_Training(i,5)+_Parameters(45,0)*_Training(i,1)*_Training(i,2)+_Parameters(42,0)*pow(_Training(i,2),2)+2*_Parameters(53,0)*_Training(i,3)*_Training(i,5)+_Parameters(46,0)*_Training(i,1)*_Training(i,3)+_Parameters(48,0)*_Training(i,2)*_Training(i,3)+_Parameters(43,0)*pow(_Training(i,3),2)+2*_Parameters(54,0)*_Training(i,4)*_Training(i,5)+_Parameters(47,0)*_Training(i,1)*_Training(i,4)+_Parameters(49,0)*_Training(i,2)*_Training(i,4)+_Parameters(50,0)*_Training(i,3)*_Training(i,4)+_Parameters(44,0)*pow(_Training(i,4),2)+_Parameters(91,0)*pow(_Training(i,1),3)+_Parameters(92,0)*pow(_Training(i,1),2)*_Training(i,2)+_Parameters(98,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(103,0)*pow(_Training(i,2),3)+_Parameters(93,0)*pow(_Training(i,1),2)*_Training(i,3)+_Parameters(95,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(102,0)*pow(_Training(i,2),2)*_Training(i,3)+_Parameters(99,0)*_Training(i,1)*pow(_Training(i,3),2)+_Parameters(105,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(107,0)*pow(_Training(i,3),3)+_Parameters(94,0)*pow(_Training(i,1),2)*_Training(i,4)+_Parameters(96,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+_Parameters(97,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+_Parameters(103,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(104,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(108,0)*pow(_Training(i,3),2)*_Training(i,4)+_Parameters(100,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(106,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(109,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(110,0)*pow(_Training(i,4),3)+2*_Parameters(111,0)*pow(_Training(i,1),2)*_Training(i,5)+2*_Parameters(112,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+2*_Parameters(113,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+2*_Parameters(114,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+2*_Parameters(116,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+2*_Parameters(117,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+2*_Parameters(119,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+2*_Parameters(115,0)*pow(_Training(i,2),2)*_Training(i,5)+2*_Parameters(118,0)*pow(_Training(i,3),2)*_Training(i,5)+2*_Parameters(120,0)*pow(_Training(i,4),2)*_Training(i,5)+3*_Parameters(121,0)*_Training(i,1)*pow(_Training(i,5),2)+3*_Parameters(132,0)*_Training(i,2)*pow(_Training(i,5),2)+3*_Parameters(123,0)*_Training(i,3)*pow(_Training(i,5),2)+3*_Parameters(124,0)*_Training(i,4)*pow(_Training(i,5),2)+4*_Parameters(125,0)*pow(_Training(i,5),3);
                            

                        
                        break;
                    
                    case 6:
                    
                        _Jacobian(5,0) = _Parameters(1,0)+2*_Parameters(7,0)*_Training(i,1)+_Parameters(8,0)*_Training(i,2)+_Parameters(10,0)*_Training(i,3)+_Parameters(13,0)*_Training(i,4)+_Parameters(17,0)*_Training(i,5)+_Parameters(22,0)*_Training(i,6)+3*_Parameters(28,0)*pow(_Training(i,1),2)+2*_Parameters(29,0)*_Training(i,1)*_Training(i,2)+_Parameters(30,0)*pow(_Training(i,2),2)+2*_Parameters(32,0)*_Training(i,1)*_Training(i,2)+_Parameters(34,0)*_Training(i,2)*_Training(i,3)+_Parameters(35,0)*pow(_Training(i,3),2)+2*_Parameters(38,0)*_Training(i,1)*_Training(i,4)+_Parameters(41,0)*_Training(i,2)*_Training(i,4)+_Parameters(42,0)*_Training(i,3)*_Training(i,4)+_Parameters(44,0)*pow(_Training(i,4),2)+2*_Parameters(48,0)*_Training(i,1)*_Training(i,5)+_Parameters(52,0)*_Training(i,2)*_Training(i,5)+_Parameters(53,0)*_Training(i,3)*_Training(i,5)+_Parameters(54,0)*_Training(i,4)*_Training(i,5)+_Parameters(58,0)*pow(_Training(i,5),2)+2*_Parameters(63,0)*_Training(i,1)*_Training(i,6)+_Parameters(68,0)*_Training(i,2)*_Training(i,6)+_Parameters(69,0)*_Training(i,3)*_Training(i,6)+_Parameters(70,0)*_Training(i,4)*_Training(i,6)+_Parameters(71,0)*_Training(i,5)*_Training(i,6)+_Parameters(78,0)*pow(_Training(i,6),2)
                            
                            +4*_Parameters(56,0)*pow(_Training(i,1),3)+3*_Parameters(57,0)*pow(_Training(i,1),2)*_Training(i,2)+2*_Parameters(58,0)*_Training(i,1)*pow(_Training(i,2),2)+_Parameters(59,0)*pow(_Training(i,2),3)+3*_Parameters(61,0)*pow(_Training(i,1),2)*_Training(i,3)+2*_Parameters(63,0)*_Training(i,1)*_Training(i,2)*_Training(i,3)+_Parameters(64,0)*pow(_Training(i,2),2)*_Training(i,3)+2*_Parameters(65,0)*_Training(i,1)*pow(_Training(i,3),2)+_Parameters(66,0)*_Training(i,2)*pow(_Training(i,3),2)+_Parameters(68,0)*pow(_Training(i,3),3)+3*_Parameters(71,0)*pow(_Training(i,1),2)*_Training(i,4)+2*_Parameters(72,0)*_Training(i,1)*_Training(i,2)*_Training(i,4)+2*_Parameters(73,0)*_Training(i,1)*_Training(i,3)*_Training(i,4)+2*_Parameters(74,0)*pow(_Training(i,2),2)*_Training(i,4)+_Parameters(75,0)*_Training(i,2)*_Training(i,3)*_Training(i,4)+_Parameters(76,0)*pow(_Training(i,3),2)*_Training(i,4)+2*_Parameters(81,0)*_Training(i,1)*pow(_Training(i,4),2)+_Parameters(82,0)*_Training(i,2)*pow(_Training(i,4),2)+_Parameters(83,0)*_Training(i,3)*pow(_Training(i,4),2)+_Parameters(87,0)*pow(_Training(i,4),3)+3*_Parameters(91,0)*pow(_Training(i,1),2)*_Training(i,5)+2*_Parameters(92,0)*_Training(i,1)*_Training(i,2)*_Training(i,5)+2*_Parameters(93,0)*_Training(i,1)*_Training(i,3)*_Training(i,5)+2*_Parameters(94,0)*_Training(i,1)*_Training(i,4)*_Training(i,5)+_Parameters(95,0)*_Training(i,2)*_Training(i,3)*_Training(i,5)+_Parameters(96,0)*_Training(i,2)*_Training(i,4)*_Training(i,5)+_Parameters(97,0)*_Training(i,3)*_Training(i,4)*_Training(i,5)+_Parameters(98,0)*pow(_Training(i,2),2)*_Training(i,5)+_Parameters(99,0)*pow(_Training(i,3),2)*_Training(i,5)+_Parameters(100,0)*pow(_Training(i,4),2)*_Training(i,5)+2*_Parameters(111,0)*_Training(i,1)*pow(_Training(i,5),2)+_Parameters(112,0)*_Training(i,2)*pow(_Training(i,5),2)+_Parameters(113,0)*_Training(i,3)*pow(_Training(i,5),2)+_Parameters(114,0)*_Training(i,4)*pow(_Training(i,5),2)+_Parameters(121,0)*pow(_Training(i,5),3)
                            
                            
                            
                            
                            
                            
                            
                            ;

                            
                            
                            
                            
                            
                        
                        _Jacobian(5,1) = _Parameters(2,0)+2*_Parameters(9,0)*_Training(i,2)+_Parameters(8,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6)+3*_Parameters(31,0)*pow(_Training(i,2),2)+2*_Parameters(30,0)*_Training(i,1)*_Training(i,2)+_Parameters(29,0)*pow(_Training(i,1),2)+2*_Parameters(33,0)*_Training(i,2)*_Training(i,3)+_Parameters(34,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*pow(_Training(i,3),2)+2*_Parameters(39,0)*_Training(i,2)*_Training(i,4)+_Parameters(41,0)*_Training(i,1)*_Training(i,4)+_Parameters(43,0)*_Training(i,3)*_Training(i,4)+_Parameters(45,0)*pow(_Training(i,4),2)+2*_Parameters(49,0)*_Training(i,2)*_Training(i,5)+_Parameters(52,0)*_Training(i,1)*_Training(i,5)+_Parameters(55,0)*_Training(i,3)*_Training(i,5)+_Parameters(56,0)*_Training(i,4)*_Training(i,5)+_Parameters(59,0)*pow(_Training(i,5),2)+2*_Parameters(64,0)*_Training(i,2)*_Training(i,6)+_Parameters(68,0)*_Training(i,1)*_Training(i,6)+_Parameters(72,0)*_Training(i,3)*_Training(i,6)+_Parameters(73,0)*_Training(i,4)*_Training(i,6)+_Parameters(74,0)*_Training(i,5)*_Training(i,6)+_Parameters(79,0)*pow(_Training(i,6),2);
                        
                        _Jacobian(5,2) = _Parameters(3,0)+_Parameters(10,0)*_Training(i,1)+_Parameters(11,0)*_Training(i,2)+2*_Parameters(12,0)*_Training(i,3)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6)+_Parameters(32,0)*pow(_Training(i,1),2)+_Parameters(33,0)*pow(_Training(i,2),2)+_Parameters(34,0)*_Training(i,1)*_Training(i,2)+2*_Parameters(35,0)*_Training(i,1)*_Training(i,3)+_Parameters(36,0)*_Training(i,2)*_Training(i,3)+3*_Parameters(37,0)*pow(_Training(i,3),2)+2*_Parameters(40,0)*_Training(i,3)*_Training(i,4)+_Parameters(42,0)*_Training(i,1)*_Training(i,4)+_Parameters(43,0)*_Training(i,2)*_Training(i,4)+_Parameters(46,0)*pow(_Training(i,4),2)+2*_Parameters(50,0)*_Training(i,3)*_Training(i,5)+_Parameters(53,0)*_Training(i,1)*_Training(i,5)+_Parameters(55,0)*_Training(i,2)*_Training(i,5)+_Parameters(57,0)*_Training(i,4)*_Training(i,5)+_Parameters(60,0)*pow(_Training(i,5),2)+2*_Parameters(65,0)*_Training(i,3)*_Training(i,6)+_Parameters(69,0)*_Training(i,1)*_Training(i,6)+_Parameters(72,0)*_Training(i,2)*_Training(i,6)+_Parameters(75,0)*_Training(i,4)*_Training(i,6)+_Parameters(76,0)*_Training(i,5)*_Training(i,6)+_Parameters(80,0)*pow(_Training(i,6),2);
                        
                        _Jacobian(5,3) = _Parameters(4,0)+_Parameters(13,0)*_Training(i,1)+_Parameters(14,0)*_Training(i,2)+_Parameters(15,0)*_Training(i,3)+2*_Parameters(16,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6)+_Parameters(38,0)*pow(_Training(i,1),2)+_Parameters(39,0)*pow(_Training(i,2),2)+_Parameters(40,0)*pow(_Training(i,3),2)+_Parameters(41,0)*_Training(i,1)*_Training(i,2)+_Parameters(42,0)*_Training(i,1)*_Training(i,3)+_Parameters(43,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(44,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(45,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(46,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(47,0)*pow(_Training(i,4),2)+2*_Parameters(51,0)*_Training(i,4)*_Training(i,5)+_Parameters(54,0)*_Training(i,1)*_Training(i,5)+_Parameters(56,0)*_Training(i,2)*_Training(i,5)+_Parameters(57,0)*_Training(i,3)*_Training(i,5)+_Parameters(61,0)*pow(_Training(i,5),2)+2*_Parameters(66,0)*_Training(i,4)*_Training(i,6)+_Parameters(70,0)*_Training(i,1)*_Training(i,6)+_Parameters(73,0)*_Training(i,2)*_Training(i,6)+_Parameters(75,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,5)*_Training(i,6)+_Parameters(81,0)*pow(_Training(i,6),2);
                        
                        _Jacobian(5,4) = _Parameters(5,0)+_Parameters(17,0)*_Training(i,1)+_Parameters(18,0)*_Training(i,2)+_Parameters(19,0)*_Training(i,3)+2*_Parameters(21,0)*_Training(i,5)+_Parameters(20,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,6)+_Parameters(48,0)*pow(_Training(i,1),2)+_Parameters(49,0)*pow(_Training(i,2),2)+_Parameters(50,0)*pow(_Training(i,3),2)+_Parameters(51,0)*pow(_Training(i,4),2)+_Parameters(52,0)*_Training(i,1)*_Training(i,2)+_Parameters(53,0)*_Training(i,1)*_Training(i,3)+_Parameters(54,0)*_Training(i,1)*_Training(i,4)+_Parameters(55,0)*_Training(i,2)*_Training(i,3)+_Parameters(56,0)*_Training(i,2)*_Training(i,4)+_Parameters(57,0)*_Training(i,3)*_Training(i,4)+2*_Parameters(58,0)*_Training(i,1)*_Training(i,5)+2*_Parameters(59,0)*_Training(i,2)*_Training(i,5)+2*_Parameters(60,0)*_Training(i,3)*_Training(i,5)+2*_Parameters(61,0)*_Training(i,4)*_Training(i,5)+3*_Parameters(62,0)*pow(_Training(i,5),2)+2*_Parameters(67,0)*_Training(i,5)*_Training(i,6)+_Parameters(71,0)*_Training(i,1)*_Training(i,6)+_Parameters(74,0)*_Training(i,2)*_Training(i,6)+_Parameters(76,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,4)*_Training(i,6)+_Parameters(82,0)*pow(_Training(i,6),2);
                        
                        _Jacobian(5,5) = _Parameters(6,0)+_Parameters(22,0)*_Training(i,1)+_Parameters(23,0)*_Training(i,2)+_Parameters(24,0)*_Training(i,3)+2*_Parameters(27,0)*_Training(i,6)+_Parameters(25,0)*_Training(i,4)+_Parameters(26,0)*_Training(i,5)+_Parameters(63,0)*pow(_Training(i,1),2)+_Parameters(64,0)*pow(_Training(i,2),2)+_Parameters(65,0)*pow(_Training(i,3),2)+_Parameters(66,0)*pow(_Training(i,4),2)+_Parameters(68,0)*_Training(i,1)*_Training(i,2)+_Parameters(69,0)*_Training(i,1)*_Training(i,3)+_Parameters(70,0)*_Training(i,1)*_Training(i,4)+_Parameters(71,0)*_Training(i,1)*_Training(i,5)+_Parameters(72,0)*_Training(i,2)*_Training(i,3)+_Parameters(73,0)*_Training(i,2)*_Training(i,4)+_Parameters(74,0)*_Training(i,2)*_Training(i,5)+_Parameters(75,0)*_Training(i,3)*_Training(i,4)+_Parameters(76,0)*_Training(i,3)*_Training(i,5)+_Parameters(77,0)*_Training(i,4)*_Training(i,5)+2*_Parameters(78,0)*_Training(i,1)*_Training(i,6)+2*_Parameters(79,0)*_Training(i,2)*_Training(i,6)+2*_Parameters(80,0)*_Training(i,3)*_Training(i,6)+2*_Parameters(81,0)*_Training(i,4)*_Training(i,6)+2*_Parameters(85,0)*_Training(i,5)*_Training(i,6)+2*_Parameters(67,0)*pow(_Training(i,5),2)+3*_Parameters(83,0)*pow(_Training(i,6),2);
                        
                            
                        break;
                            
                    case 7:
                            
                        _Jacobian(6,0) = _Parameters(1,0)+2*_Parameters(8,0)*_Training(i,1)+_Parameters(9,0)*_Training(i,2)+_Parameters(11,0)*_Training(i,3)+_Parameters(14,0)*_Training(i,4)+_Parameters(18,0)*_Training(i,5)+_Parameters(23,0)*_Training(i,6)+_Parameters(29,0)*_Training(i,7)+3*_Parameters(36,0)*pow(_Training(i,1),2)+2*_Parameters(37,0)*_Training(i,1)*_Training(i,2)+_Parameters(38,0)*pow(_Training(i,2),2)+2*_Parameters(40,0)*_Training(i,1)*_Training(i,3)+_Parameters(42,0)*_Training(i,2)*_Training(i,3)+_Parameters(43,0)*pow(_Training(i,3),2)+2*_Parameters(46,0)*_Training(i,1)*_Training(i,4)+_Parameters(49,0)*_Training(i,2)*_Training(i,4)+_Parameters(50,0)*_Training(i,3)*_Training(i,4)+_Parameters(52,0)*pow(_Training(i,4),2)+2*_Parameters(56,0)*_Training(i,1)*_Training(i,5)+_Parameters(60,0)*_Training(i,2)*_Training(i,5)+_Parameters(61,0)*_Training(i,3)*_Training(i,5)+_Parameters(62,0)*_Training(i,4)*_Training(i,5)+_Parameters(66,0)*pow(_Training(i,5),2)+2*_Parameters(71,0)*_Training(i,1)*_Training(i,6)+_Parameters(76,0)*_Training(i,2)*_Training(i,6)+_Parameters(77,0)*_Training(i,3)*_Training(i,6)+_Parameters(78,0)*_Training(i,4)*_Training(i,6)+_Parameters(79,0)*_Training(i,5)*_Training(i,6)+_Parameters(86,0)*pow(_Training(i,6),2)+2*_Parameters(92,0)*_Training(i,1)*_Training(i,7)+_Parameters(98,0)*_Training(i,2)*_Training(i,7)+_Parameters(99,0)*_Training(i,3)*_Training(i,7)+_Parameters(100,0)*_Training(i,4)*_Training(i,7)+_Parameters(101,0)*_Training(i,5)*_Training(i,7)+_Parameters(102,0)*_Training(i,6)*_Training(i,7)+_Parameters(113,0)*pow(_Training(i,7),2);
                        
                        _Jacobian(6,1) = _Parameters(2,0)+2*_Parameters(10,0)*_Training(i,2)+_Parameters(9,0)*_Training(i,1)+_Parameters(12,0)*_Training(i,3)+_Parameters(15,0)*_Training(i,4)+_Parameters(19,0)*_Training(i,5)+_Parameters(24,0)*_Training(i,6)+_Parameters(30,0)*_Training(i,7)+3*_Parameters(39,0)*pow(_Training(i,2),2)+2*_Parameters(38,0)*_Training(i,1)*_Training(i,2)+_Parameters(37,0)*pow(_Training(i,1),2)+2*_Parameters(41,0)*_Training(i,2)*_Training(i,3)+_Parameters(42,0)*_Training(i,1)*_Training(i,3)+_Parameters(44,0)*pow(_Training(i,3),2)+2*_Parameters(47,0)*_Training(i,2)*_Training(i,4)+_Parameters(49,0)*_Training(i,1)*_Training(i,4)+_Parameters(51,0)*_Training(i,3)*_Training(i,4)+_Parameters(53,0)*pow(_Training(i,4),2)+2*_Parameters(57,0)*_Training(i,2)*_Training(i,5)+_Parameters(60,0)*_Training(i,1)*_Training(i,5)+_Parameters(63,0)*_Training(i,3)*_Training(i,5)+_Parameters(64,0)*_Training(i,4)*_Training(i,5)+_Parameters(67,0)*pow(_Training(i,5),2)+2*_Parameters(72,0)*_Training(i,2)*_Training(i,6)+_Parameters(76,0)*_Training(i,1)*_Training(i,6)+_Parameters(80,0)*_Training(i,3)*_Training(i,6)+_Parameters(81,0)*_Training(i,4)*_Training(i,6)+_Parameters(82,0)*_Training(i,5)*_Training(i,6)+_Parameters(87,0)*pow(_Training(i,6),2)+2*_Parameters(93,0)*_Training(i,2)*_Training(i,7)+_Parameters(98,0)*_Training(i,1)*_Training(i,7)+_Parameters(103,0)*_Training(i,3)*_Training(i,7)+_Parameters(104,0)*_Training(i,4)*_Training(i,7)+_Parameters(105,0)*_Training(i,5)*_Training(i,7)+_Parameters(106,0)*_Training(i,6)*_Training(i,7)+_Parameters(114,0)*pow(_Training(i,7),2);
                        
                        ;
                        _Jacobian(6,2) = _Parameters(3,0)+2*_Parameters(13,0)*_Training(i,3)+_Parameters(12,0)*_Training(i,2)+_Parameters(11,0)*_Training(i,1)+_Parameters(16,0)*_Training(i,4)+_Parameters(20,0)*_Training(i,5)+_Parameters(25,0)*_Training(i,6)+_Parameters(31,0)*_Training(i,7)+3*_Parameters(45,0)*pow(_Training(i,3),2)+2*_Parameters(44,0)*_Training(i,3)*_Training(i,2)+_Parameters(41,0)*pow(_Training(i,2),2)+2*_Parameters(43,0)*_Training(i,1)*_Training(i,3)+_Parameters(42,0)*_Training(i,1)*_Training(i,2)+_Parameters(40,0)*pow(_Training(i,1),2)+2*_Parameters(48,0)*_Training(i,3)*_Training(i,4)+_Parameters(50,0)*_Training(i,1)*_Training(i,4)+_Parameters(51,0)*_Training(i,2)*_Training(i,4)+_Parameters(54,0)*pow(_Training(i,4),2)+2*_Parameters(58,0)*_Training(i,3)*_Training(i,5)+_Parameters(61,0)*_Training(i,1)*_Training(i,5)+_Parameters(63,0)*_Training(i,2)*_Training(i,5)+_Parameters(65,0)*_Training(i,4)*_Training(i,5)+_Parameters(68,0)*pow(_Training(i,5),2)+2*_Parameters(73,0)*_Training(i,3)*_Training(i,6)+_Parameters(77,0)*_Training(i,1)*_Training(i,6)+_Parameters(80,0)*_Training(i,2)*_Training(i,6)+_Parameters(83,0)*_Training(i,4)*_Training(i,6)+_Parameters(84,0)*_Training(i,5)*_Training(i,6)+_Parameters(88,0)*pow(_Training(i,6),2)+2*_Parameters(94,0)*_Training(i,3)*_Training(i,7)+_Parameters(99,0)*_Training(i,1)*_Training(i,7)+_Parameters(103,0)*_Training(i,2)*_Training(i,7)+_Parameters(107,0)*_Training(i,4)*_Training(i,7)+_Parameters(108,0)*_Training(i,5)*_Training(i,7)+_Parameters(109,0)*_Training(i,6)*_Training(i,7)+_Parameters(115,0)*pow(_Training(i,7),2);
                        
                        
                        _Jacobian(6,3) = _Parameters(4,0)+_Parameters(14,0)*_Training(i,1)+_Parameters(15,0)*_Training(i,2)+_Parameters(16,0)*_Training(i,3)+2*_Parameters(17,0)*_Training(i,4)+_Parameters(21,0)*_Training(i,5)+_Parameters(26,0)*_Training(i,6)+_Parameters(32,0)*_Training(i,7)+_Parameters(46,0)*pow(_Training(i,1),2)+_Parameters(47,0)*pow(_Training(i,2),2)+_Parameters(48,0)*pow(_Training(i,3),2)+_Parameters(49,0)*_Training(i,1)*_Training(i,2)+_Parameters(50,0)*_Training(i,1)*_Training(i,3)+_Parameters(51,0)*_Training(i,2)*_Training(i,3)+2*_Parameters(52,0)*_Training(i,1)*_Training(i,4)+2*_Parameters(53,0)*_Training(i,2)*_Training(i,4)+2*_Parameters(53,0)*_Training(i,3)*_Training(i,4)+3*_Parameters(55,0)*pow(_Training(i,4),2)+2*_Parameters(59,0)*_Training(i,4)*_Training(i,5)+_Parameters(62,0)*_Training(i,1)*_Training(i,5)+_Parameters(64,0)*_Training(i,2)*_Training(i,5)+_Parameters(65,0)*_Training(i,3)*_Training(i,5)+_Parameters(69,0)*pow(_Training(i,5),2)+2*_Parameters(74,0)*_Training(i,4)*_Training(i,6)+_Parameters(78,0)*_Training(i,1)*_Training(i,6)+_Parameters(81,0)*_Training(i,2)*_Training(i,6)+_Parameters(83,0)*_Training(i,3)*_Training(i,6)+_Parameters(85,0)*_Training(i,5)*_Training(i,6)+_Parameters(89,0)*pow(_Training(i,6),2)+2*_Parameters(95,0)*_Training(i,1)*_Training(i,7)+_Parameters(100,0)*_Training(i,1)*_Training(i,7)+_Parameters(104,0)*_Training(i,2)*_Training(i,7)+_Parameters(107,0)*_Training(i,3)*_Training(i,7)+_Parameters(110,0)*_Training(i,5)*_Training(i,7)+_Parameters(111,0)*_Training(i,6)*_Training(i,7)+_Parameters(116,0)*pow(_Training(i,7),2);
                        
                        _Jacobian(6,4) = _Parameters(5,0)+_Parameters(18,0)*_Training(i,1)+_Parameters(19,0)*_Training(i,2)+_Parameters(20,0)*_Training(i,3)+_Parameters(21,0)*_Training(i,4)+2*_Parameters(22,0)*_Training(i,5)+_Parameters(27,0)*_Training(i,6)+_Parameters(33,0)*_Training(i,7)+_Parameters(56,0)*pow(_Training(i,1),2)+_Parameters(57,0)*pow(_Training(i,2),2)+_Parameters(58,0)*pow(_Training(i,3),2)+_Parameters(59,0)*pow(_Training(i,4),2)+_Parameters(60,0)*_Training(i,1)*_Training(i,2)+_Parameters(61,0)*_Training(i,1)*_Training(i,3)+_Parameters(62,0)*_Training(i,1)*_Training(i,4)+_Parameters(63,0)*_Training(i,2)*_Training(i,3)+_Parameters(64,0)*_Training(i,2)*_Training(i,4)+_Parameters(65,0)*_Training(i,3)*_Training(i,4)+2*_Parameters(66,0)*_Training(i,1)*_Training(i,5)+2*_Parameters(67,0)*_Training(i,2)*_Training(i,5)+2*_Parameters(68,0)*_Training(i,3)*_Training(i,5)+2*_Parameters(69,0)*_Training(i,4)*_Training(i,5)+3*_Parameters(70,0)*pow(_Training(i,5),2)+2*_Parameters(75,0)*_Training(i,5)*_Training(i,6)+_Parameters(79,0)*_Training(i,1)*_Training(i,6)+_Parameters(82,0)*_Training(i,2)*_Training(i,6)+_Parameters(84,0)*_Training(i,3)*_Training(i,6)+_Parameters(85,0)*_Training(i,4)*_Training(i,6)+_Parameters(90,0)*pow(_Training(i,6),2)+2*_Parameters(96,0)*_Training(i,5)*_Training(i,7)+_Parameters(101,0)*_Training(i,1)*_Training(i,7)+_Parameters(105,0)*_Training(i,2)*_Training(i,7)+_Parameters(108,0)*_Training(i,3)*_Training(i,7)+_Parameters(110,0)*_Training(i,4)*_Training(i,7)+_Parameters(112,0)*_Training(i,6)*_Training(i,7)+_Parameters(117,0)*pow(_Training(i,7),2);
                        
                        _Jacobian(6,5) = _Parameters(6,0)+_Parameters(23,0)*_Training(i,1)+_Parameters(24,0)*_Training(i,2)+_Parameters(25,0)*_Training(i,3)+_Parameters(26,0)*_Training(i,4)+_Parameters(27,0)*_Training(i,5)+2*_Parameters(28,0)*_Training(i,6)+_Parameters(34,0)*_Training(i,7)+_Parameters(71,0)*pow(_Training(i,1),2)+_Parameters(72,0)*pow(_Training(i,2),2)+_Parameters(73,0)*pow(_Training(i,3),2)+_Parameters(74,0)*pow(_Training(i,4),2)+_Parameters(75,0)*pow(_Training(i,5),2)+_Parameters(76,0)*_Training(i,1)*_Training(i,2)+_Parameters(77,0)*_Training(i,1)*_Training(i,3)+_Parameters(78,0)*_Training(i,1)*_Training(i,4)+_Parameters(79,0)*_Training(i,1)*_Training(i,5)+_Parameters(80,0)*_Training(i,2)*_Training(i,3)+_Parameters(81,0)*_Training(i,2)*_Training(i,4)+_Parameters(82,0)*_Training(i,2)*_Training(i,5)+_Parameters(83,0)*_Training(i,3)*_Training(i,4)+_Parameters(84,0)*_Training(i,3)*_Training(i,5)+_Parameters(85,0)*_Training(i,4)*_Training(i,5)+2*_Parameters(86,0)*_Training(i,1)*_Training(i,6)+2*_Parameters(87,0)*_Training(i,2)*_Training(i,6)+2*_Parameters(88,0)*_Training(i,3)*_Training(i,6)+2*_Parameters(89,0)*_Training(i,4)*_Training(i,6)+2*_Parameters(90,0)*_Training(i,5)*_Training(i,6)+3*_Parameters(91,0)*pow(_Training(i,6),2)+2*_Parameters(97,0)*_Training(i,6)*_Training(i,7)+_Parameters(102,0)*_Training(i,1)*_Training(i,7)+_Parameters(106,0)*_Training(i,2)*_Training(i,7)+_Parameters(109,0)*_Training(i,3)*_Training(i,7)+_Parameters(111,0)*_Training(i,4)*_Training(i,7)+_Parameters(112,0)*_Training(i,5)*_Training(i,7)+_Parameters(118,0)*pow(_Training(i,7),2);
                        
                        _Jacobian(6,6) = _Parameters(7,0)+_Parameters(29,0)*_Training(i,1)+_Parameters(30,0)*_Training(i,2)+_Parameters(31,0)*_Training(i,3)+_Parameters(32,0)*_Training(i,4)+_Parameters(33,0)*_Training(i,5)+_Parameters(34,0)*_Training(i,6)+2*_Parameters(35,0)*_Training(i,7)+_Parameters(92,0)*pow(_Training(i,1),2)+_Parameters(93,0)*pow(_Training(i,2),2)+_Parameters(94,0)*pow(_Training(i,3),2)+_Parameters(95,0)*pow(_Training(i,4),2)+_Parameters(96,0)*pow(_Training(i,5),2)+_Parameters(97,0)*pow(_Training(i,6),2)+_Parameters(98,0)*_Training(i,1)*_Training(i,2)+_Parameters(99,0)*_Training(i,1)*_Training(i,3)+_Parameters(100,0)*_Training(i,1)*_Training(i,4)+_Parameters(101,0)*_Training(i,1)*_Training(i,5)+_Parameters(102,0)*_Training(i,1)*_Training(i,6)+_Parameters(103,0)*_Training(i,2)*_Training(i,3)+_Parameters(104,0)*_Training(i,2)*_Training(i,4)+_Parameters(105,0)*_Training(i,2)*_Training(i,5)+_Parameters(106,0)*_Training(i,2)*_Training(i,6)+_Parameters(107,0)*_Training(i,3)*_Training(i,4)+_Parameters(108,0)*_Training(i,3)*_Training(i,5)+_Parameters(109,0)*_Training(i,3)*_Training(i,6)+_Parameters(110,0)*_Training(i,4)*_Training(i,5)+_Parameters(111,0)*_Training(i,4)*_Training(i,6)+_Parameters(112,0)*_Training(i,5)*_Training(i,6)+2*_Parameters(113,0)*_Training(i,1)*_Training(i,7)+2*_Parameters(114,0)*_Training(i,2)*_Training(i,7)+2*_Parameters(115,0)*_Training(i,3)*_Training(i,7)+2*_Parameters(116,0)*_Training(i,4)*_Training(i,7)+2*_Parameters(117,0)*_Training(i,5)*_Training(i,7)+2*_Parameters(118,0)*_Training(i,6)*_Training(i,7)+3*_Parameters(119,0)*pow(_Training(i,7),2);
                        
                        break;
                    }
                    break;
                    
            }
            
        }
        
        
        double f = 0.0;
        return f;
    }
    
    
    
    // Estimating the Lyapunov exponent second approach.
	template <typename Embedding, typename Lyapunov, typename EA>
	double lyapunov_estimation(Embedding& d , Lyapunov& l, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        // input data can be found here (defined in config file or command line):
        
        matrix_type _InputMatrix;
        matrix_type _InitialMatrix;
        vector_type_distance _EuclideanDistance;
        _InputMatrix.resize(_IntegerInput.size() - d + 1, d);
        _InitialMatrix.resize(_IntegerInput.size() - d + 1, d);
        _EuclideanDistance.resize(_IntegerInput.size() - d + 1);
        _EuclideanDistance(0) = std::numeric_limits<double>::max();
        
        for (int i = 0 ; i < d ; i++)
        {
            matrix_type ShiftedTemp(_IntegerInput.begin() + i, _IntegerInput.end() - d + i + 1 , 1);
            matrix_type Temp(_IntegerInput.begin() , _IntegerInput.end() - d + 1 , 1 , _IntegerInput(i));
            column(_InputMatrix, i) = column(ShiftedTemp , 0);
            column(_InitialMatrix, i) = column(Temp , 0);
        }
        
        
        unsigned Index;
        l = 0; //Largest LE
        double MinDistance = std::numeric_limits<double>::max();
        
        
        double EvolZero;
        double EvolPrime;
        
        for (unsigned i = 0 ; i < _EuclideanDistance.size() - 1; i++)
        {
            
            row_type InitialRow(_InitialMatrix, i);
            
            for (int j = 0 ; j < _EuclideanDistance.size() ; j++)
            {
                row_type InputRow(_InputMatrix, j);
                _EuclideanDistance(j) = norm_2(InputRow - InitialRow);
            }
            
            _EuclideanDistance(i) = MinDistance;
            
            Index = 0;
            for (unsigned j = 0 ; j < _EuclideanDistance.size() ; j++)
            {
                if (_EuclideanDistance(j) < MinDistance)
                {
                    MinDistance = _EuclideanDistance(j);
                    Index = j;
                }
            }
            

            
            row_type EvolRowZero(_InputMatrix, i);
            row_type EvolRowOne (_InputMatrix, Index);
            
            EvolZero = norm_2(EvolRowOne - EvolRowZero);
            
            
            row_type EvolRowPrime    (_InputMatrix, i + 1);
            row_type EvolRowPrimeOne (_InputMatrix, Index + 1);
            
            EvolPrime = norm_2(EvolRowPrimeOne - EvolRowPrime);
            
            l += log2(EvolPrime/EvolZero)/(_EuclideanDistance.size() - 1);

        }
        
        return l;
    }
    
    
    
    // Estimating the prediction horizon.
	template <typename PredictionHorizon, typename Lyapunov, typename EA>
	unsigned prediction_horizon_estimation(PredictionHorizon& h , Lyapunov& l, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        // input data can be found here (defined in config file or command line):
        h = unsigned(1.0 / l);
        return h;
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
