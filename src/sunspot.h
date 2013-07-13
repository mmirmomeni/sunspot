#ifndef _SUNSPOT_H_
#define _SUNSPOT_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cmath>
#include <ea/fitness_function.h>
#include <ea/meta_data.h>
#include <mkv/parse.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace ealib;

// Common meta-data needed for sunspot prediction.
LIBEA_MD_DECL(SUNSPOT_INPUT, "sunspot.input", std::string);
LIBEA_MD_DECL(SUNSPOT_TEST_INPUT, "sunspot.test_input", std::string);
LIBEA_MD_DECL(SUNSPOT_PREDICTION_HORIZON, "sunspot.prediction_horizon", std::size_t);

/*! Fitness function for sunspot number prediction.
 */
struct sunspot_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef boost::numeric::ublas::matrix<int> matrix_type; //!< Type for matrix that will store raw sunspot numbers.
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::matrix_column<matrix_type> column_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::vector<int> vector_type; //!< Type for a vector of sunspot numbers.
    typedef boost::numeric::ublas::vector<double> dvector_type; //!< Type for a vector of sunspot numbers.
    
    mkv::markov_network::desc_type _desc; //!< Description of the markov network we'll be using.

    // training data:
    matrix_type _train_input; //!< All historical sunspot number data used during fitness evaluation (inputs to MKV network).
    vector_type _train_observed; //!< Observed (real) historical sunspot numbers.

    // testing data:
    matrix_type _test_input; //!< All historical sunspot number data used during fitness evaluation (inputs to MKV network).
    vector_type _test_observed; //!< Observed (real) historical sunspot numbers.
    
    //! Initialize this fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        mkv::parse_desc(get<MKV_DESC>(ea), _desc);
        
        load_file(get<SUNSPOT_INPUT>(ea), _train_input, _train_observed);
        load_file(get<SUNSPOT_TEST_INPUT>(ea), _test_input, _test_observed);
    }
    
    void load_file(const std::string& filename, matrix_type& input, vector_type& observed) {
        // read in the training data:
        std::ifstream infile(filename.c_str());
        std::string line;
        int nrow=0;
        if(infile.is_open()) {
            for(int i=0; i<4; ++i) {
                getline(infile,line);
            }
            infile >> nrow;
        } else {
            throw file_io_exception("sunspot.h::initialize: could not open " + filename);
        }
        
<<<<<<< HEAD
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
            _Training(i,48)  = pow(_IntegerObservedED(i + MAX_ED - 3),2)*_IntegerObservedED(i + MAX_ED - 4);  //###########//
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
=======
        const int ncol = 9; // fixed number of bits per ssn

        // _input: a matrix where each row vector i is assumed to be the complete
        // **binary input vector** to the MKV network at time i.
        input.resize(nrow,ncol);
>>>>>>> 2d55879d34291d42095ec9562117a3ac98d3d0cd
        
        // _observed: a vector where element i corresponds to the observed (real) sunspot
        // number corresponding to row i in _inputs.
        observed.resize(nrow);

        for(int i=0; i<nrow; ++i) {
            int t0=0;
            infile >> t0 >> observed(i);

            std::bitset<ncol> t0b = ~std::bitset<ncol>(t0);
            for(int j=0; j<ncol; ++j) {
                input(i,j) = t0b[ncol - j - 1];
            }

            // logic above is a bit strange.  consider instead:
            // std::bitset<ncol> ssn(t0);
            // for(int j=0; j<ncol; ++j) {
            //     _input(i,j) = ssn[j];
            // }
        }
    }
    
    /*! Test an individual for multiple predictions.
     
     Here the output of the Markov network is interpreted for each of n time steps.
     */
    template <typename Individual, typename RNG, typename EA>
	dvector_type eval(Individual& ind, matrix_type& output, matrix_type& input, vector_type& observed, RNG& rng, EA& ea, bool recurse=false) {
        namespace bnu=boost::numeric::ublas;
        mkv::markov_network net = mkv::make_markov_network(_desc, ind.repr().begin(), ind.repr().end(), rng.seed(), ea);
        
        // outputs from the MKV network, initialized to the same size as the number
        // of observations by the fixed prediction horizon of 8 time steps.
        std::size_t ph=get<SUNSPOT_PREDICTION_HORIZON>(ea);
        output.resize(input.size1(), ph);
        
        bool first=true;
        
        // run each row of _inputs through the MKV network for a single update,
        // place the results in the output matrix:
        for(std::size_t i=0; i<input.size1(); ++i) {
            if(recurse) {
                std::vector<int> v;
                if(first) {
                    row_type r(input,i);
                    std::copy(r.begin(), r.end(), std::back_inserter(v));
                    first = false;
                } else {
                    algorithm::range_pair2int(net.begin_output(),
                                              net.begin_output()+2*input.size2(),
                                              std::back_inserter(v));
                }
                mkv::update(net, 1, v.begin());
            } else {
                row_type r(input,i);
                mkv::update(net, 1, r.begin());
            }
            
            for(std::size_t j=0; j<ph; ++j) {
                output(i,j) = algorithm::range_pair2int(net.begin_output()+j*2*input.size2(),
                                                        net.begin_output()+(j+1)*2*input.size2());
            }
        }
        
        // fitness == 1.0/(1.0 + sum_{i=1}^{prediction horizon} RMSE_i)
        dvector_type rmse(get<SUNSPOT_PREDICTION_HORIZON>(ea));
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            column_type c(output,i);
            assert(c.size() == observed.size());
            
            // given:
            // input | observed | output (predictions)
            //                  | ph0 ph1 ph2 ph3...
            // -------------------------------------
            // i0    | i1       | o1,1 o1,2 o1,3 o1,4
            // i1    | i2       | o2,1 o2,2 o2,3 o2,4
            // i2    | i3       | o3,1 o3,2 o3,3 o3,4
            //
            // we need to match the columns of the output matrix with the
            // appropriate range in the "observed" column, e.g.:
            //
            // err += observed(i1,i2,i3...) - output(o1,1, o2,1, o3,1...)
            // err += observed(i2,i3,i4...) - output(o1,2, o2,2, o3,2...)
            //
            // note that as i increases, we're dropping elements off the front
            // of the observed, and the back of the output matrix.
            vector_type err =
            bnu::vector_range<vector_type>(observed, bnu::range(i,observed.size()))
            - bnu::vector_range<column_type>(c, bnu::range(0,c.size()-i));
            
            rmse(i) = sqrt(static_cast<double>(bnu::inner_prod(err,err)) / static_cast<double>(err.size()));
        }

        return rmse;
    };
    
    template <typename Individual, typename RNG, typename EA>
    dvector_type train(Individual& ind, matrix_type& output, RNG& rng, EA& ea) {
        return eval(ind, output, _train_input, _train_observed, rng, ea);
    }

    template <typename Individual, typename RNG, typename EA>
    dvector_type test(Individual& ind, matrix_type& output, RNG& rng, EA& ea, bool recurse=false) {
        return eval(ind, output, _test_input, _test_observed, rng, ea, recurse);
    }

    //! Calculate fitness of an individual.
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        typedef bnu::vector<int> int_vector_type;
        matrix_type output;
        vector_type rmse = train(ind, output, rng, ea);
        
        return 100.0 / (1.0 + std::accumulate(rmse.begin(), rmse.end(), 0.0));
    }
};

template <typename EA>
struct sunspot_detail : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_detail";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_detail.dat");
        df.add_field("observed");
        
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().train(ind, output, ea.rng(), ea);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._train_observed(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }

            df.endl();
        }
    }
};

template <typename EA>
struct sunspot_test : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_test";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_test.dat");
        df.add_field("observed");
        
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().test(ind, output, ea.rng(), ea);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._test_observed(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }
            
            df.endl();
        }
    }
};

template <typename EA>
struct sunspot_test_rmse : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_test_rmse";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_test_rmse.dat");
        df.add_field("total_rmse");
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        sunspot_fitness::dvector_type rmse = ea.fitness_function().test(ind, output, ea.rng(), ea);
        
        df.write(std::accumulate(rmse.begin(),rmse.end(),0.0)).write_all(rmse.begin(), rmse.end()).endl();
    }
};

template <typename EA>
struct sunspot_recursive_test : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_recursive_test";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_recursive_test.dat");
        df.add_field("observed");
        
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().test(ind, output, ea.rng(), ea, true);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._test_observed(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }
            
            df.endl();
        }
    }
};

#endif
