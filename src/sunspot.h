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
#include <cmath>
#include <ea/fitness_function.h>
#include <ea/meta_data.h>
#include <mkv/parse.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace ea;

// Common meta-data needed for sunspot prediction.
LIBEA_MD_DECL(SUNSPOT_INPUT, "sunspot.input", std::string);

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
