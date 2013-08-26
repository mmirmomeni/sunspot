#ifndef _SUNSPOT_H_
#define _SUNSPOT_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <cmath>
#include <ea/fitness_function.h>
#include <ea/meta_data.h>

#include <fstream>
#include <string>
#include <iostream>

using namespace ealib;

// training data filename:
LIBEA_MD_DECL(SUNSPOT_TRAIN, "sunspot.train_filename", std::string);
// testing data filename:
LIBEA_MD_DECL(SUNSPOT_TEST, "sunspot.test_filename", std::string);
// # of bits left and right of radix point (for fixed-point representations):
LIBEA_MD_DECL(SUNSPOT_INTEGER_BITS, "sunspot.integer_bits", std::size_t);
LIBEA_MD_DECL(SUNSPOT_FRACTIONAL_BITS, "sunspot.fractional_bits", std::size_t);
// how many time steps into the future we generate predictions for:
LIBEA_MD_DECL(SUNSPOT_PREDICTION_HORIZON, "sunspot.prediction_horizon", std::size_t);
// how many lags are to be included in the inputs?
LIBEA_MD_DECL(SUNSPOT_INPUT_LAGS, "sunspot.input_lags", std::size_t);
// limit the size of the dataset?
LIBEA_MD_DECL(SUNSPOT_LIMIT, "sunspot.limit", std::size_t);


/*! Fitness function for sunspot number prediction.
 */
struct sunspot_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef boost::numeric::ublas::matrix<int> matrix_type; //!< Type for matrix that will store raw sunspot numbers.
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::matrix_column<matrix_type> column_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::vector<int> vector_type; //!< Type for a vector of sunspot numbers.
    typedef boost::numeric::ublas::vector<double> dvector_type; //!< Type for a vector of sunspot numbers.
    
    // training data:
    matrix_type _train_input; //!< Inputs to the Markov network during fitness evaluation.
    vector_type _train_t; //!< Values at time t (corresponds to input).
    vector_type _train_tplus1; //!< Values at time t+1 (i.e., the value to be predicted).
    
    // testing data:
    matrix_type _test_input; //!< Inputs to the Markov network during fitness evaluation.
    vector_type _test_t; //!< Values at time t (corresponds to input).
    vector_type _test_tplus1; //!< Values at time t+1 (i.e., the value to be predicted).
    
    enum { IDX_T=0, IDX_X=1 }; // indices into the data
    
    //! Initialize this fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        load_file(get<SUNSPOT_TRAIN>(ea), _train_input, _train_t, _train_tplus1, ea);
        load_file(get<SUNSPOT_TEST>(ea), _test_input, _test_t, _test_tplus1, ea);
    }
    
    
    /*! Load training or testing data from a file.
     
     The file format is as follows:
     
     # commented lines begin with a '#'
     # the first non-commented line is assumed to contain the header
     time x
     # all following non-comment, non-empty lines are to contain at least two
     # decimal numbers, the first being the time, and the second the observed value
     # of the time seriers
     0 75
     1 23
     # in the case of a decimal number, e.g.:
     2 87.4
     # ... x will be converted to a fixed-point number.
     
     This loader is a bit more complicated than it otherwise could be, but we're 
     also transparently handling integer vs. fixed-point numbers and trying to be
     very conservative with regard to the datafile format.
     
     \note this has been tested to produce the same ouput as the original loader,
     with two changes: the original form of this data was (a) inverted, and 
     (b) reversed the MSB and LSB.  we "correct" those here, but may have to undo 
     them depending on performance.

     \todo rewrite using boost::spirit?
     */
    template <typename EA>
    void load_file(const std::string& filename, matrix_type& input, vector_type& t, vector_type& tplus1, EA& ea) {
        using namespace boost;
        namespace bnu=boost::numeric::ublas;

        // read in the raw training data, split by field, store it in varsize string matrix:
        std::ifstream infile(filename.c_str());
        
        if(!infile.is_open()) {
            throw file_io_exception("sunspot.h::initialize: could not open " + filename);
        }
        
        static const regex comment("^#.*");
        static const regex decimal("\\d+\\.\\d+");
        bool header=false;
        std::string line;
        std::vector<std::vector<std::string> > smat;

        while(getline(infile, line)) {
            // remove leading & trailing whitespace:
            trim(line);
            // only consider lines with length > 0:
            if(line.empty()) {
                continue;
            }            
            // skip all comments:
            if(regex_match(line, comment)) {
                continue;
            }
            // the first non-comment line is assumed to be the header:
            if(!header) {
                header = true;
                continue;
            }
            
            // split all remaining lines into fields, add them to the string matrix:
            std::vector<std::string> fields;
            split(fields, line, is_space());
            smat.push_back(fields);
        }
        
        // are we limiting the size of the dataset?
        if(get<SUNSPOT_LIMIT>(ea,0)>0) {
            smat.resize(std::min(smat.size(), get<SUNSPOT_LIMIT>(ea)));
        }
        
        // now convert the elements in the string matrix to our input matrix
        // and observed vector:
        std::size_t nbits = get<SUNSPOT_INTEGER_BITS>(ea) + get<SUNSPOT_FRACTIONAL_BITS>(ea);
        std::size_t lags = get<SUNSPOT_INPUT_LAGS>(ea);
        long factor = 1 << get<SUNSPOT_FRACTIONAL_BITS>(ea);
        long maxval = 1 << nbits;

        vector_type raw(smat.size());
        
        for(std::size_t i=0; i<smat.size(); ++i) {
            // we're going to start just by looking at the value (2nd col, IDX_X):
            long x=0;
            
            // if we have a decimal in the string, we assume that we're reading a double:
            if(regex_match(smat[i][IDX_X], decimal)) {
                x = static_cast<long>(boost::lexical_cast<double>(smat[i][IDX_X]) * factor);
            } else {
                // it's already an integer; we're done
                x = boost::lexical_cast<long>(smat[i][IDX_X]) << get<SUNSPOT_FRACTIONAL_BITS>(ea);
            }
            
            // check for a bad value:
            if(x > maxval) {
                throw bad_argument_exception("sunspot.h::load: bad value in datafile, " + boost::lexical_cast<std::string>(x));
            }
            
            raw(i) = x;
        }
        
        // build the t and t+1 vectors from the raw data:
        t = bnu::vector_range<vector_type>(raw, bnu::range(lags-1,raw.size()-1));
        tplus1 = bnu::vector_range<vector_type>(raw, bnu::range(lags,raw.size()));
        assert(t.size() == tplus1.size());
        
        // now we need to construct the input matrix, where each row of the matrix
        // is an input to the markov network.  we'll start by building the integer
        // version of the matrix, like so:
        //
        // lags==1:
        // i=1: imat.col(0) = raw(0,size-1)
        //
        // lags==2:
        // i=1: imat.col(0) = raw(1,size-1)
        // i=2: imat.col(1) = raw(0,size-2)
        //
        // lags==3:
        // i=1: imat.col(0) = raw(2,size-1)
        // i=2: imat.col(1) = raw(1,size-2)
        // i=3: imat.col(2) = raw(0,size-3)
        // ...
        matrix_type imat(raw.size()-lags,lags); // lags >= 1
        for(std::size_t i=1; i<=lags; ++i) {
            column_type c(imat,i-1);
            c = bnu::vector_range<vector_type>(raw, bnu::range(lags-i, raw.size()-i));
        }

        // now convert this thing to binary, and we're all done:
        input.resize(imat.size1(), nbits*imat.size2());
        for(std::size_t i=0; i<imat.size1(); ++i) {
            for(std::size_t j=0; j<imat.size2(); ++j) {
                std::bitset<sizeof(long)*8> bits(imat(i,j));
                for(std::size_t k=0; k<nbits; ++k) {
                    input(i,j*nbits+k) = bits[k];
                }
            }
        }
        
        assert(input.size1() == t.size());
    }
    
    /*! Test an individual for multiple predictions.
     
     Here the output of the Markov network is interpreted for each of n time steps.
     */
    template <typename Individual, typename RNG, typename EA>
	dvector_type eval(Individual& ind, matrix_type& output, matrix_type& input, vector_type& observed, RNG& rng, EA& ea, bool recurse=false) {
        namespace bnu=boost::numeric::ublas;
        mkv::markov_network& net = ealib::phenotype(ind,rng,ea);

        // number of network updates per sample:
        const std::size_t nupdates = get<mkv::MKV_UPDATE_N>(ea);
        
        // outputs from the MKV network, initialized to the same size as the number
        // of observations by the fixed prediction horizon of 8 time steps.
        const std::size_t ph=get<SUNSPOT_PREDICTION_HORIZON>(ea);
        output.resize(input.size1(), ph);
        
        // size of input and output variables, in bits:
        const std::size_t nbits = get<SUNSPOT_INTEGER_BITS>(ea) + get<SUNSPOT_FRACTIONAL_BITS>(ea);

        // run each row of _inputs through the MKV network:
        for(std::size_t i=0; i<input.size1(); ++i) {
            row_type r(input,i);
            mkv::update(net, 1, nupdates, r);
            
            for(std::size_t j=0; j<ph; ++j) {
                output(i,j) = algorithm::range_pair2int(net.begin_output()+j*2*nbits,
                                                        net.begin_output()+(j+1)*2*nbits);
            }
        }
        
        // fitness == 1.0/(1.0 + sum_{i=1}^{prediction horizon} RMSE_i)
        dvector_type rmse(ph);
        for(std::size_t i=0; i<ph; ++i) {
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
            dvector_type err =
            bnu::vector_range<vector_type>(observed, bnu::range(i,observed.size()))
            - bnu::vector_range<column_type>(c, bnu::range(0,c.size()-i));
            
            // err is a possibly-fixed-point integer vector; convert it to floating
            // point, and calculate rmse:
            long factor = 1 << get<SUNSPOT_FRACTIONAL_BITS>(ea);
            dvector_type derr = err / static_cast<double>(factor);
            
            double r = sqrt(bnu::inner_prod(derr,derr) / static_cast<double>(derr.size()));
            rmse(i) = r;
        }
        
        return rmse;
    };
    
    template <typename Individual, typename RNG, typename EA>
    dvector_type train(Individual& ind, matrix_type& output, RNG& rng, EA& ea) {
        return eval(ind, output, _train_input, _train_tplus1, rng, ea);
    }
    
    template <typename Individual, typename RNG, typename EA>
    dvector_type test(Individual& ind, matrix_type& output, RNG& rng, EA& ea, bool recurse=false) {
        return eval(ind, output, _test_input, _test_tplus1, rng, ea, recurse);
    }
    
    //! Calculate fitness of an individual.
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        matrix_type output;
        dvector_type rmse = train(ind, output, rng, ea);
        return 100.0 / (1.0 + std::accumulate(rmse.begin(), rmse.end(), 0.0));
    }
};


#include "analysis.h"

#endif
