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

using namespace ealib;

// Common meta-data needed for sunspot prediction.
LIBEA_MD_DECL(SUNSPOT_INPUT, "sunspot.input", std::string);
LIBEA_MD_DECL(SUNSPOT_PREDICTION_HORIZON, "sunspot.prediction_horizon", std::size_t);

/*! Fitness function for sunspot number prediction.
 */
struct sunspot_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef boost::numeric::ublas::matrix<int> matrix_type; //!< Type for matrix that will store raw sunspot numbers.
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
    typedef boost::numeric::ublas::matrix_column<matrix_type> column_type; //!< Row type for the matrix.
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
    
    /*! Test an individual for multiple predictions.
     
     Here the output of the Markov network is interpreted for each of n time steps.
     */
    template <typename Individual, typename RNG, typename EA>
	void test(Individual& ind, matrix_type& output, RNG& rng, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        mkv::markov_network net = mkv::make_markov_network(_desc, ind.repr().begin(), ind.repr().end(), rng.seed(), ea);
        
        // outputs from the MKV network, initialized to the same size as the number
        // of observations by the fixed prediction horizon of 8 time steps.
        std::size_t ph=get<SUNSPOT_PREDICTION_HORIZON>(ea);
        output.resize(_observed.size1(),ph);
        
        // run each row of _inputs through the MKV network for a single update,
        // place the results in the output matrix:
        for(std::size_t i=0; i<_input.size1(); ++i) {
            row_type r(_input,i);
            mkv::update(net, 1, r.begin());
            
            for(std::size_t j=0; j<ph; ++j) {
                output(i,j) = static_cast<double>(algorithm::range_pair2int(net.begin_output()+j*_input.size2(),
                                                                            net.begin_output()+(j+1)*_input.size2()));
            }
        }
    };
    
    //! Calculate fitness of an individual.
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        namespace bnu=boost::numeric::ublas;
        matrix_type output;
        test(ind, output, rng, ea);
        
        // fitness == 1.0/(1.0 + sum_{i=1}^{prediction horizon} RMSE_i)
        double rmse=0.0;
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            column_type c(output,i);
            
            // careful; have to lag these correctly...
            bnu::vector<int> err = _IntegerObserved - c;
            
            rmse += sqrt(bnu::inner_prod(err,err) / static_cast<double>(err.size()));
        }
        
        return 100.0 / (1.0 + rmse);
    }
};

/*! Save the detailed graph of the dominant individual in graphviz format.
 */
template <typename EA>
struct test_sunspot : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "test_sunspot";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("test_sunspot.dat");
        df.add_field("observed");
        
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("predicted_tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().test(ind, output, ea.rng(), ea);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._IntegerObserved(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }

            df.endl();
        }
    }
};

#endif
