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
        
        const int ncol = 9; // fixed number of bits per ssn

        // _input: a matrix where each row vector i is assumed to be the complete
        // **binary input vector** to the MKV network at time i.
        input.resize(nrow,ncol);
        
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
