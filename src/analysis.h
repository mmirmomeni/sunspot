#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

template <typename EA>
struct sunspot_data : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_data";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        
        typename EA::fitness_function_type::matrix_type& M=ea.fitness_function()._train_input;
        typename EA::fitness_function_type::vector_type& T=ea.fitness_function()._train_t;
        typename EA::fitness_function_type::vector_type& Tp1=ea.fitness_function()._train_tplus1;
        
        datafile df("sunspot_data.dat");
        df.comment("WARNING: bits are reversed (i0=LSB) and multiple numbers may be encoded per line");
        df.add_field("t");
        for(std::size_t i=0; i<M.size2(); ++i) {
            df.add_field("i" + boost::lexical_cast<std::string>(i));
        }
        df.add_field("tplus1");
        
        for(std::size_t i=0; i<M.size1(); ++i) {
            typename EA::fitness_function_type::row_type r(M,i);
            df.write(T(i))
            .write_all(r.begin(), r.end())
            .write(Tp1(i))
            .endl();
        }
    }
};


template <typename EA>
struct sunspot_train_detail : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_train_detail";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_train_detail.dat");
        df.add_field("observed");
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().train(ind, output, ea.rng(), ea);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._train_tplus1(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }
            
            df.endl();
        }
    }
};

template <typename EA>
struct sunspot_test_detail : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_test_detail";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_test_detail.dat");
        df.add_field("observed");
        
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        ea.fitness_function().test(ind, output, ea.rng(), ea);
        
        for(std::size_t i=0; i<output.size1(); ++i) {
            df.write(ea.fitness_function()._test_tplus1(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }
            
            df.endl();
        }
    }
};

template <typename EA>
struct sunspot_train_rmse : public ealib::analysis::unary_function<EA> {
    static const char* name() { return "sunspot_train_rmse";}
    
    virtual void operator()(EA& ea) {
        using namespace ealib;
        using namespace ealib::analysis;
        typename EA::individual_type& ind = analysis::find_dominant(ea);
        
        datafile df("sunspot_train_rmse.dat");
        df.add_field("total_rmse");
        for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
            df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
        }
        
        sunspot_fitness::matrix_type output;
        sunspot_fitness::dvector_type rmse = ea.fitness_function().train(ind, output, ea.rng(), ea);
        
        df.write(std::accumulate(rmse.begin(),rmse.end(),0.0)).write_all(rmse.begin(), rmse.end()).endl();
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
            df.write(ea.fitness_function()._test_tplus1(i));
            
            for(std::size_t j=0; j<output.size2(); ++j) {
                df.write(output(i,j));
            }
            
            df.endl();
        }
    }
};

#endif
