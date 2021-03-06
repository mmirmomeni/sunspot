/* sunspot.cpp
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
#include <ea/evolutionary_algorithm.h>
#include <ea/generational_models/death_birth_process.h>
#include <ea/selection/elitism.h>
#include <ea/datafiles/fitness.h>
#include <ea/cmdline_interface.h>
#include <ea/direct_markov_network.h>
#include "sunspot.h"
using namespace ealib;
using namespace mkv;

template <typename EA>
struct sunspot_configuration : public markov_network_configuration<EA> {
    typedef markov_network_configuration<EA> base_type;
    
    //! Called as the final step of EA initialization.
    virtual void initialize(EA& ea) {
        base_type::initialize(ea);
        
        std::size_t nbits = get<SUNSPOT_INTEGER_BITS>(ea) + get<SUNSPOT_FRACTIONAL_BITS>(ea);
        
        // we're currently limiting the number of bits that we're using to == long:
        assert(nbits < (sizeof(long)*8));
        
        // make sure we have at least 1 lag:
        assert(get<SUNSPOT_INPUT_LAGS>(ea) >= 1);
        
        boost::get<mkv::markov_network::IN>(base_type::mkv_desc) = get<SUNSPOT_INPUT_LAGS>(ea) * nbits;
        boost::get<mkv::markov_network::OUT>(base_type::mkv_desc) = 2 * nbits * get<SUNSPOT_PREDICTION_HORIZON>(ea);
    }
};


typedef evolutionary_algorithm<
mkv::representation_type,
mkv::mutation_type,
sunspot_fitness,
sunspot_configuration,
recombination::asexual,
generational_models::death_birth_process<selection::proportionate< >, selection::elitism<selection::random> >
> ea_type;


/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        // markov network options
        add_option<MKV_DESC>(this);
        add_option<MKV_UPDATE_N>(this);
        add_option<MKV_GATE_TYPES>(this);
        add_option<MKV_INITIAL_GATES>(this);
        add_option<MKV_MAX_GATES>(this);
        add_option<MKV_MIN_GATES>(this);
        add_option<GATE_INPUT_LIMIT>(this);
        add_option<GATE_INPUT_FLOOR>(this);
        add_option<GATE_OUTPUT_LIMIT>(this);
        add_option<GATE_OUTPUT_FLOOR>(this);
        add_option<GATE_HISTORY_LIMIT>(this);
        add_option<GATE_HISTORY_FLOOR>(this);
        add_option<GATE_WV_STEPS>(this);
        
        // ea options
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DUPLICATION_P>(this);
        
        add_option<REPLACEMENT_RATE_P>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<ELITISM_N>(this);
        
        // sunspot options
        add_option<SUNSPOT_TRAIN>(this);
        add_option<SUNSPOT_TEST>(this);
        add_option<SUNSPOT_INTEGER_BITS>(this);
        add_option<SUNSPOT_FRACTIONAL_BITS>(this);
        add_option<SUNSPOT_PREDICTION_HORIZON>(this);
        add_option<SUNSPOT_LIMIT>(this);
        add_option<SUNSPOT_INPUT_LAGS>(this);
    }
    
    virtual void gather_tools() {
        add_tool<sunspot_data>(this);
        add_tool<sunspot_train_detail>(this);
        add_tool<sunspot_test_detail>(this);
        add_tool<sunspot_train_rmse>(this);
        add_tool<sunspot_test_rmse>(this);
        add_tool<sunspot_train_predictions>(this);
        add_tool<sunspot_test_predictions>(this);
        add_tool<mkv::causal_graph>(this);
        add_tool<mkv::reduced_graph>(this);
        add_tool<mkv::genetic_graph>(this);
        add_tool<mkv::network_statistics>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness>(this, ea);
    }
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
