/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/autodiff/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/core/utility/share_obj.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>


namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    std::cout << "\n================    Test program for fully implicit three-phase black-oil flow     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    if (!use_deck) {
        OPM_THROW(std::runtime_error, "This program must be run with an input deck. "
                  "Specify the deck with deck_filename=deckname.data (for example).");
    }
    std::shared_ptr<GridManager> grid;
    std::shared_ptr<BlackoilPropertiesInterface> props;
    std::shared_ptr<BlackoilPropsAdInterface> new_props;
    std::shared_ptr<RockCompressibility> rock_comp;
    BlackoilState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    std::string deck_filename = param.get<std::string>("deck_filename");

#define USE_NEW_PARSER 1
#if USE_NEW_PARSER
    Opm::ParserPtr newParser(new Opm::Parser() );
    Opm::DeckConstPtr newParserDeck = newParser->parseFile( deck_filename );
#else
    std::shared_ptr<EclipseGridParser> deck;
    deck.reset(new EclipseGridParser(deck_filename));
#endif

    // Grid init
#if USE_NEW_PARSER
    grid.reset(new GridManager(newParserDeck));
#else
    grid.reset(new GridManager(*deck));
#endif

#if USE_NEW_PARSER
    Opm::EclipseWriter outputWriter(param, newParserDeck, share_obj(*grid->c_grid()));
#else
    Opm::EclipseWriter outputWriter(param, deck, share_obj(*grid->c_grid()));
#endif

    // Rock and fluid init
#if USE_NEW_PARSER
    props.reset(new BlackoilPropertiesFromDeck(newParserDeck, *grid->c_grid(), param));
    new_props.reset(new BlackoilPropsAdFromDeck(newParserDeck, *grid->c_grid()));
#else
    props.reset(new BlackoilPropertiesFromDeck(*deck, *grid->c_grid(), param));
    new_props.reset(new BlackoilPropsAdFromDeck(*deck, *grid->c_grid()));
#endif

    // check_well_controls = param.getDefault("check_well_controls", false);
    // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
    // Rock compressibility.
#if USE_NEW_PARSER
    rock_comp.reset(new RockCompressibility(newParserDeck));
#else
    rock_comp.reset(new RockCompressibility(*deck));
#endif
    // Gravity.
#if USE_NEW_PARSER
    gravity[2] = newParserDeck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;
#else
    gravity[2] = deck->hasField("NOGRAV") ? 0.0 : unit::gravity;
#endif
    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
        initBlackoilSurfvol(*grid->c_grid(), *props, state);
        enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
        const PhaseUsage pu = props->phaseUsage();
        if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
            const int np = props->numPhases();
            const int nc = grid->c_grid()->number_of_cells;
            for (int c = 0; c < nc; ++c) {
                state.gasoilratio()[c] = state.surfacevol()[c*np + pu.phase_pos[Gas]]
                    / state.surfacevol()[c*np + pu.phase_pos[Oil]];
            }
        }
    } else {
#if USE_NEW_PARSER
        initBlackoilStateFromDeck(*grid->c_grid(), *props, newParserDeck, gravity[2], state);
#else
        initBlackoilStateFromDeck(*grid->c_grid(), *props, *deck, gravity[2], state);
#endif
    }

    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::ofstream outStream;
    std::string output_dir;
    if (output) {
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        std::string filename = output_dir + "/timing.param";
        outStream.open(filename.c_str(), std::fstream::trunc | std::fstream::out);
        // open file to clean it. The file is appended to in SimulatorTwophase
        filename = output_dir + "/step_timing.param";
        std::fstream step_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
        step_os.close();
        param.writeParam(output_dir + "/simulation.param");
    }

#if USE_NEW_PARSER
    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;

    WellStateFullyImplicitBlackoil well_state;
    Opm::TimeMapPtr timeMap(new Opm::TimeMap(newParserDeck));
    SimulatorTimer simtimer;
    std::shared_ptr<EclipseState> eclipseState(new EclipseState(newParserDeck));

    // initialize variables
    simtimer.init(timeMap, /*beginReportStepIdx=*/0, /*endReportStepIdx=*/0);

    SimulatorReport fullReport;
    for (size_t reportStepIdx = 0; reportStepIdx < timeMap->numTimesteps(); ++reportStepIdx) {
        // Report on start of a report step.
        std::cout << "\n"
                  << "---------------------------------------------------------------\n"
                  << "--------------    Starting report step " << reportStepIdx << "    --------------\n"
                  << "---------------------------------------------------------------\n"
                  << "\n";

        WellsManager wells(eclipseState,
                           reportStepIdx,
                           *grid->c_grid(),
                           props->permeability());

        if (reportStepIdx == 0) {
            // @@@ HACK: we should really make a new well state and
            // properly transfer old well state to it every epoch,
            // since number of wells may change etc.
            well_state.init(wells.c_wells(), state);
        }

        simtimer.init(timeMap,
                      /*beginReportStepIdx=*/reportStepIdx,
                      /*endReportStepIdx=*/reportStepIdx + 1);

        if (reportStepIdx == 0)
            outputWriter.writeInit(simtimer, state, well_state.basicWellState());

        // Create and run simulator.
        SimulatorFullyImplicitBlackoil simulator(param,
                                                 *grid->c_grid(),
                                                 *new_props,
                                                 rock_comp->isActive() ? rock_comp.get() : 0,
                                                 wells,
                                                 linsolver,
                                                 grav);
        SimulatorReport episodeReport = simulator.run(simtimer, state, well_state);

        outputWriter.writeTimeStep(simtimer, state, well_state.basicWellState());
        fullReport += episodeReport;
    }

    std::cout << "\n\n================    End of simulation     ===============\n\n";
    fullReport.report(std::cout);

    if (output) {
        std::string filename = output_dir + "/walltime.param";
        std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
        fullReport.reportParam(tot_os);
        warnIfUnusedParams(param);
    }
#else
    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << "                        (number of epochs: "
              << (deck->numberOfEpochs()) << ")\n\n" << std::flush;

    SimulatorReport rep;
    // With a deck, we may have more epochs etc.
    WellStateFullyImplicitBlackoil well_state;
    int step = 0;
    SimulatorTimer simtimer;
    // Use timer for last epoch to obtain total time.
    deck->setCurrentEpoch(deck->numberOfEpochs() - 1);
    simtimer.init(*deck);

    const double total_time = simtimer.totalTime();
    for (int epoch = 0; epoch < deck->numberOfEpochs(); ++epoch) {
        // Set epoch index.
        deck->setCurrentEpoch(epoch);

        // Update the timer.
        if (deck->hasField("TSTEP")) {
            simtimer.init(*deck);
        } else {
            if (epoch != 0) {
                OPM_THROW(std::runtime_error, "No TSTEP in deck for epoch " << epoch);
            }
            simtimer.init(param);
        }
        simtimer.setCurrentStepNum(step);
        simtimer.setTotalTime(total_time);

        // Report on start of epoch.
        std::cout << "\n\n--------------    Starting epoch " << epoch << "    --------------"
                  << "\n                  (number of steps: "
                  << simtimer.numSteps() - step << ")\n\n" << std::flush;

        // Create new wells, well_state
        WellsManager wells(*deck, *grid->c_grid(), props->permeability());
        // @@@ HACK: we should really make a new well state and
        // properly transfer old well state to it every epoch,
        // since number of wells may change etc.
        if (epoch == 0) {
            well_state.init(wells.c_wells(), state);
        }

        if (epoch == 0)
            outputWriter.writeInit(simtimer, state, well_state.basicWellState());

        // Create and run simulator.
        SimulatorFullyImplicitBlackoil simulator(param,
                                                 *grid->c_grid(),
                                                 *new_props,
                                                 rock_comp->isActive() ? rock_comp.get() : 0,
                                                 wells,
                                                 linsolver,
                                                 grav);
        outputWriter.writeTimeStep(simtimer, state, well_state.basicWellState());

        if (epoch == 0) {
            warnIfUnusedParams(param);
        }
        SimulatorReport epoch_rep = simulator.run(simtimer, state, well_state);
        if (output) {
            epoch_rep.reportParam(outStream);
        }
        // Update total timing report and remember step number.
        rep += epoch_rep;
        step = simtimer.currentStepNum();
    }

    std::cout << "\n\n================    End of simulation     ===============\n\n";
    rep.report(std::cout);

    if (output) {
        std::string filename = output_dir + "/walltime.param";
        std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
        rep.reportParam(tot_os);
    }
#endif
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

