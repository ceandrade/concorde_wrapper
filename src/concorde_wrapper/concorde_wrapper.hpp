/******************************************************************************
 * concorde_wrapper.hpp: Interface for Concorde Functions Wrapper.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 07, 2012 by andrade
 * Last update: Jul 23, 2012 by andrade
 *
 * This code is released under LICENSE.md.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#ifndef CONCORDE_WRAPPER_HPP
#define CONCORDE_WRAPPER_HPP

// Concorde inclusions
extern "C" {
    #include "machdefs.h"
    #include "util.h"
    #include "linkern.h"
    #include "kdtree.h"
    #include "edgegen.h"
    #include "macrorus.h"
}

#include <vector>
#include <lemon/full_graph.h>
#include <lemon/dim2.h>
using namespace lemon;

/**
 * \brief Concorde Functions Wrapper.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date Jun 16, 2012
 *
 * This class contains wrappers for the Concorde Lin Kernighan functions.
 * Can be use in multi-thread programs that use OpenMP.
 */
class ConcordeWrapper {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor.
         * \param seed the seed to be used by Concorde.
         * \param max_threads Used to instantiate thread-safe data.
         * \param graph the full graph where will be extracted subcycles.
         * \param dist the edges weigts.
         * \param time_bound Time limit to linkern run.
         * \param length_bound stops when linkern reaches a cycle such cost <= length_bound.
         * \param kick_type Initialization method to linkern. See linkern.h for more details.
         * \param stallcount The max number of 4-swaps without progress.
         */
        ConcordeWrapper(const int seed = 0,
                        const unsigned max_threads = 1,
                        const FullGraph *graph = NULL,
                        const FullGraph::EdgeMap<int> *dist = NULL,
                        const double time_bound = -1.0,
                        const double length_bound = -1.0,
                        const int kick_type = CC_LK_WALK_KICK,
                        const int stallcount = 1000);

        /** \brief Destructor. */
        ~ConcordeWrapper();
        //@}

        /** \name Load and settings methods. */
        //@{
        /** \brief Set a new graph to wrapper.
         * \param graph a reference to the graph.
         * \param dist a reference to edge weights.
         */
        void updateData(const FullGraph *graph,
                        const FullGraph::EdgeMap<int> *dist);

        /** \brief Set the number of threads and rebuild the data structures.
         * \param max_threads Maximum number of threads.
         */
        void setThreadsNum(const unsigned max_threads);

        /// Get the number of threads.
        unsigned getThreadsNum();

        /** \brief Load a TSP Lib instance to the graph.
         * \param[in] instance_file the file to be loaded.
         * \param[out] graph the graph loaded.
         * \param[out] dist the edge weights loaded.
         * \param[out] coords the 2D coordinates.
         * \param[out] loaded_coords it will be set to false if the instance has
         * not geometric data. In this case, the loaded_coords map will be empty.
         * It will set to true otherwise.
         *
         * \warning The caller MUST destroy _graph and _dist, to free memory.
         */
        void loadTSP(char *instance_file,
                     FullGraph *graph,
                     FullGraph::EdgeMap<int> *dist,
                     FullGraph::NodeMap< dim2::Point<double> > *coords,
                     bool *loaded_coords);
        //@}

        /**
         * \brief From the list of vertices, return a improved cycle using
         * the Lin Kernighan algorithm.
         *
         * \param[in] nodes A list of nodes to try to obtain a good cycle. The
         *                  order of these nodes defines the initial cycle.
         * \param[out] cycle A list of nodes representing the cycle.
         * \return the cycle value.
         */
        double getLinKernTour(const std::vector<lemon::FullGraph::Node> &nodes,
                              std::vector<lemon::FullGraph::Node> *cycle);

    public:
        /** \name External constant data */
        //@{
        /// A reference to the full graph from where we extracted the subcycles.
        const lemon::FullGraph *graph;

        /// A reference to te edges cost.
        const lemon::FullGraph::EdgeMap<int> *dist;
        //@}

        /** \name Lin Kernighan options.
         * See linkern.c for more details.
         * */
        //@{
        /// Time limit to linkern run.
        double time_bound;

        /// Stops when linkern reaches a cycle such cost <= length_bound.
        double length_bound;

        /// Initialization method to linkern.
        int kick_type;

        /// The max number of 4-swaps without progress.
        int stallcount;
        //@}

    protected:
        /** \name Update methods. */
        //@{
        /** \brief Clean up Concorde data structures. */
        void cleanupData();
        //@}

    protected:
        /// Hold the state of the random generator from Concorde.
        CCrandstate rstate;

        /** \name Concorde data structures for each thread. */
        //@{
        /// Hold the graph data and other stuff.
        std::vector<CCdatagroup*> concorde_data;

        /// Hold the initial cycles.
        std::vector<int*> incycles;

        /// Hold the result cycles.
        std::vector<int*> outcycles;

        /// Hold the temp cycles.
        std::vector<int*> templists;
        //@}

        /// Verify the wrapper status.
        bool init;

        /// Verify the wrapper status.
        unsigned max_threads;
};

#endif // CONCORDE_WRAPPER_HPP

