/******************************************************************************
 * test.cpp: Main Process.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 16, 2012 by andrade
 * Last update: Jun 16, 2012 by andrade
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cstdlib>
#include <vector>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <lemon/full_graph.h>
#include <lemon/lgf_writer.h>
#include <lemon/dim2.h>

using namespace std;
using namespace lemon;

#include "concorde_wrapper.hpp"

double printCycle(const vector<FullGraph::Node> &cycle, const FullGraph &graph,
                  const FullGraph::EdgeMap<int> &dist);

//-------------------------------[ Main ]------------------------------------//

int main(int argc, char* argv[]) {
    if(argc < 2) {
        cerr << "> Please, give the TSP instance file" << endl;
        return 1;
    }

    FullGraph graph;
    FullGraph::EdgeMap<int> dist(graph);
    FullGraph::NodeMap< dim2::Point<double> > coords(graph);
    bool loaded_coords;

    try {
        const unsigned long SEED = 152269023;
        const unsigned long MAX_THR = 2;

        ConcordeWrapper concorde(SEED, MAX_THR);
        concorde.loadTSP(argv[1], &graph, &dist, &coords, &loaded_coords);
        concorde.updateData(&graph, &dist);

        if(loaded_coords)
            for(FullGraph::NodeIt it(graph); it != INVALID; ++it)
                cout << "\n> " << graph.id(it) << ": "<< coords[it];
        else
            cout << "\nWe couldn't load the geometric information" << endl;

        vector<FullGraph::Node> nodes;
        vector<FullGraph::Node> *cycle = new vector<FullGraph::Node>();

        for(FullGraph::NodeIt it(graph); it != INVALID; ++it)
            nodes.push_back(it);

        concorde.getLinKernTour(nodes, cycle);

        double cost = 0.0;
        cout << "\n>> Cycle: ";
        cost = printCycle(nodes, graph, dist);
        cout << "\n|| Cost: " << cost << endl;

        delete cycle;
    }
    catch(exception& e) {
        cerr << "\n***********************************************************"
             << "\n****  Exception Occured: " << e.what()
             << "\n***********************************************************"
             << endl;
        return 70; // BSD software internal error code
    }
    return 0;
}

//----------------------------[ Print Cycle ]---------------------------------//

double printCycle(const vector<FullGraph::Node> &cycle, const FullGraph &graph,
                  const FullGraph::EdgeMap<int> &dist) {
    double cost = 0.0;
    vector<FullGraph::Node>::const_iterator it1, it2;

    for(it2 = cycle.begin(); it2 != cycle.end(); ) {
        it1 = it2++;
        cout << graph.id(*it1) << " ";

        if(it2 != cycle.end())
            cost += dist[graph.edge(*it1,*it2)];
    }
    cost += dist[graph.edge(*it1,*(cycle.begin()))];

    return cost;
}

