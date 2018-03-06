/******************************************************************************
 * concorde_wrapper.cpp: Implementation for Concorde Functions Wrapper.
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

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <lemon/lgf_writer.h>
#endif

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include "concorde_wrapper.hpp"
using namespace lemon;
using namespace std;

//--------------------------[ Default Constructor ]---------------------------//

ConcordeWrapper::ConcordeWrapper(const int seed, const unsigned _max_threads,
    const FullGraph *_graph, const FullGraph::EdgeMap<int> *_dist,
    const double _time_bound, const double _length_bound, const int _kick_type,
    const int _stallcount):

    graph(_graph), dist(_dist),
    time_bound(_time_bound), length_bound(_length_bound),
    kick_type(_kick_type), stallcount(_stallcount),

    rstate(),
    concorde_data(_max_threads, (CCdatagroup*)NULL),
    incycles(_max_threads, (int*)NULL),
    outcycles(_max_threads, (int*)NULL),
    templists(_max_threads, (int*)NULL),
    init(false),
    max_threads(_max_threads)
{
    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Building the wrapper"
         << endl;

    if(_graph != NULL)
        cout << "> Graph size: " << graph->nodeNum() << "\n";
    else
        cout << "> Graph: null\n";

    cout << "> Seed: " << seed << "\n"
         << "> # of threads: " << max_threads
         << endl;
    #endif

    if(seed != 0)
        CCutil_sprand(seed, &rstate);
    else
        CCutil_sprand((int)CCutil_real_zeit(), &rstate);

    if(graph != NULL && dist != NULL)
        updateData(graph, dist);

    #ifdef DEBUG
    cout << "------------------------------------------------------" << endl;
    #endif
}

//------------------------------[ Destructor ]--------------------------------//

ConcordeWrapper::~ConcordeWrapper() {
    cleanupData();
}

//-----------------------------[ Update Data ]--------------------------------//

void ConcordeWrapper::updateData(const FullGraph *_graph,
                                 const FullGraph::EdgeMap<int> *_dist) {

    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Updating graph data on wrapper"
         << endl;
    #endif

    if(_graph == NULL || _dist == NULL) {
        stringstream error_msg;
        error_msg << "\n*** Error on updating data"
                  << " at " << __PRETTY_FUNCTION__
                  << "\n*** _graph and _dist MUST BE valid references";
        throw runtime_error(error_msg.str());
    }

    // If init, clean up first.
    if(init)
        cleanupData();

    graph = _graph;
    dist = _dist;

    // Initializing Concorde structures
    const int ncount = graph->nodeNum();

    concorde_data.resize(max_threads, (CCdatagroup*)NULL);
    incycles.resize(max_threads, (int*)NULL);
    outcycles.resize(max_threads, (int*)NULL);
    templists.resize(max_threads, (int*)NULL);

    for(unsigned i = 0; i < max_threads; ++i) {
        CCdatagroup *dat = new CCdatagroup();
        CCutil_init_datagroup(dat);

        // First, we set the norm to CC_MATRIXNORM because
        // we will pass the edge lengths explicitly.
        CCutil_dat_setnorm(dat, CC_MATRIXNORM);

        // Allocating space in dat.
        dat->adj = CC_SAFE_MALLOC(ncount, int*);
        dat->adjspace = CC_SAFE_MALLOC((ncount) * (ncount+1)/2, int);

        concorde_data[i] = dat;

        // Allocate the initial and result cycles.
        incycles[i] = CC_SAFE_MALLOC(ncount, int);
        outcycles[i] = CC_SAFE_MALLOC(ncount, int);
        templists[i] = CC_SAFE_MALLOC(2 * graph->edgeNum(), int);
    }

    init = true;

    #ifdef DEBUG
    cout << "------------------------------------------------------" << endl;
    #endif
}

//----------------------------[ Clean Up Data ]-------------------------------//

void ConcordeWrapper::cleanupData() {
    if(init) {
        for(unsigned i = 0; i < concorde_data.size(); ++i) {
            CCdatagroup *dat = concorde_data[i];
            CCutil_freedatagroup(dat);
            delete dat;

            CC_IFFREE(incycles[i], int);
            CC_IFFREE(outcycles[i], int);
            CC_IFFREE(templists[i], int);
        }
    }

    concorde_data.clear();
    incycles.clear();
    outcycles.clear();
    templists.clear();

    init = false;
}

//-------------------------[ Set Threads Number ]-----------------------------//

void ConcordeWrapper::setThreadsNum(const unsigned _max_threads) {
    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Setting the number of threads to "
         << (_max_threads > 0? "_max_threads" : "1" )
         << endl;
    #endif

    if(_max_threads == 0)
        max_threads = 1;
    else
        max_threads = _max_threads;

    if(init)
        cleanupData();

    this->updateData(graph, dist);
}

//-------------------------[ Get Threads Number ]-----------------------------//

unsigned ConcordeWrapper::getThreadsNum() {
    return this->max_threads;
}

//------------------------------[ Load TSP ]----------------------------------//

void ConcordeWrapper::loadTSP(char *instance_file,
             FullGraph *_graph, FullGraph::EdgeMap<int> *_dist,
             FullGraph::NodeMap< dim2::Point<double> > *_coords, bool *loaded_coords) {
    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Loading TSP Lib Instance: " << instance_file
         << endl;
    #endif

    if(_graph == NULL || _dist == NULL) {
        stringstream error_msg;
        error_msg << "\n*** Error on trying load TSPLIB instance "
                  << instance_file
                  << " at " << __PRETTY_FUNCTION__
                  << "\n*** _graph and _dist MUST BE valid references";
        throw runtime_error(error_msg.str());
    }

    // Fist, we need to create Concorde structures to be loaded.
    CCdatagroup dat;
    int ncount;
    CCutil_init_datagroup(&dat);

    // Load TSPLin instance file.
    if(CCutil_gettsplib(instance_file, &ncount, &dat)) {
        CCutil_freedatagroup(&dat);
        stringstream error_msg;
        error_msg << "\n*** Error on trying load TSPLIB instance "
                  << instance_file
                  << " at " << __PRETTY_FUNCTION__
                  << " \n*** File couldn't be loaded ***";
        throw runtime_error(error_msg.str());
    }

    #ifdef FULLDEBUG
    cout << "\nConcorde loaded data: "
         << "\n- Number of vertices: " << ncount
         << endl;
    #endif

    _graph->resize(ncount);

    // Copy edge lengths.
    for(int i = 0; i < ncount; ++i)
        for(int j = i+1; j < ncount; ++j)
            (*_dist)[_graph->edge(graph->nodeFromId(i), graph->nodeFromId(j))] = CCutil_dat_edgelen(i, j, &dat);

    // Copy the coordinates, if exist.
    if((dat.norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE ||
       (dat.norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        for(int i = 0; i < ncount; ++i)
            (*_coords)[_graph->nodeFromId(i)] = dim2::Point<double>(dat.x[i], dat.y[i]);

        *loaded_coords = true;
    }
    else if((dat.norm & CC_MATRIX_NORM_SIZE) == CC_MATRIX_NORM_SIZE) {
        ifstream file(instance_file, ios::in);

        // Locate DISPLAY_DATA_SECTION
        string line;
        while(!(line == "DISPLAY_DATA_SECTION" || line == "EOF")) {
            getline(file, line);
            line.erase(line.find_last_not_of(" \n\r\t")+1);
        }

        if(line == "DISPLAY_DATA_SECTION") {
            int node;
            double x, y;
            for(int i = 0; i < ncount; ++i) {
                file >> node >> x >> y;
                (*_coords)[_graph->nodeFromId(i)] = dim2::Point<double>(x, y);
            }
            *loaded_coords = true;
        }
        else {
            *loaded_coords = false;
        }

        file.close();
    }
    // What is this norm? If we cannot read, then we don't load the coordinates.
    else
        *loaded_coords = false;

    // Clean up
    CCutil_freedatagroup(&dat);

    #ifdef FULLDEBUG
    cout << "\nLoaded graph:\n";
    graphWriter(*_graph)
        .edgeMap("Cost", *_dist)
        .nodeMap("Coord", *_coords)
        .run();

    if(!loaded_coords)
        cout << "\n\n This instance has not coordinate/geometric data.";

    cout << "------------------------------------------------------" << endl;
    #endif
}

//-----------------------[ Get Lin Kernighan Tour ]---------------------------//

double ConcordeWrapper::getLinKernTour(const vector<FullGraph::Node> &nodes,
                                      vector<FullGraph::Node> *cycle) {
    #ifdef DEBUG
    {
    cout << "\n------------------------------------------------------\n"
         << "> Get Lin Kernighan Tour";

    #ifndef _OPENMP
    cout << " (thread 0)";
    #else
    cout << " (thread " << omp_get_thread_num() << ")";
    #endif
    #endif

    if(!init) {
        stringstream error_msg;
        error_msg << "\n*** Error"
                  << " at " << __PRETTY_FUNCTION__
                  << "\n*** Concorde Wrapper MUST HAVE a valid graph"
                  << " (use the updateData() at least once before call getLinKernTour())";
        throw runtime_error(error_msg.str());
    }

    #ifdef DEBUG
    cout << "\n> Original cycle: ";
    double cost = 0.0;
    vector<FullGraph::Node>::const_iterator it1, it2;
    for(it2 = nodes.begin(); it2 != nodes.end(); ) {
        it1 = it2++;
        cout << graph->id(*it1) << " ";

        if(it2 != nodes.end())
            cost += (*dist)[graph->edge(*it1,*it2)];
    }
    cost += (*dist)[graph->edge(*it1,*(nodes.begin()))];

    cout << "\n> Cost: " << cost
         << endl;
    }
    #endif

    if(cycle == NULL) {
        stringstream error_msg;
        error_msg << "\n*** Error "
                  << " at " << __PRETTY_FUNCTION__
                  << "\n*** cycle MUST BE valid reference";
        throw runtime_error(error_msg.str());
    }

    if(nodes.size() < 4) {
        cycle->clear();
        cycle->reserve(nodes.size());
        double cycle_value = 0.0;

        vector<FullGraph::Node>::const_iterator it1, it2;
        for(it2 = nodes.begin(); it2 != nodes.end(); ) {
            cycle->push_back(*it2);
            it1 = it2++;
            if(it2 != nodes.end())
                cycle_value += (*dist)[graph->edge(*it1,*it2)];
        }
        cycle_value += (*dist)[graph->edge(*it1,*(nodes.begin()))];

        #ifdef DEBUG
        cout << "> Nothing to do. Return the original cycle"
             << "\n------------------------------------------------------" << endl;
        #endif

        return cycle_value;
    }

    #ifndef _OPENMP
    const unsigned thread = 0;
    #else
    const unsigned thread = omp_get_thread_num();;
    #endif

    // Set data structures according to the current thread.
    CCdatagroup *dat = concorde_data[thread];
    int *incycle = incycles[thread];
    int *outcycle = outcycles[thread];
    int *templist = templists[thread];

    const int ncount = nodes.size();

    // Setting the pointers of the edges in the matrix.
    for(int i = 0, j = 0; i < ncount; ++i) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }

    // Load the data in Concorde. We use the MATRIX LOWER DIAG ROW format.
    for(int i = 0; i < ncount; ++i)
        for (int j = 0; j <= i; ++j) {
            if(i != j)
                dat->adj[i][j] = (int) (*dist)[graph->edge(nodes[i],nodes[j])];
            else
                dat->adj[i][j] = 0;
        }

    #ifdef FULLDEBUG
    cout << "\nConcorde loaded data: ";
    copy(dat->adjspace, dat->adjspace + ((ncount) * (ncount+1)/2),
         ostream_iterator<int>(cout, " "));

    cout << endl;
    #endif

    // Create the initial cycle. Simple put the vertices in the receiving order.
    for(int i = 0; i < ncount; ++i)
        incycle[i] = i;


    #ifdef DEBUG
    const int run_silently = 0;
    #else
    const int run_silently = 1;
    #endif

    int tempcount = ((ncount) * (ncount-1)/2);

    // This functions is nice when we divide the graph by regions.
    // ATTENTION: using templist isn't thread safe here. To do it, we need modify
    // the templist initialization, i.e., build and destroy then in this function.
    // See the constructor for more details.
//     int nearnum = 8;
//     CCedgegen_junk_k_nearest (ncount, nearnum, dat, (double *) NULL, 1,
//                              &tempcount, &templist, run_silently);

    // If we don't use some function like above, we must build the templist by hand.
    for(int j = 0, count = 0; j < ncount; ++j)
        for(int k = j+1; k < ncount; ++k) {
            templist[count++] = j;
            templist[count++] = k;
        }


     #ifdef FULLDEBUG
    cout << "\nTemp Edges: ";
    copy(templist, templist + (tempcount * 2), ostream_iterator<int>(cout, " "));
    cout << endl;
    #endif

    // Run Lin Kernighan
    double cycle_value = 0.0;
    int repeatcount = 0; //(ncount < 10? 0 : ncount);

    int result = CClinkern_tour(ncount, dat, tempcount, templist, stallcount,
                                repeatcount, incycle, outcycle, &cycle_value,
                                run_silently, time_bound, length_bound,
                                (char *) NULL, kick_type, &rstate);
    if(result) {
        throw runtime_error("\n*** Somethig is wrong with CClinkern_tour() ***\n");
    }

    cycle->clear();
    cycle->reserve(ncount);
    for(int i = 0; i < ncount; ++i)
        cycle->push_back(nodes[outcycle[i]]);

    #ifdef DEBUG
    cout << "\n\n> Cycle found:"
         << "\n>> Concorde numbering: ";
    copy(outcycle, outcycle+ncount, ostream_iterator<int>(cout, " "));
    cout << "\n|| Cost: " << cycle_value;

    cout << "\n>> Original numbering: ";
    {
        double cost = 0.0;
        vector<FullGraph::Node>::const_iterator it1, it2;
        for(it2 = cycle->begin(); it2 != cycle->end(); ) {
            it1 = it2++;
            cout << graph->id(*it1) << " ";

            if(it2 != cycle->end())
                cost += (*dist)[graph->edge(*it1,*it2)];
        }
        cost += (*dist)[graph->edge(*it1,*(cycle->begin()))];

        cout << "\n|| Cost: " << cost;
    }
    cout << "\n------------------------------------------------------" << endl;
    #endif

    return cycle_value;
}

