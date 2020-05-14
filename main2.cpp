#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <chrono>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <cfloat>

using namespace std;

struct PathToGoal {
    vector<int> route;
    long double length{0};

    static void printPathTaken(PathToGoal path) {
        /// Useful for debugging
        for (int i = path.route.size() - 1; i >= 0; i--) {
            cout << path.route[i] << " ";
        }
        cout << endl;
    }

    PathToGoal operator+=(PathToGoal const &rhs) {
        /// Note: Only adds the path to the end of the +=
        /// if the last node is equal to the first node in the
        /// next path
        if (rhs.route[0] == this->route[this->route.size() - 1]) {
            this->length += rhs.length;
            for (int i = 1; i < rhs.route.size(); i++) {
                // since the first one is the same we skip it
                this->route.push_back(rhs.route[i]);
            }
            return *this;
        } else {
            cerr << "Attempted to connect unconnectable paths" << endl;
            return PathToGoal{};
        }
    }

    bool operator==(PathToGoal const &other) {
        // If the lengths are exactly equal they're *probably* the same paths.
        // Huge speed boost for a risk
        return this->length == other.length;
    }

};

struct OptimalPathMap {
    OptimalPathMap() = default;

    OptimalPathMap(int size) : nextNode{vector<int>(size, -1)},
                               lengthToGoal{vector<double>(size, LDBL_MIN)} {}

    vector<int> nextNode;
    vector<double> lengthToGoal;
};

struct shortestPathLengthComparator {
    bool operator()(const pair<int, long double> &p1,
                    const pair<int, long double> &p2) {
        return p1.second > p2.second;
    }
};

struct shortestPathComparator {
    bool operator()(PathToGoal const &p1, PathToGoal const &p2) {
        return p1.length < p2.length;
    }
};

long double timeSince(
        chrono::time_point<std::chrono::high_resolution_clock> start) {

    long double time =
            chrono::duration_cast<chrono::nanoseconds>
                    (chrono::high_resolution_clock::now() - start).count();
    time /= 1000000000;
    return time;
}

class Searcher {
public:
    explicit Searcher(string const &filePath) {
        int numberOfNodes;
        int numberOfEdges;
        ifstream ifs(filePath);
        ifs >> numberOfNodes >> numberOfEdges;

        for (int i = 0; i < numberOfNodes; i++) {
            graph.emplace_back();
        }

        // Load in graph information
        for (int i = 0; i < numberOfEdges; i++) {
            int from, to;
            long double length;
            ifs >> from >> to >> length;
            graph[from][to] = length;
//            graph[from].insert({to, length}); // use this to find the 1035 number
        }
        // Grab the source and goal
        ifs >> source >> goal >> K;
        deadSpurs = vector<bool>(graph.size(), false);
        oPM = OptimalPathMap(numberOfNodes);
        fillOPM(); // just do this straight away
    }

    PathToGoal buildPath(int fromNode) {
        /// uses a reverse node map to build a path from the optimal paths
        PathToGoal path;
        path.route.push_back(fromNode);
        path.length = oPM.lengthToGoal[fromNode];
        if (path.length == DBL_MAX) {
            cerr << "Attempted to build an invalid path" << endl;
        }
        while (path.route[path.route.size() - 1] != goal) {
            path.route.push_back(oPM.nextNode[path.route[path.route.size() - 1]]);
        }
        return path;
    }

    //=======================Searching Algorithms=============================

    void fillOPM() {
        /// Generate all the lengths to the goal and optimal paths to there.
        /// Essentially complete Dijkstra's
        priority_queue<pair<int, long double>,
                vector<pair<int, long double>>,
                shortestPathLengthComparator> pQueue{};
        pQueue.push({goal, 0});
        oPM.lengthToGoal[goal] = 0;
        oPM.nextNode[goal] = goal; // just make goal loop into itself since
        // it is the final goal

        while (!pQueue.empty()) {
            auto topPath = pQueue.top();
            pQueue.pop();
            for (auto &vertex : graph[topPath.first]) {
                // Update it if it improve the grid
                if (oPM.nextNode[vertex.first] != -1 &&
                    oPM.lengthToGoal[vertex.first] <
                    oPM.lengthToGoal[topPath.first] + vertex.second) {
                    continue;
                }

                oPM.lengthToGoal[vertex.first] = oPM.lengthToGoal[topPath.first]
                                                 + vertex.second;

                oPM.nextNode[vertex.first] = topPath.first;

                pQueue.push(pair<int, long double>{
                        vertex.first, oPM.lengthToGoal[vertex.first]
                });
            }
        }
    }

    PathToGoal generatePathToGoal(int from, unordered_set<int> &bannedNodes, unordered_set<int> &bannedInitialEdges) {
        // since we have the oPM we actually just need to find one valid move
        // outside of bannedNodes, then we can use oPM to fill in the rest.
        // Since the previous path won't contain the anything in the optimal
        // path, this should be good. We can, however, verify this case
        if (graph[from].empty()) {
            // an empty path will symbolise its impossible
            return PathToGoal{};
        }
        priority_queue<pair<int, double>,
                vector<pair<int, double>>,
                shortestPathLengthComparator> options;
        for (auto &vertex :  graph[from]) {
            // split each into its own if statement for readability
            if (oPM.lengthToGoal[vertex.first] == DBL_MAX) {
                // no path existed anyways
                continue;
            }
            // Could have merged these into a single set, however, this just makes
            // the outside more intuitive and allows us to replace this funciton with
            // a Dijkstra in the future as well
            if (bannedNodes.find(vertex.first) != bannedNodes.end()) {
                continue; // its a banned node
            }
            if (bannedInitialEdges.find(vertex.first) != bannedNodes.end()) {
                continue; // its a banned node
            }
            options.push({vertex.first, oPM.lengthToGoal[vertex.first] + graph[from][vertex.first]});
        }

        PathToGoal path;
        bool valid = false;
        while (!options.empty() && !valid) {
            // we'll first set the length once we know its valid
            path.route.clear();
            path.route.push_back(from);
            path.route.push_back(options.top().first);
            path.length = options.top().second;
            options.pop();

            valid = true;
            while (valid && path.route[path.route.size() - 1] != goal) {
                path.route.push_back(oPM.nextNode[path.route[path.route.size() - 1]]);
                if (bannedNodes.find(path.route[path.route.size() - 1]) != bannedNodes.end()) {
                    // if its a banned node this path is no longer valid
                    valid = false;
                }
                // since our length won't actually be the optimal length we
                // have to add the distance from the first node to the
                // second node
            }
        }
        return path;
    }

    //==========================TOOLS FOR YEN=================================

    PathToGoal slicePathUpToAndEqualTo(PathToGoal &path, int upToAndIncluding) {
        /// Slices a path from the start of it up to the end (excluding the end)
        /// Includes the end
        /// Mostly here for access to the graph
        PathToGoal newPath{};
        newPath.route.reserve(upToAndIncluding);
        for (int i = 0; i <= upToAndIncluding; i++) {
            newPath.route.push_back(path.route[i]);
            if (i != 0)
                newPath.length += graph[path.route[i - 1]][path.route[i]];
        }
        return newPath;
    }

    static bool
    checkEqualToRoot(vector<int> &rootPath, vector<int> &pathToCheck) {
        for (int i = 0; i < rootPath.size() &&
                        i < pathToCheck.size(); i++) {
            if (rootPath[i] != pathToCheck[i]) return false;
        }
        return true;
    }

    void
    removeEdgeFromGraph(int from, int to, vector<tuple<int, int, long double>>
    &graphChanges) {
        auto iter = graph[from].find(to);
        // occasionally we end up trying to remove the same edge twice, so we
        // have to prevent that from happening
        if (iter != graph[from].end()) {
            graphChanges.emplace_back(from, to, iter->second);
            graph[from].erase(iter);
            removeEdgeFromGraph(to, from, graphChanges); // remove the
            // opposite direction as well
        }
    }

    // =======================YEN VARIATIONS==================================
    vector<PathToGoal> yen(int numberOfPaths) {
        /// if no function pointer is passed the basic version is ran
        using namespace std::placeholders;
        auto fp = bind(&Searcher::getPotentialPaths_basic, *this, _1, _2);
        return yen(numberOfPaths, fp);
    }

    vector<PathToGoal> yen() {
        /// if no function pointer is passed the basic version is ran
        using namespace std::placeholders;
        auto fp = bind(&Searcher::getPotentialPaths_basic, *this, _1, _2);
        return yen(K, fp);
    }

    vector<PathToGoal> yenApproximation() {
        using namespace std::placeholders;
        auto fp = bind(&Searcher::getPotentialPaths_random, *this, _1, _2);
        sorted = false;
        return yen(K, fp);
    }

    vector<PathToGoal> yenApproximation(int numberOfPaths) {
        using namespace std::placeholders;
        auto fp = bind(&Searcher::getPotentialPaths_random, *this, _1, _2);
        sorted = false;
        return yen(numberOfPaths, fp);
    }

    vector<PathToGoal> yen(int numberOfPaths,
                           std::function<set<PathToGoal, shortestPathComparator> &
                                   (vector<PathToGoal> &,
                                    set<PathToGoal, shortestPathComparator> &, int)>
                           &&fun) {
        auto startTime = chrono::high_resolution_clock::now();
        /// Returns a vector of k shortest loopless paths
        vector<PathToGoal> confirmedPaths{buildPath(source)};
        // since set has an O(1) lookup and keeps the highest element at the
        // front, its the easiest to use. An unordered set is another option,
        // but requires us to define a custom hash function
        set<PathToGoal, shortestPathComparator> potentialPaths;
        // Instead of fully recovering the graph, we can just 'recover' it by
        // recording what we've done so far to it

        // since the first one was already found we start from the second x
        for (int i = 0; i < numberOfPaths - 1; i++) {
            // The possible nodes go from the second node up to the second
            // last of the previously found shortest path
            potentialPaths = fun(confirmedPaths, potentialPaths, i);

            if (potentialPaths.empty()) {
                cout << "ran out of paths" << endl;
                break;
            }

            // our potential paths are in a set so they should already be sorted
            confirmedPaths.push_back(*potentialPaths.begin());
            potentialPaths.erase(potentialPaths.begin());
        }

        if (!sorted) // if the approximation algorithm was used the result wouldn't be right
            sort(confirmedPaths.begin(), confirmedPaths.end(), shortestPathComparator());

        // For some reason we aren't finding paths in order
        int k = 0;
        cout << "Took " << timeSince(startTime) << " seconds" << endl;
        for (auto &p : confirmedPaths) {
            cout << "K-" << ++k << ":\t" << setprecision(20) << p.length << endl;
//            if (k >= 10) break; // option to only print top k if wanted
        }
        return confirmedPaths;
    }

    set<PathToGoal, shortestPathComparator> &getPotentialPaths_basic(
            vector<PathToGoal> &confirmedPaths, set<PathToGoal, shortestPathComparator> &potentialPaths) {
        /// This just runs the normal Yen's core
        /// which iterates through all the nodes in a given path and
        /// defines confirmedPaths
        PathToGoal &lastPath = confirmedPaths[confirmedPaths.size() - 1];
        for (int j = 0; j < lastPath.route.size() - 1; j++) {
            PathToGoal newPath = yenGenerateSpurPath(confirmedPaths, j);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
//                cout << j + 1 << '\t'
//                     << newPath.length << endl;
            }
        }
        return potentialPaths;
    }

    set<PathToGoal, shortestPathComparator> &getPotentialPaths_random(
            vector<PathToGoal> &confirmedPaths,
            set<PathToGoal, shortestPathComparator> &potentialPaths) {
        /// Instead of Yen's core, this function iterates through a
        /// fraction of the route randomly
        // TODO: actually make this function
        PathToGoal &lastPath = confirmedPaths[confirmedPaths.size() - 1];
        auto rng = default_random_engine(chrono::system_clock::now()
                                                 .time_since_epoch().count());
        vector<int> indexes;
        indexes.reserve(lastPath.route.size() - 1);
        for (int x = 0; x < lastPath.route.size() - 1; ++x) {
            indexes.push_back(x);
        }
        shuffle(indexes.begin(), indexes.end(), rng);

        int nodesToSolveFor = sqrt(lastPath.route.size());
        int pathIndex = 0;
        while (pathIndex < nodesToSolveFor && pathIndex < indexes.size()) {
            int index = indexes[pathIndex];
//            if (deadSpurs[lastPath.route[index]]) continue; // skip likely dead spurs
            PathToGoal newPath = yenGenerateSpurPath(confirmedPaths, pathIndex);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
            } else if (newPath.route[newPath.route.size() - 1] != goal) {
                // we've confirmed that we've fully explored this spur so
                // don't try anymore
                deadSpurs[lastPath.route[index]] = true;
                nodesToSolveFor++; // since this spur didn't count we add another attempt
            }
            pathIndex++;
        }
        return potentialPaths;
    }

    PathToGoal yenGenerateSpurPath(vector<PathToGoal> &confirmedPaths,
                                   int spurNodeIndex) {
        /// Generates a spur path by cutting off the last confirmed path up to that
        /// node, and then

        // Side note about this function: Yes, I could generate both the banned
        // edges and nodes as a single variable, however, making them into distinct
        // variables makes it more flexible
        PathToGoal &lastPath = confirmedPaths[confirmedPaths.size() - 1];
        int spurNode = lastPath.route[spurNodeIndex];
        PathToGoal basePath = slicePathUpToAndEqualTo(lastPath, spurNodeIndex);
        unordered_set<int> bannedStartingEdges;
        unordered_set<int> bannedNodes;

        // Filter out future nodes that were explored
        for (auto &p : confirmedPaths) {
            if (checkEqualToRoot(basePath.route, p.route))
                bannedStartingEdges.insert(p.route[spurNodeIndex + 1]);
        }

        // Filter out the base path
        bannedNodes = unordered_set<int>(basePath.route.begin(), basePath.route.end());

        // Find the path and combine it with the base
        PathToGoal spurPath = generatePathToGoal(
                spurNode, bannedNodes, bannedStartingEdges);
        if (!spurPath.route.empty()
            && spurPath.route[spurPath.route.size() - 1] == goal) {
            // if we actually found a valid path
            basePath += spurPath;
        } else {
            basePath.length = numeric_limits<long double>::max();
        }
        return basePath;
    }

private:
    // contains the int x then the distance from the node
    vector<unordered_map<int, long double>> graph;
    vector<bool> deadSpurs;
    OptimalPathMap oPM;
    int source{}, goal{}, K{};
    bool sorted{true}; // this gets toggled if the approximation algorithm is used
};

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        cout << "Error: an input and an output file is required."
             << endl
             << "Call the program like `program.exe input.txt`"
             << endl;
    }
    Searcher s{argv[1]};
    cout << "-----------------True values-----------------------" << endl;
    s.yen();
    cout << "----------------Approximation-----------------------" << endl;
    s.yenApproximation();
    return 0;
}
