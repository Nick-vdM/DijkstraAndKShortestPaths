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

using namespace std;

struct Path {
    vector<int> route;
    double length{0};

    static void printPathTaken(Path path) {
        /// Useful for debugging
        for (int i = path.route.size() - 1; i >= 0; i--) {
            cout << path.route[i] << " ";
        }
        cout << endl;
    }

    Path operator+=(Path const &rhs) {
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
            return Path{};
        }
    }

    bool operator==(Path const &other) {
        return this->route == other.route;
    }

};

struct shortestPathLengthComparator {
    bool operator()(const pair<int, double> &p1, const pair<int, double> &p2) {
        return p1.second > p2.second;
    }
};

struct shortestPathComparator {
    bool operator()(Path const &p1, Path const &p2) {
        return p1.length < p2.length;
    }
};

double timeSince(
        chrono::time_point<std::chrono::high_resolution_clock> start) {

    double time =
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
            double length;
            ifs >> from >> to >> length;
            graph[from][to] = length;
        }
        // Grab the source and goal
        ifs >> source >> goal >> K;
    }

    // =======================GENERAL TOOLS=================================
    array<int, 3> setToLongPath(double maxSearchTime) {
        /// This function sets the source and goal to be as far apart as
        /// possible in the given amount of time
        // TODO: fix this? lol
        auto startClock = chrono::high_resolution_clock::now();
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> uni(0, 50);
        double timeSinceStart = 0;
        array<int, 3> startGoalLength{-1, -1, numeric_limits<int>::min()};
        // note: the number of nodes is more impactful
        // to the k-1 search statespace than the actual length of the path
        while (timeSinceStart < maxSearchTime) {
            array<int, 3> newStartGoalLength{-1, -1, -1};
            for(int i = 0; i < 2; i++) newStartGoalLength[i] = uni(rng);
            if (newStartGoalLength[0] == newStartGoalLength[1]) continue;
            newStartGoalLength[2] = singleGoalDijkstra(
                    newStartGoalLength[0],
                    newStartGoalLength[1]
            ).length;

            if (newStartGoalLength[3] > startGoalLength[3]) {
                startGoalLength = newStartGoalLength;
            }
            timeSinceStart = timeSince(startClock);
        }
        if (startGoalLength[0] != -1) {
            source = startGoalLength[0];
            goal = startGoalLength[1];
        }
        // in case if they're needed elsewhere, return it as well
        return startGoalLength;
    }

    // =======================TOOLS FOR DIJ==================================
    Path buildPath(vector<int> &nodeBefore, double &length, int localSource) {
        /// uses a reverse node map to build a path
        Path pathTaken;
        pathTaken.length = length;
        pathTaken.route.push_back(goal);
        pathTaken.route.reserve(nodeBefore.size());

        int currentNode = goal;
        while (currentNode != localSource) {
            currentNode = nodeBefore[currentNode];
            pathTaken.route.push_back(currentNode);
        }
        reverse(pathTaken.route.begin(), pathTaken.route.end());
        return pathTaken;
    }

    // ==============================DIJ======================================
    Path singleGoalDijkstra() {
        // without a specified function source it uses the class's source
        return singleGoalDijkstra(source, goal, unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource) {
        return singleGoalDijkstra(localSource, goal, unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource, int localGoal) {
        return singleGoalDijkstra(localSource, localGoal,
                                  unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource, unordered_set<int> bannedNodes) {
        return singleGoalDijkstra(localSource, goal, bannedNodes);
    }

    Path singleGoalDijkstra(int localSource, int localGoal, unordered_set<int>
    bannedNodes) {
        /// returns the path that was taken
        /// uses the class goal instead as unique goals are never used
        /// returns a large negative empty path if it fails to find it
        /// (double min value)

        priority_queue<pair<int, double>,
                vector<pair<int, double>>,
                shortestPathLengthComparator> pQueue{};
        pQueue.push(pair<int, double>{localSource, 0});
        int pushes = 0;
        // we could use our graph structure here, however, this has a faster
        // lookup (but still O(1) and uses O(n) space anyways. So we get a tiny
        // speed boost in exchange for a tiny amount of extra space, and we
        // are valuing the speed higher than the space
        vector<double> distance(graph.size(), numeric_limits<double>::max());
        distance[localSource] = 0;
        vector<int> nodeBefore(graph.size(), -1);
        double bestPathLength{numeric_limits<double>::min()};
        while (!pQueue.empty()) {
            auto topPath = pQueue.top();
            pQueue.pop();

            // we will always compare to the overall goal
            if (topPath.first == localGoal) {
                bestPathLength = topPath.second;
                break;
            }

            for (auto &vertex : graph[topPath.first]) {
                if (nodeBefore[vertex.first] != -1 ||
                    bannedNodes.find(vertex.first) != bannedNodes.end()) {
                    // if its already been explored or banned from exploring,
                    // skip it
                    continue;
                }
                distance[vertex.first] = distance[topPath.first]
                                         + vertex.second;
                nodeBefore[vertex.first] = topPath.first;

                pQueue.push(pair<int, double>{vertex.first,
                                              distance[vertex.first]});
                pushes++;
                if(pushes % 10000 == 0) cout << pushes << endl;
            }
        }
        if (bestPathLength == numeric_limits<double>::min()) {
            // failed to find a path
            Path junk;
            junk.length = bestPathLength;
            return junk;
        }
        return buildPath(nodeBefore, bestPathLength, localSource);
    }
    //==========================TOOLS FOR YEN=================================

    Path slicePathUpToAndEqualTo(Path &path, int upToAndIncluding) {
        /// Slices a path from the start of it up to the end (excluding the end)
        /// Includes the end
        /// Mostly here for access to the graph
        Path newPath{};
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

    void removeEdgeFromGraph(int from, int to, vector<tuple<int, int, double>>
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

    void resetGraphChanges(
            vector<tuple<int, int, double>> &graphChanges) {
        // Essentially graphChanges stores {from, to, length}
        for (auto item : graphChanges) {
            graph[get<0>(item)][get<1>(item)] = get<2>(item);
        }
        graphChanges.clear();
    }

    // =======================YEN VARIATIONS==================================
    vector<Path> basicYen(int numberOfPaths) {
        /// Just does yen from start to finish sequentially
        auto startTime = chrono::high_resolution_clock::now();
        /// Returns a vector of k shortest loopless paths
        vector<Path> confirmedPaths{singleGoalDijkstra()};
        // since set has an O(1) lookup and keeps the highest element at the
        // front, its the easiest to use. An unordered set is another option,
        // but requires us to define a custom hash function
        set<Path, shortestPathComparator> potentialPaths;
        // Instead of fully recovering the graph, we can just 'recover' it by
        // recording what we've done so far to it
        vector<tuple<int, int, double>> graphChanges;

        // since the first one was already found we start from the second index
        for (int i = 0; i < numberOfPaths - 1; i++) {
            // The possible nodes go from the second node up to the second
            // last of the previously found shortest path
            for (int j = 0; j < confirmedPaths[i].route.size() - 1; j++) {
                yenGetNextPath(confirmedPaths, potentialPaths,
                               graphChanges, i, j);
            }

            if (potentialPaths.empty()) {
                break;
            }

            // our potential paths are in a set so they should already be sorted
            confirmedPaths.push_back(*potentialPaths.begin());
            potentialPaths.erase(potentialPaths.begin());
        }

        // For some reason we aren't finding paths in order
        sort(confirmedPaths.begin(), confirmedPaths.end(),
             shortestPathComparator());
        int k = 0;
        cout << "Took " << timeSince(startTime) << "seconds" << endl;
        for (auto &p : confirmedPaths) {
            cout << "K-" << k++ << ":\t" << confirmLength(p) << endl;
        }
        return confirmedPaths;
    }

    vector<Path> randomSelectionYen(int numberOfPaths) {
        /// Essentially we just
        auto startTime = chrono::high_resolution_clock::now();
        vector<Path> confirmedPaths{singleGoalDijkstra()};
        set<Path, shortestPathComparator> potentialPaths;
        vector<tuple<int, int, double>> graphChanges;

        for (int i = 0; i < numberOfPaths - 1; i++) {
            vector<int> unexploredNodes{};
            unexploredNodes.reserve(confirmedPaths[i].route.size());
            for (int j = 1; j < confirmedPaths[i].route.size(); ++j) {
                unexploredNodes.push_back(j);
            }
            shuffle(unexploredNodes.begin(), unexploredNodes.end(),
                    mt19937(random_device()()));

            for (int j = 0; j < (confirmedPaths[i].route.size() - 2) / 10;
                 j++) {
                yenGetNextPath(
                        confirmedPaths, potentialPaths,
                        graphChanges, i,
                        unexploredNodes[j]);
            }

            if (potentialPaths.empty()) {
                cout << "Didn't find anything" << endl;
                break;
            }

            confirmedPaths.push_back(*potentialPaths.begin());
            potentialPaths.erase(potentialPaths.begin());
        }

        sort(confirmedPaths.begin(), confirmedPaths.end(),
             shortestPathComparator());
        cout << "Took " << timeSince(startTime) << "seconds" << endl;
        int k = 0;
        for (auto &p : confirmedPaths) {
            cout << "K-" << k++ << ":\t" << confirmLength(p) << endl;
        }
        return confirmedPaths;
    }

    vector<Path> heatingYen(int numberOfPaths) {
        /// Essentially, the first few confirmed paths are
        /// seen as 'worth' more as the shorter paths are more likely to come
        /// from them. So here, we start by taking 100% of the early paths, and
        /// and gradually decrease the percentage that is taken; in other words,
        /// the opposite of simulated annealing

    }

    vector<Path> SGDYen(int numberOfPaths) {
        /// Since Yen's graph has a pattern on it, we can try to use Stochastic
        /// Gradient Descent to find where to go next. Ideally, we can find
        /// the local minimum (or most local minimums) of the current path
        /// instead of processing through everything. This is especially
        /// useful in enormous paths that were possibly impossible to navigate
        /// before

    }

    void yenGetNextPath(
            vector<Path> &confirmedPaths,
            set<Path, shortestPathComparator> &potentialPaths,
            vector<tuple<int, int, double>> &graphChanges,
            int pathNumber, int startingNodeIndex) {

        int nodeToStartFrom = confirmedPaths[pathNumber].route[startingNodeIndex];
        Path basePath = slicePathUpToAndEqualTo(
                confirmedPaths[pathNumber], startingNodeIndex
        );

        // Disconnect all previously explored connections
        for (auto &path : confirmedPaths) {
            if (checkEqualToRoot(basePath.route,
                                 path.route)) {
                // If they're the same up to that index, remove the
                // next option that path took
                removeEdgeFromGraph(path.route[startingNodeIndex],
                                    path.route[startingNodeIndex + 1],
                                    graphChanges);
            }
        }

        // Generate the path from the new start to the goal
        Path spurPath =
                singleGoalDijkstra(nodeToStartFrom,
                                   unordered_set<int>(
                                           basePath.route.begin(),
                                           basePath.route.end()));

        // If it isn't a garbage path, add it to the potential paths
        // *note: garbage paths have a -double_minimum length
        if (spurPath.length != numeric_limits<double>::min()) {
            basePath += spurPath;
            if (potentialPaths.find(basePath) == potentialPaths.end()) {
                // if basePath not in potential paths, add it
                potentialPaths.insert(basePath);
            }
        }

        // Fix the graph back up
        resetGraphChanges(graphChanges);
    }

    double confirmLength(Path const &p) {
        /// debugging function; if the length that we get from this is
        /// different from what the path has, something is very wrong
        double length = 0;
        for (int i = 0; i < p.route.size() - 1; i++) {
            length += graph[p.route[i]][p.route[i + 1]];
        }
        return length;
    }

private:
    // contains the int index then the distance from the node
    vector<unordered_map<int, double>> graph;
    int source{}, goal{}, K{};
};


int main(int argc, char *argv[]) {
    if (argc <= 1) {
        cout << "Error: an input and an output file is required."
             << endl
             << "Call the program like `program.exe input.txt`"
             << endl;
    }
    Searcher s{argv[1]};

    auto result = s.randomSelectionYen(10);
//    for (int i = 0; i < result.size(); ++i) {
//        cout << i + 1 << ": " << result[i].length << endl;
//    }
    cout << "------------------------------------" << endl;
    result = s.basicYen(10);
    return 0;
}
