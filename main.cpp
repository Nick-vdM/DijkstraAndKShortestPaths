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

using namespace std;

struct Path {
    vector<int> route;
    long double length{0};

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
    bool operator()(const pair<int, long double> &p1,
                    const pair<int, long double> &p2) {
        return p1.second > p2.second;
    }
};

struct shortestPathComparator {
    bool operator()(Path const &p1, Path const &p2) {
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
            if (graph[from].find(to) == graph[from].end()) {
                graph[from][to] = length;
            }
        }
        // Grab the source and goal
        ifs >> source >> goal >> K;
    }

    // =======================GENERAL TOOLS=================================
    array<int, 3> setToLongPath(long double maxSearchTime) {
        /// This function sets the source and goal to be as far apart as
        /// possible in the given amount of time
        // TODO: fix this? lol
        auto startClock = chrono::high_resolution_clock::now();
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> uni(0, 50);
        long double timeSinceStart = 0;
        array<int, 3> startGoalLength{-1, -1, numeric_limits<int>::min()};
        // note: the number of nodes is more impactful
        // to the k-1 search statespace than the actual length of the path
        while (timeSinceStart < maxSearchTime) {
            array<int, 3> newStartGoalLength{-1, -1, -1};
            for (int i = 0; i < 2; i++) newStartGoalLength[i] = uni(rng);
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
    Path
    buildPath(vector<int> &nodeBefore, long double &length, int localSource) {
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
        return singleGoalDijkstra(source, goal, unordered_set<int>{},
                                  unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource) {
        return singleGoalDijkstra(localSource, goal, unordered_set<int>{},
                                  unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource, int localGoal) {
        return singleGoalDijkstra(localSource, localGoal,
                                  unordered_set<int>{}, unordered_set<int>{});
    }

    Path singleGoalDijkstra(int localSource, unordered_set<int> bannedNodes,
                            unordered_set<int> bannedEdgesFromStart) {
        return singleGoalDijkstra(localSource, goal, bannedNodes,
                                  bannedEdgesFromStart);
    }

    Path singleGoalDijkstra(int localSource, int localGoal, unordered_set<int>
    bannedNodes, unordered_set<int> bannedEdgesFromStart) {
        /// returns the path that was taken
        /// uses the class goal instead as unique goals are never used
        /// returns a large negative empty path if it fails to find it
        /// (long double min value)

        priority_queue<pair<int, long double>,
                vector<pair<int, long double>>,
                shortestPathLengthComparator> pQueue{};
        pQueue.push(pair<int, long double>{localSource, 0});
        // we could use our graph structure here, however, this has a faster
        // lookup (but still O(1) and uses O(n) space anyways. So we get a tiny
        // speed boost in exchange for a tiny amount of extra space, and we
        // are valuing the speed higher than the space
        vector<long double> distance(graph.size(),
                                     numeric_limits<long double>::max());
        distance[localSource] = 0;
        vector<int> nodeBefore(graph.size(), -1);
        long double bestPathLength{numeric_limits<long double>::min()};
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
                    //                        distance[vertex.first] <
                    //                    distance[topPath.first] + vertex.second ||
                    bannedNodes.find(vertex.first) != bannedNodes.end()) {
                    // if its already been explored or banned from exploring,
                    // skip it
                    continue;
                }
                // if its our first topPath, check whether the vertex
                // is a banned edge
                if (topPath.first == localSource &&
                    bannedEdgesFromStart.find(vertex.first)
                    != bannedEdgesFromStart.end()) {
                    continue;
                }
                distance[vertex.first] = distance[topPath.first]
                                         + vertex.second;
                nodeBefore[vertex.first] = topPath.first;

                pQueue.push(pair<int, long double>{vertex.first,
                                                   distance[vertex.first]});
            }
        }
        if (bestPathLength == numeric_limits<long double>::min()) {
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

    void resetGraphChanges(
            vector<tuple<int, int, long double>> &graphChanges) {
        // Essentially graphChanges stores {from, to, length}
        for (auto item : graphChanges) {
            graph[get<0>(item)][get<1>(item)] = get<2>(item);
        }
        graphChanges.clear();
    }

    // =======================YEN VARIATIONS==================================

    vector<Path> yen(int numberOfPaths, std::function<vector<Path> &
            (vector<Path> &,
             set<Path, shortestPathComparator> &, int)> fun) {
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
        vector<tuple<int, int, long double>> graphChanges;

        // since the first one was already found we start from the second index
        for (int i = 0; i < numberOfPaths - 1; i++) {
            // The possible nodes go from the second node up to the second
            // last of the previously found shortest path
            fun(confirmedPaths, potentialPaths, i);

            if (potentialPaths.empty()) {
                cout << "ran out of paths" << endl;
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
            cout << "K-" << ++k << ":\t" << setprecision(7) << confirmLength
                    (p) << endl;
            if (k >= 10) break;
        }
        return confirmedPaths;
    }

    vector<Path> &getPotentialPaths_basic(
            vector<Path> &confirmedPaths,
            set<Path, shortestPathComparator> &potentialPaths, int i) {
        /// This just runs the normal Yen's core
        /// which iterates through all the nodes in a given path and
        /// defines confirmedPaths
        for (int j = 0; j < confirmedPaths[i].route.size() - 1; j++) {
            Path newPath = yenGenerateSpurPath(confirmedPaths, j);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
            }
        }
        return confirmedPaths;
    }

    vector<Path> &getPotentialPaths_random(
            vector<Path> &confirmedPaths,
            set<Path, shortestPathComparator> &potentialPaths, int i, double
            const fractionToTest) {
        /// Instead of Yen's core, this function iterates through a
        /// fraction of the route randomly
        // TODO: actually make this function
        for (int j = 0; j < confirmedPaths[i].route.size() - 1; j++) {
            Path newPath = yenGenerateSpurPath(confirmedPaths, j);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
            }
        }
        return confirmedPaths;
    }

    vector<Path> &getPotentialPaths_SGD(
            vector<Path> &confirmedPaths,
            set<Path, shortestPathComparator> &potentialPaths, int i) {
        /// Here we use stochastic gradient descent to attempt to find
        /// a single local minimum
        // TODO: actually make this function
        for (int j = 0; j < confirmedPaths[i].route.size() - 1; j++) {
            Path newPath = yenGenerateSpurPath(confirmedPaths, j);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
            }
        }
        return confirmedPaths;
    }

    vector<Path> &getPotentialPaths_SGDComplete(
            vector<Path> &confirmedPaths,
            set<Path, shortestPathComparator> &potentialPaths, int i) {
        /// Here we use stochastic gradient descent to attempt to find
        /// all of the local minimums. Essentially, we find the second
        /// derivative of the function and if its concave up we speed up the
        /// learning rate, and if its concave down we slow down the learning
        /// rate. Ideally it only explores the concave downs of the function

        // TODO: actually make this function
        for (int j = 0; j < confirmedPaths[i].route.size() - 1; j++) {
            Path newPath = yenGenerateSpurPath(confirmedPaths, j);
            if (potentialPaths.find(newPath) == potentialPaths.end() &&
                newPath.route[newPath.route.size() - 1] == goal) {
                potentialPaths.insert(newPath);
            }
        }
        return confirmedPaths;
    }

    Path yenGenerateSpurPath(vector<Path> &confirmedPaths,
                             int spurNodeIndex) {
        Path &lastPath = confirmedPaths[confirmedPaths.size() - 1];
        int spurNode = lastPath.route[spurNodeIndex];
        Path basePath = slicePathUpToAndEqualTo(lastPath, spurNodeIndex);
        unordered_set<int> bannedStartingEdges;
        unordered_set<int> bannedNodes;

        // Filter out future nodes that were explored
        for (auto &p : confirmedPaths) {
            if (checkEqualToRoot(basePath.route, p.route))
                bannedStartingEdges.insert(p.route[spurNodeIndex + 1]);
        }

        // Filter out the base path
        for (int i : basePath.route) {
            bannedNodes.insert(i);
        }

        // Find the path and combine it with the base
        Path spurPath = singleGoalDijkstra(
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

    long double confirmLength(Path const &p) {
        /// debugging function; if the length that we get from this is
        /// different from what the path has, something is very wrong
        long double length = 0;
        for (int i = 0; i < p.route.size() - 1; i++) {
            length += graph[p.route[i]][p.route[i + 1]];
        }
        return length;
    }

private:
    // contains the int index then the distance from the node
    vector<unordered_map<int, long double>> graph;
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

    using namespace std::placeholders;
    auto fp = bind(&Searcher::getPotentialPaths_basic, s, _1, _2, _3);
    auto result = s.yen(10, fp);

    return 0;
}
