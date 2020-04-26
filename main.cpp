#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <queue>

using namespace std;

struct Graph {
    map<pair<int, int>, double> edges;
    unordered_map<int, vector<int>> connectionLookup;
};

struct Path {
    unordered_set<int> visitedNodes{};
    int currentNode{0};
    double length{0};

    Path copyAndAppend(int newNode, double newLength) {
        auto newPath = *this;
        newPath.currentNode = newNode;
        newPath.length = newLength;
        newPath.visitedNodes.insert(newNode);
        return newPath;
    }
};

struct greatestPath{
    bool operator()(const Path & p1, const Path & p2){
        return p1.length > p2.length;
    }
};

template<class A, class B>
struct greatestSecondPairValue {
    bool operator()(const pair<A, B> & p1, const pair<A, B> & p2) {
        return p1.second > p2.second;
    }
};

class Searcher {
public:
    Searcher(string const &filePath) {
        int numberOfNodes;
        int numberOfEdges;
        ifstream ifs(filePath);
        ifs >> numberOfNodes >> numberOfEdges;

        // Load in graph information
        for (int i = 0; i < numberOfEdges; i++) {
            int from, to;
            double length;
            ifs >> from >> to >> length;
            graph.edges[pair<int, int>{from, to}] = length;

            auto it = graph.connectionLookup.find(from);

            if (it == graph.connectionLookup.end()) {
                graph.connectionLookup[from] = vector<int>{to};
            } else {
                // append it to the values
                it->second.push_back(to);
            }
        }
        // Grab the source and goal
        ifs >> source >> goal >> K;
    }

    vector<double> plainDijkstra() {
        // just finds a single path
        auto startTime = chrono::high_resolution_clock::now();

        priority_queue<pair<int, double>,
                vector<pair<int, double>>,
                greatestSecondPairValue<int, double>> pQueue{};
        pQueue.push(pair<int, double>{0, 0});
        unordered_set<int> explored;
        vector<double> pathValues;
        int pushed = 0;
        while (!pQueue.empty()) {
            auto topPath = pQueue.top();
            pQueue.pop();
            if (topPath.first == goal) {
                pathValues.push_back(topPath.second);
                if (pathValues.size() == 1) break; // enough paths were found
            } else
                for (auto vertex : graph.connectionLookup[topPath.first]) {
                    if (explored.find(vertex) !=
                        explored.end()) {
                        // if its already been explored, skip it
                        continue;
                    }
                    explored.insert(vertex);

                    double newLength = topPath.second + graph.edges[pair<int,
                            int>{topPath.first, vertex}];
                    pQueue.push(pair<int, double>{vertex, newLength});
                    pushed++;
                }
        }

        auto endTime = chrono::high_resolution_clock::now();
        double processTime =
                chrono::duration_cast<chrono::microseconds>
                        (endTime - startTime).count();
        cout << "Took " << fixed << processTime / 1000000
             << " seconds to find and " << pushed << " push operations"
             << endl;
        return pathValues;
    }

    vector<double> DijkstraKSP() {
        // This approach doesn't work and is gonna get canned
        auto startTime = chrono::high_resolution_clock::now();

        priority_queue<Path, vector<Path>, greatestPath> pQueue;

        Path beginning{};
        beginning.copyAndAppend(0, 0);
        pQueue.push(beginning);
        int pushed = 0;

        vector<double> pathValues;
        while (!pQueue.empty()) {
            auto topPath = pQueue.top();
            pQueue.pop();
            cout << topPath.length << "\t" << topPath.visitedNodes.size()
                 << "\t" << pushed << " " << pQueue.size() << endl;

            if (topPath.currentNode == goal) {
                pathValues.push_back(topPath.length);
                if (pathValues.size() == K) break; // enough paths were found
            } else
                for (auto vertex : graph.connectionLookup[topPath.currentNode]) {
                    if (topPath.visitedNodes.find(vertex) !=
                        topPath.visitedNodes.end()) {
                        // if its already been found, skip it
                        continue;
                    }

                    double newLength = topPath.length + graph.edges[pair<int,
                            int>{topPath.currentNode, vertex}];
                    pQueue.push(topPath.copyAndAppend(vertex, newLength));
                    pushed++;
                }
        }

        auto endTime = chrono::high_resolution_clock::now();
        double processTime =
                chrono::duration_cast<chrono::microseconds>
                        (endTime - startTime).count();
        cout << "Took " << fixed << processTime / 1000000
             << " seconds to find "
             << endl;
        return pathValues;
    };

private:
    Graph graph;
    int source, goal, K;
};


int main(int argc, char *argv[]) {
    if (argc <= 1) {
        cout << "Error: an input and an output file is required."
             << endl
             << "Call the program like `program.exe input.txt`"
             << endl;
    }
    Searcher s{argv[1]};
    auto result = s.plainDijkstra();
    for (auto v : result) {
        cout << v << " ";
    }
    result = s.DijkstraKSP();
    for (auto v : result) {
        cout << v << " ";
    }
    cout << endl;
    return 0;
}
