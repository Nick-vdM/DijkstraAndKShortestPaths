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

struct greatestPath {
    bool operator()(const Path &p1, const Path &p2) {
        return p1.length > p2.length;
    }
};

template<class A, class B>
struct shortestPathComparator {
    bool operator()(const pair<A, B> &p1, const pair<A, B> &p2) {
        return p1.second > p2.second;
    }
};

/*
struct shortestPathComparator {
public:
    explicit shortestPathComparator(vector<double *> &lengths) : lengths{lengths} {}

    bool operator()(int &node1, int &node2) {
        return *lengths[node1] > *lengths[node2];
    }

private:
    vector<double *> lengths;
};
 */

class Searcher {
public:
    Searcher(string const &filePath) {
        int numberOfNodes;
        int numberOfEdges;
        ifstream ifs(filePath);
        ifs >> numberOfNodes >> numberOfEdges;

        graph.reserve(numberOfNodes);
        for (int i = 0; i < numberOfNodes; i++) {
            graph.emplace_back(numberOfNodes, numeric_limits<double>::max());
            connectionLookup.emplace_back();
        }

        // Load in graph information
        for (int i = 0; i < numberOfEdges; i++) {
            int from, to;
            double length;
            ifs >> from >> to >> length;
            graph[from][to] = length;
            connectionLookup[from].push_back(to);
        }
        // Grab the source and goal
        ifs >> source >> goal >> K;
    }

    double dijkstra() {
        // just finds a single path
        auto startTime = chrono::high_resolution_clock::now();
        /*
        vector<double*> lengths;
        lengths.reserve(graph.size());
        for(int i = 0; i < graph.size(); i++){
            lengths.push_back(new double(numeric_limits<double>::max()));
        }
        shortestPathComparator comparator{lengths};
         */

        priority_queue<pair<int, double>,
                vector<pair<int, double>>,
                shortestPathComparator<int, double>> pQueue{};
        pQueue.push(pair<int, double>{source, 0});
        vector<int> nodeBefore(graph.size(), -1);
        double bestPathLength{0};
        while (!pQueue.empty()) {
            auto topPath = pQueue.top();
            pQueue.pop();

            for (auto vertex : connectionLookup[topPath.first]) {
                if (nodeBefore[vertex] != -1) {
                    // if its already been explored, skip it
                    continue;
                }
                nodeBefore[vertex] = topPath.first;

                double newLength = topPath.second +
                                   graph[topPath.first][vertex];

                if (vertex == goal) {
                    bestPathLength = newLength;
                    break;
                }
                pQueue.push(pair<int, double>{vertex, newLength});
            }
        }

        auto endTime = chrono::high_resolution_clock::now();
        double processTime =
                chrono::duration_cast<chrono::microseconds>
                        (endTime - startTime).count();
        cout << "Took " << fixed << processTime / 1000000
             << " seconds to find"
             << endl;
//        printPathTaken(nodeBefore);
        return bestPathLength;
    }

    void printPathTaken(vector<int> &pathMap) {
        int currentNode = goal;
        vector<int> pathTaken{goal};
        while (currentNode != source) {
            currentNode = pathMap[currentNode];
            pathTaken.push_back(currentNode);
        }
        for (int i = pathTaken.size() - 1; i >= 0; i--) {
            cout << pathTaken[i] << " ";
        }
        cout << endl;
        double length = 0;
        for (int i = pathTaken.size() - 1; i > 0; i--) {
            length += graph[pathTaken[i]][pathTaken[i - 1]];
            cout << "From " << pathTaken[i] << " to " << pathTaken[i - 1]
                 << " = " << graph[pathTaken[i]][pathTaken[i - 1]]
                 << "\t Total = " << length << endl;
        }
    }

private:
    vector<vector<double>> graph;
    vector<vector<int>> connectionLookup;
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
    auto result = s.dijkstra();
    cout << result << endl;
//    auto results = s.DijkstraKSP();
//    for(auto r : results){
//        cout << r << " ";
//    }
    return 0;
}
