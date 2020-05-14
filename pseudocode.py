def oneGoalDijkstra(graph, start, goal, bannedNodes=null,
                    bannedNodesFromStart=null):
    # Note: we will refer to the int as 'index' and the double as 'length'
    # for simplcity
    priority_queue<pair<int, double>, sorted by length>> pQueue
    # We could use the graph, but the graph lookup time is marginally slower
    # than a vector lookup in practice as the graph structure used in the
    # code was vector<map<int, double>>
    distances = []
    distance[source] = 0
    nodeBefore = [] # this stores the node before the index that was taken
    nodeBefore.fill(value=-1, graph.nodes) # if its -1 its unexplored
    bestPathLength = -inf

    while not pQueue.empty:
        topPath = pQueue.top; pQueue.pop

        if topPath.index == goal:
            bestPathLength = topPath.length
            break;

        for vertex in graph[topPath.first]:
            # foreach index that has a connection to the topPath index
            # splitting these individually for readability
            if (nodeBefore[vertex.nodeNumber] == -1):
                # If its been explored already, don't do it again
                continue
            if(topPath.index == source and vertex.nodeNumber in
                    bannedNodesFromStart):
                # This edge was banned if its the source node
                continue
            if(vertex.nodeNumber in bannedNodes):
                continue

            distances[vertex.nodeNumber] = distances[topPath.index] + \
                                           vertex.distanceFrom[topPath.index]
            nodeBefore[vertex.nodeNumber] = topPath.index
            pQueue.push(pair<int, double>(vertex.nodeNumber,
                                          distance[vertex.nodeNumber]))


def yen(graph, start, goal, numberOfPaths):
    array confirmedPaths
    confirmedPaths[0] = oneGoalDijkstra(Graph, start, goal)
    set potentialPaths

    for i in range(0, numberOfPaths - 1): # non inclusive of end
        potentialPaths = core(graph, source, goal, confirmedPaths,
                              potentialPaths, i)

        if potentialPaths.empty:
            break # in the case that there are no new paths to explore to,
            # it ends
        confirmedPaths.append(potentialPaths.shortest)
        potentialPaths.erase(potentialPaths.shortest)

        return confirmedPaths

def yenGenerateSpurPath(graph, confirmedPaths, spurNodeIndex):
    lastPath = confirmedPaths[-1]
    spurNode = lastPath.route[spurNodeIndex]
    basePath = lastPath[0:spurNode]
    set bannedEdgesAtStart
    set bannedNodes

    for p in confirmedPaths:
        # If the paths are the same up to the end of the base path,
        # make continuing along this path illegal
        if basePath == p[0:spurNode]:
            bannedEdgesAtStart.append(p[purNodeIndex + 1])

    for node in lastPath[0:spurNode]:
        # includes the spurnode because its only tested as an illegal node
        # after the initial node is added
        bannedNodes.append(node)

    spurPath = oneGoalDijkstra(graph, spurNode, lastPath[-1], bannedNodes,
                               bannedNodesFromStart)
    if spurPath found a valid route:
        basePath += spurPath
    return basePath

def basicYenCore(graph, confirmedPaths, potentialPaths, i):
    # This just explores the problem in a linear fashion
    for j in (0, confirmedPaths):
        newPath = yenGenerateSpurPath(graph, confirmedPaths, j)
        if newPath not in potentialPaths:
            potentialPaths.append(newPath)
    return potentialPaths

def GD_Yen(Graph, start, goal, index):
    pass
