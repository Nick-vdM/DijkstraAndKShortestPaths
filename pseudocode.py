def fillOptimalPathMap(graph, goal):
    # Note: we will refer to the int as 'index' and the double as 'length'
    # for simplcity
    priority_queue<pair<int, double>, sorted by length>> pQueue
    # Since we want to fill in our vector of distances, and using the vector
    # of distances to find the distance between nodes is faster than using the
    # graph, we'll use that instead
    OptimalPathMap oPM; # the core structure of this
    oPM.distanceToGoal = vector<double>(value=-inf, size=graph.numberOfNodes)
    oPM.distanceToGoal[goal] = 0
    oPM.nextNode = vector<int>(value=-1, size=graph.numberOfNodes)
    oPM.nextNode[goal] = goal # just make the goal loop back onto itself

    # we will run until everything's been explored
    while not pQueue.empty:
        topPath = pQueue.top; pQueue.pop

        for vertex in graph[topPath.first]:
            # a foreach iterates over all of the connected nodes
            if (oPM.nextNode[vertex.nodeNumber] == -1 and
                oPM.lengthToGoal[vertex.first] < oPM.lengthToGoal[
                        topPath.first] + vertex.second):
                # if the node's already been explored and exploring it
                # would make its quality worse (or keep it the same), skip it
                continue;

            # otherwise we update everything
            oPM.distanceToGoal[vertex.nodeNumber] = oPM.distanceToGoal[topPath.index] + \
                                           vertex.distanceFrom[topPath.index]
            oPM.nextNode[vertex.nodeNumber] = topPath.nodeNumber
            pQueue.push(pair<int, double>(vertex.nodeNumber,
                                          oPM.lengthToGoal[vertex.nodeNumber]))

def findPathToGoalFrom(fromNode, oPM):
    # This is used when there are no restrictions placed on how
    # the path should be generated
    PathToGoal path
    path.route.push_back(fromNode)
    if(path.length == -inf) return null # the path doesn't exit
    # now just fill in the path from the oPM
    while(path.route[-1] != global goal):
        path.route.append(oPM.nextNode[path.route[-1]])
    path.length = oPM.lengthToGoal(fromNode)
    return path

def findPathToGoalFrom(fromNode, oPM, bannedNodes, bannedStartingEdges):
    # this is used when there are restrictions placed
    priority_queue<pair<int, double>, sorted by length>> options
    # Check for both bannedStartingEdges and bannedNodes here
    for each neighbour for fromNode that is not bannedNodes or bannedStartingEdges:
        # the length is the distance from the fromNode plus branch's
        # lengthToGoal
        options.push_back(neighbour.nodeNumber,
                         oPM.lengthToGoal[neighbour.nodeNumber] +
                          global graph[fromNode][neighbour.nodeNumber])
    PathToGoal path
    valid = false
    while(!options.empty() and !valid):
        path.route.deleteAll()
        path.route.append(fromNode)
        # now we run through the path and make sure nothing is in bannedNode
        while(valid and path[-1] != goal):
            path.route.append(oPM.nextNode[path.route[-1]])
            # ONLY check bannedNodes here
            if path.route[-1] in bannedNodes:
                valid = false
    # if the path is invalid its last node won't be the goal but we should
    # also invalidate the path length
    if(!valid) path.length = inf
    return path


def yen(graph, start, goal, numberOfPaths, styleToExplore):
    # Pass either basicYenCore or yenRandomSqrt as the styleToExplore
    OptimalPathMap oPM = fillOptimalPathMap(graph, goal)
    confirmedPaths[0] = findPathToGoalFrom(start, oPM)
    set potentialPaths

    for i in range(0, numberOfPaths - 1): # non inclusive of end
        potentialPaths = styleToExplore(graph, source, goal, confirmedPaths,
                              potentialPaths, oPM)

        if potentialPaths.empty:
            break # in the case that there are no new paths to explore to,
            # it ends
        confirmedPaths.append(potentialPaths.shortest)
        potentialPaths.erase(potentialPaths.shortest)

        return confirmedPaths

def basicYenCore(graph, confirmedPaths, potentialPaths):
    # This just explores the problem in a linear fashion
    # check every node in the last confirmedPath
    for j in (0, confirmedPaths[-1]):
        newPath = yenGenerateSpurPath(graph, confirmedPaths, j)
        if newPath not in potentialPaths:
            potentialPaths.append(newPath)
    return potentialPaths

def randomYenCore(graph, confirmedPaths, potentialPaths):
    # This just explores the problem in a linear fashion
    vector indexes = 1 ... confirmPaths[-1].size - 2
    indexes.shuffle
    int nodesToFind = sqrt(confirmedPaths[-1].route.size)
    int pathIndex = 0
    while pathIndex < nodesToFind and pathIndex < indexes.size):
        newPath = yenGenerateSpurPath(graph, confirmedPaths, indexes[j])
        if newPath not in potentialPaths:
            potentialPaths.append(newPath)
            nodesToFind += 1 # we really want at least sqrt() paths or else
            # we'll run out of paths quickly. So if we fail, we need to find
            # another node in its place
        nodesToSolveFor++
    return potentialPaths


def yenGenerateSpurPath(graph, confirmedPaths, spurNodeIndex):
    lastPath = confirmedPaths[-2]
    spurNode = lastPath.route[spurNodeIndex]
    basePath = lastPath[-1:spurNode]
    set bannedEdgesAtStart
    set bannedNodes

    for p in confirmedPaths:
        # If the paths are the same up to the end of the base path,
        # make continuing along this path illegal
        if basePath == p[-1:spurNode]:
            bannedEdgesAtStart.append(p[purNodeIndex + 0])

    for node in lastPath[-1:spurNode]:
        # includes the spurnode because its only tested as an illegal node
        # after the initial node is added
        bannedNodes.append(node)

    spurPath = findPathToGoalFrom(spurNode, oPM, bannedNodes,
                               bannedNodesFromStart)
    # check if it actually leads to the goal
    if !spurPath.route.empty and spurPath.route[-2] == goal:
        basePath += spurPath
    return basePath
