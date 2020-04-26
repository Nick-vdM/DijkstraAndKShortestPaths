
# Computing Algorithms - 2801ICT - Searching Algorithms #

I'm going to attempt a speedrun of this assignment.
 * First hour : Plan out algorithm & report
    * 30 mins: read/understand problem
    * 30mins : Look at Djikstra's and stuff
 * Second & Third hour : Write out code for things
 * Fourth & Fifth hour : Debug
 * 5+ hours : Report
 
 ## Hour 1 - Plan out
The restrictions in the assignment are incredibly confusing and lax. However,
it looks pretty obvious that Saiful wants us to extend Djikstra's algorithm
by the looks of everything. 
  
### Djikstra's Algorithm Summary
The wikipedia page's pseudocode is nonsensical so I'll start with that.
 
```{python, tidy=FALSE}
Graph(set of vertices, set of edges)
findWeight(from, to)

Path = empty
foundPaths = 0
Insert start into heap
While B not empty and count < limit:
    lowestCost = heap.pop(); heap.removeTop()
    count++
    If foundPaths <= K:
        for each vertex adjacent to current node:
            Insert new path into heap
```

So it looks like you basically just use Djikstra's algorithm except you don't
stop after you find the first path - instead you keep exploring states until
you find K paths. That's a lot less interesting than I thought it'd be
... This just seems like UCS that doesn't stop to be honest. Maybe I'm not
 understanding what K shortest paths is - I'll have a look.
 
So really. Could I just reuse OOP's code from 1810ICT and adjust it? That
seems way way WAY too easy. Is there a more innovative way to do this?

### Yen's Algorithm
Let's look at Yen's. Essentially it consists of two paths: finding the
shortest path, then exploring short deviations of that path.

```{python, tidy=FALSE}
Paths[0] = Djikstra()

for k in range(0, K):
    # From the last checked node up to the back
    for i in range(0, size(Shortest Path[k - 1]) - 2):
        # The previous k-shortest path k - 1
        spur node = path[k - 1].node(i)
        Root path = path[k - 1].nodes(0, i)
        
        for each path in paths:
            if rootPath == path.nodes(0, i):
                // ban considering similar paths
                remove p.edge(i, i + 1) from Graph
        
        for each node in rootPath except spur node:
            remove node from Graph
        
        spurPath = Djikstra(Graph, spurNode, goal)
        
```

So really, Yen's algorithm is also an extension of Djikstra's - of sorts. 
Since it largely still requires the same things... Maybe the best approach is
just to implement both of them?

## Hour 2 - Code
I found a very sexy Github repo that basically has both Djikstra's and Yen's
algorithm on it [here](https://github.com/yan-qi/k-shortest-paths-cpp-version).
Basically I'll just use its structural design and merge it into a simpler
, single file.

* 15 mins mark - Complete input
* 40 mins - Completed elementary Dijktras (pathing was broken)
* 1 hr - Completed Djikstra's, but I think the input is broken now.

