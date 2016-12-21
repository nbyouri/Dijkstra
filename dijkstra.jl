using Base.Collections

# TODO: load gml

""" Graph Edge """
type Edge
  from::Int
  to::Int
  weight::Int

  toString::Function

  function Edge(from::Int, to::Int, weight::Int)
    this = new()
    this.from = from
    this.to = to
    this.weight = weight

    # Print Edge fields
    this.toString = function ()
      println(this.from, " <-> ", this.to, " = ",  this.weight)
    end

    return this
  end
end

""" Edge Weighted undirected graph """
type EdgeWeightedGraph
  V::Int
  E::Int
  adj::Array{Array{Edge}}

  toString::Function

  function EdgeWeightedGraph(V::Int, E::Int, adj::Array{Array{Edge}})
    this = new()
    this.V = V
    this.E = E
    this.adj = adj

    # Print adjacency matrix for the graph
    this.toString = function ()
      for i = 1 : V
        for j = 1 : length(adj[i])
          adj[i][j].toString()
        end
      end
    end

    return this
  end
end

""" Shortest path Dijkstra algorithm """
type DijkstraSP
  pq::PriorityQueue
  distTo::Array{Int}
  edgeTo::Array{Edge}

  relax::Function
  pathTo::Function
  dist::Function
  hasPathTo::Function

  function DijkstraSP(g::EdgeWeightedGraph, s::Int)
    this = new()
    pq = PriorityQueue()
    # double ?
    distTo = Array{Int}(g.V)
    edgeTo = Array{Edge}(g.V)

    # Core of the algorithm, walk the graph
    this.relax = function (v::Int)
      for i = 1 : length(g.adj[v])
        e = g.adj[v][i]
        w = e.to
        alt = distTo[v] + e.weight
        if (distTo[w] > alt)
          distTo[w] = alt
          edgeTo[w] = e
          # XXX change key? -- check complexity of dequeue and enqueue
          try
            dequeue!(pq, w)
          catch e
          end
          enqueue!(pq, w, distTo[w])
        end
      end
    end

    # Returns true if there is a path to destination
    this.hasPathTo = function (to::Int)
      return pathTo[to] < 100 #  XXX
    end

    # Path to a vertex
    # returns a vector of nodes to the destination
    this.pathTo = function (to::Int)
      if !hasPathTo(to)
        return
      end

      #for i = 1 : edgeTo
      #end
    end

    this.dist = function (to::Int)
      return distTo[to]
    end

    # Initialise paths to infinity
    for v = 1 : length(distTo)
      distTo[v] = 100 # XXX
    end
    # Initialise source to 0
    distTo[s] = 0

    # Add the source to the queue
    enqueue!(pq, s, 0)

    # Walk through the graph to update distTo
    while (!isempty(pq))
      this.relax(dequeue!(pq))
    end

    return this
  end
end

"""
test undirected weighted graph
------------------------------

A          5            B
 +----------------------XXXXXXXXX   1
 |                      |       XXXXXXX
 |                      |             XXXX
 |                      |                XXXXX
 |                      |                   XXX E
1|                      |1                XXX
 |                      |               XXX
 |                      |          XXXXX
 |                      |       XXXX 3
 |                      |   XXXXX
 +---------------------XXXXXX
C           0            D

distTo[E] should be 3 and the path [A, C, D, B, E]
"""
ab = Edge(1, 2, 5)
ac = Edge(1, 3, 1)

ba = Edge(2, 1, 5)
bd = Edge(2, 4, 1)
be = Edge(2, 5, 1)

ca = Edge(3, 1, 1)
cd = Edge(3, 4, 0)

dc = Edge(4, 3, 0)
db = Edge(4, 2, 1)
de = Edge(4, 5, 3)

eb = Edge(5, 2, 1)
ed = Edge(5, 4, 3)

# adjacency lists
adj = Array{Edge}[]
adjA = Edge[]
push!(adjA, ab) # A <-> B
push!(adjA, ac) # A <-> C

adjB = Edge[]
push!(adjB, ba) # B <-> A
push!(adjB, bd) # B <-> D
push!(adjB, be) # B <-> E

adjC = Edge[]
push!(adjC, ca) # C <-> A
push!(adjC, cd) # C <-> D

adjD = Edge[]
push!(adjD, dc) # D <-> C
push!(adjD, db) # D <-> B
push!(adjD, de) # D <-> E

adjE = Edge[]
push!(adjE, eb) # E <-> B
push!(adjE, ed) # E <-> D

push!(adj, adjA)
push!(adj, adjB)
push!(adj, adjC)
push!(adj, adjD)
push!(adj, adjE)

# graph with 5 vertices, 6 edges
g = EdgeWeightedGraph(5, 6, adj)
d = DijkstraSP(g, 1)
println(d.dist(5))
# print adjacency list
