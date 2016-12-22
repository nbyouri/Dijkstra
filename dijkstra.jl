using Base.Collections

#
# Graph Edge
#
#
#
type Edge
  # Variables
  from::Int
  to::Int
  weight::Float64

  # Constructor
  function Edge(from::Int, to::Int, weight::Float64)
    this = new()
    this.from = from
    this.to = to
    this.weight = weight

    return this
  end
end

#
# Edge Weighted undirected graph
#
#
#
type EdgeWeightedGraph
  # Variables
  V::Int
  E::Int
  adj::Array{Array{Edge}}

  # Constructor
  function EdgeWeightedGraph(V::Int, E::Int, adj::Array{Array{Edge}})
    this = new()
    this.V = V
    this.E = E
    this.adj = adj

    return this
  end
end

#
# Shortest path Dijkstra algorithm
#
#
#
type DijkstraSP
  # Variables
  pq::PriorityQueue      # minimum oriented priority queue
  distTo::Array{Float64} # array of weights to vertices
  vertexTo::Array{Edge}  # array of edges to vertices

  # Methods
  relax::Function
  pathTo::Function
  dist::Function
  hasPathTo::Function

  # Constructor
  function DijkstraSP(g::EdgeWeightedGraph, s::Int)
    this = new()

    # Core of the algorithm, walk the graph
    this.relax = function (v::Int)
      # For each vertex neighbour
      for i = 1 : length(g.adj[v])
        e = g.adj[v][i]
        w = e.to
        alt = distTo[v] + e.weight
        # If the weight path from source + the edge's weight is smaller than
        # the distance recorded to w, update it
        if (distTo[w] > alt)
          # Update distance to vertex
          distTo[w] = alt
          # Update path to vertex
          vertexTo[w] = e.from
          # Change vertex priority
          pq[w] = distTo[w]
        end
      end
    end

    # Returns true if there is a path to destination
    this.hasPathTo = function (to::Int)
      return distTo[to] < Inf
    end

    # Path to a vertex
    # returns a vector of nodes to the destination
    this.pathTo = function (to::Int)
      # Return if there is no path
      if (!this.hasPathTo(to))
        return
      end

      # Build array of vertex indices going back from destination to source
      path = Int[]
      # Get the vertex index closest to destination
      from = vertexTo[to]
      while (from != 0)
        # Push the vertex index on the array
        push!(path, from)
        # Get the parent vertex index
        from = vertexTo[from]
      end
      # Reverse to obtain natural ordering on the array
      reverse!(path)
      # Add the destination vertex index
      push!(path, to)

      return path
    end

    # Return distance to destination
    this.dist = function (to::Int)
      return distTo[to]
    end

    # Initialise fields
    pq = PriorityQueue()
    distTo = Array{Float64}(g.V)
    vertexTo = Array{Int}(g.V)

    # Initialise paths to infinity
    # Also initialise all vertex path to 0
    for v = 1 : length(distTo)
      distTo[v] = Inf
      vertexTo[v] = 0
    end

    # Initialise the minimum oriented priority queue with vertex indices
    for v = 1 : length(g.adj)
      enqueue!(pq, v, Inf)
    end

    # Initialise source to 0
    distTo[s] = 0

    # Set the priority of the source vertex to 0
    pq[s] = 0

    # Walk through the graph to update distTo thus finding the shortest path
    # for all other vertices
    while (!isempty(pq))
      this.relax(dequeue!(pq))
    end

    return this
  end
end

#
# Shortest path Floyd-Warshall algorithm
#
#
#
type FloydWarshallSP
  # Variables
  distTo::Array{Array{Float64}}
  vertexTo::Array{Array{Edge}}

  # Methods
  path::Function
  dist::Function
  hasPath::Function

  function FloydWarshallSP(g::EdgeWeightedGraph, s::Int)
    this = new()

    # Returns true if there is a path to destination
    this.hasPath = function (from::Int, to::Int)
      return distTo[from][to] < Inf
    end

    # Path as a vertex vector from source to destination
    # returns a vector of nodes indices to the destination
    this.path = function (from::Int, to::Int)
      # Return if there is no path
      if (!this.hasPath(from, to))
        return
      end

      # Build vector going back from destination to source
      path = Int[]
      # Get the parent vector closest to destination
      parent = vertexTo[from][to]
      while (parent != 0)
        # Add the vertex to the vector
        push!(path, parent)
        # Fetch parent vertex
        parent = vertexTo[from][parent]
      end
      # Reverse the vector order to get natural ordering
      reverse!(path)
      # Add destination
      push!(path, to)

      return path
    end

    # Return distance to destination
    this.dist = function (from::Int, to::Int)
      return distTo[from][to]
    end

    # Initialise fields
    distTo = Array{Array{Float64}}(g.V)
    vertexTo = Array{Array{Int}}(g.V)

    # Initialise paths to infinity
    # Also initialise all vertex path to 0
    for v = 1 : length(distTo)
      distTo[v] = Array{Float64}(g.V)
      vertexTo[v] = Array{Int}(g.V)
      for w = 1 : length(distTo[v])
        distTo[v][w] = Inf
        vertexTo[v][w] = 0
      end
    end

    # Initialise distances and paths
    for v = 1 : g.V
      for w = 1 : length(g.adj[v])
        e = g.adj[v][w]
        distTo[e.from][e.to] = e.weight
        vertexTo[e.from][e.to] = e.from
      end
      # Handle self loops
      if (distTo[v][v] >= 0.0)
        distTo[v][v] = 0.0
        vertexTo[v][v] = 0
      end
    end

    # FW updates
    for i = 1 : g.V
      for v = 1 : g.V
        # Optimize
        if (vertexTo[v][i] == 0)
          continue
        end
        for w = 1 : g.V
          # alt is the weight path to the neighbour vertex added to the weight
          # from source to that vertex
          alt = distTo[v][i] + distTo[w][i]
          # Update the weights and parent vertex if alt is smaller, meaning it
          # will cost less to pass through that edge
          if (distTo[v][w] > alt)
            distTo[v][w] = alt
            vertexTo[v][w] = vertexTo[i][w]
          end
        end
      end
    end

    return this
  end
end


#
#test with an acyclic undirected weighted graph
#----------------------------------------------

#A          5            B
# +----------------------XXXXXXXXX   1
# |                      |       XXXXXXX
# |                      |             XXXX
# |                      |                XXXXX
# |                      |                   XXX E
#1|                      |1                XXX
# |                      |               XXX
# |                      |          XXXXX
# |                      |       XXXX 3
# |                      |   XXXXX
# +---------------------XXXXXX
#C           0            D
#
# From A, distTo[E] should be 3 and the path [A, C, D, B, E] and thus the output
# vector should be [1,3,4,2,5]
#
ab = Edge(1, 2, 5.0)
ac = Edge(1, 3, 1.0)

ba = Edge(2, 1, 5.0)
bd = Edge(2, 4, 1.0)
be = Edge(2, 5, 1.0)

ca = Edge(3, 1, 1.0)
cd = Edge(3, 4, 0.0)

dc = Edge(4, 3, 0.0)
db = Edge(4, 2, 1.0)
de = Edge(4, 5, 3.0)

eb = Edge(5, 2, 1.0)
ed = Edge(5, 4, 3.0)

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
f = FloydWarshallSP(g, 1)
println(f.dist(1, 5))
println(f.path(1, 5))
