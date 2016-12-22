using Base.Collections

#
# Graph Edge
#
#
#
type Edge
  # Variables
  from::Int       # source vertex
  to::Int         # destination vertex
  weight::Float64 # edge weight

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
  V::Int                  # number of vertices
  E::Int                  # number of edges
  adj::Array{Array{Edge}} # adjacency matrix

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
      for e in g.adj[v]
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

  function FloydWarshallSP(g::EdgeWeightedGraph)
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

    # Main loop going through the matrix
    for i = 1 : g.V
      for v = 1 : g.V
        # Self loop, don't go
        if (vertexTo[v][i] == 0)
          continue
        end
        # For each neighbour of the vertex
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
# Load graph from a gml file
#
#
#
function loadgraph(file::String)
  # Open the file on disk
  f = open(file)
  # Load the contents into a string
  lines = readlines(f)

  # Load edges in an array
  edges = Edge[]
  for l in 1 : length(lines)
    # Only deal with edge tags
    if (contains(lines[l], "edge"))
      edgeLine = l + 2 # Skip edge and [
      source = ""
      target = ""
      weight = ""
      # Read the edge block until ]
      while (!contains(lines[edgeLine], "]"))
        curLine = split(lines[edgeLine])
        # Avoid incorrect gml
        if (length(curLine) > 2)
          println("Unsupported graph")
          return
        end
        # Get the tag values we need
        if (curLine[1] == "source")
          source = curLine[2]
        elseif (curLine[1] == "target")
          target = curLine[2]
        elseif (curLine[1] == "value")
          weight = curLine[2]
        end # we don't support other tags
        edgeLine += 1
      end
      # If we don't have a source, target or weight for an edge we don't support
      if ((length(source) == 0) ||
        (length(target) == 0) ||
        (length(weight) == 0))
        println("Error: Unsupported graph")
        return
      else
        # Add a new edge on the array
        push!(edges, Edge(parse(Int, source),
          parse(Int, target), parse(Float64, weight)))
      end
    end
  end

  # Load edges in a map with the source as key to ensure unicity to build
  # the adjacency matrix
  adj = Dict()
  for e in edges
    # Initialise the edge array
    if (!haskey(adj, e.from))
      adj[e.from] = Edge[]
    end
    if (!haskey(adj, e.to))
      adj[e.to] = Edge[]
    end
    # Add the edges for the source and target vertex
    push!(adj[e.from], e)
    push!(adj[e.to], Edge(e.to, e.from, e.weight))
  end

  # Convert the dict to the adjacency 2D array
  fadj = Array{Edge}[]
  # Create a sorted index vector
  adjIndices = collect(keys(adj))
  sort!(adjIndices)
  for v in adjIndices
    push!(fadj, adj[v])
  end

  # Close the file
  close(f)

  # Return the graph
  return EdgeWeightedGraph(length(adjIndices), length(edges), fadj)
end

function coût(g::EdgeWeightedGraph, i::Int, j::Int)
  d = DijkstraSP(g, i)
  fw = FloydWarshallSP(g)

  println("Dijkstra")
  println(d.dist(j))
  println(d.pathTo(j))
  println("Floyd-Warshall")
  println(fw.dist(i, j))
  println(fw.path(i, j))
end

#
# Main
#
#
#
# Load graph from a file where vertices are indexed from 1 to N
g = loadgraph("/Users/youri/Downloads/graph.gml")
# Or from an adjacency matrix, adj being Array{Array{Edge}}
# g = EdgeWeightedGraph(V, E, adj)

# Get the shortest path from 1 to 5 on g
coût(g, 1, 5)
