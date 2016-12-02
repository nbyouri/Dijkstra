∞ = Inf
g = [0 3 ∞ ∞ ∞ ∞ ∞;
     3 0 2 ∞ 4 ∞ ∞;
     ∞ 2 0 8 1 ∞ ∞;
     ∞ ∞ 8 0 ∞ 5 6;
     ∞ 4 1 ∞ 0 ∞ ∞;
     ∞ ∞ ∞ 5 ∞ 0 ∞;
     ∞ ∞ ∞ 6 ∞ ∞ 0]

function neighbours(v)
  col = find((v .!= ∞) &(v .!= 0))
  col, v[col]
end

function dijkstra(g, i, j)
  nb = size(g, 1)
  dists = fill(∞, nb, nb)
  dists[:,i] = 0
  for k in 1:nb
    """ n vectors of v neighbours """
    n = neighbours(g[k,:])
    """ size of vector is defined by first element size """
    for v in 1:length(n[1])
      println(n[1][v]," " ,n[2][v])
    end
  end
end

println(dijkstra(a, 1, 6))
