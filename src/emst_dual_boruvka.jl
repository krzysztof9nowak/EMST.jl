
using DataStructures
using Distances
using Statistics

import Base.isequal
import Base.hash
import Base.copy

struct Edge
    a::Int64
    b::Int64
    Edge(pa::Int64, pb::Int64) = (pa < pb) ? new(pa, pb) : new(pb, pa)
end
function isequal(a::Edge, b::Edge)
    return a.a == b.a && a.b == b.b
end
function hash(a::Edge)
    return xor(hash(a.a), hash(a.b))
end
function to_matrix(ee::Vector{Edge})
    m::Array{Int64,2} = Array{Int64}(undef, length(ee), 2)
    for zi = 1:length(ee)
        m[zi,1] = ee[zi].a
        m[zi,2] = ee[zi].b
    end
    return m
end

"""
  DTBNode

implements node of a KD tree
"""
mutable struct DTBNode
    # data / indidces of this kd tree node
    data::Array{Float64,2}
    subset::Vector{Int64}

    # kd tree node bounds
    lb::Vector{Float64}
    ub::Vector{Float64}

    # max_{v in this node} ( length of candidate edge of v)
    dQ::Float64

    left::DTBNode
    right::DTBNode

    # list of forest roots
    roots::Vector{Int64}

    function DTBNode(data::Array{Float64,2}, subset::Vector{Int64}, box_lb::Vector{Float64}, box_ub::Vector{Float64}, roots::Vector{Int64})
        n = new(data, subset, box_lb, box_ub)
        n.left = n
        n.right = n
        n.roots = roots
        n.dQ = 0
        return n
    end
end

function is_leaf(n::DTBNode)
    return n.left == n && n.right == n
end

"""
  kdtree(xx::Array{Float64,2})

initializes root node of KD tree (but does not compute the splits)
To build the tree use:
  kdt_root = kdtree( data )
  kdtree_split!( kdt_root , 10 ) # 10 is the node size where we stop splitting
"""
function kdtree(xx::Array{Float64,2})
    root = DTBNode(xx, 
                   collect(Int64(1):Int64(size(xx, 2))), 
                   minimum(xx, dims=2)[:, 1], 
                   maximum(xx, dims=2)[:, 1], 
                   collect(Int64(1):Int64(size(xx, 2))))
    return root
end


"""
  kdtree_split!(node::DTBNode, nmin::Int64)

computes splits of KD tree. nmin indicates the node size where we stop
splitting nodes.
"""
function kdtree_split!(node::DTBNode, nmin::Int64)
    # ensure that parent node is in the tree.
    # t[node.id] = node

    if ( size(node.data, 2) <= nmin )
        return node
    end

    mind = minimum(node.data, dims=2)
    maxd = maximum(node.data, dims=2)
    ds   = findmax(maxd - mind)
    ds   = ds[2][1]
    vs   = median(node.data[ds,:])
    bx   = node.data[ds,:] .<= vs

    range_a = node.subset[bx]
    range_b = node.subset[.!bx]

    data_a = node.data[:,bx]
    data_b = node.data[:,.!bx]

    box_lb_a = copy(node.lb)
    box_ub_a = copy(node.ub)
    box_ub_a[ds] = vs
    box_lb_b = copy(node.lb)
    box_lb_b[ds] = vs
    box_ub_b     = copy(node.ub)

    node.left = DTBNode(data_a, range_a, box_lb_a, box_ub_a, range_a)
    node.right = DTBNode(data_b, range_b, box_lb_b, box_ub_b, range_b)

    kdtree_split!(node.left, nmin)
    kdtree_split!(node.right, nmin)
end

"""
  compute_emst(data::Array{Float64,2};nmin::Int64=64)

Computes EMST for the given data (where columns are samples).
nmin is the max number of elements in kd-tree node.
"""
function compute_emst(data::Array{Float64,2};nmin::Int=64, plot::Bool = false)
    nmin64 = Int64(nmin)
    root = kdtree(data)
    kdtree_split!(root, nmin64)
    edges = dtb(root, IntDisjointSets(size(data, 2)); plot=plot)

    length :: Float64 = 0
    for edge in edges
        A = data[edge.a]
        B = data[edge.b]
        length += sqrt(sum((A .- B).^2))
    end

    @show length

    return to_matrix(collect(edges))
end


mutable struct CandidateEdge
    distance :: Float64
    a :: Int64
    b :: Int64
end


"""
  dtb(q::DTBNode,e::IntDisjointSets)

implements the dual-tree Boruvka algorithm which computes the EMST for the
given KD-tree.
"""
function dtb(Q::DTBNode, forest::IntDisjointSets; plot::Bool = false)
    edges = Set{Edge}()

    while forest.ngroups > 1
        ngroups = forest.ngroups;
        println("--> ngroups: $ngroups")

        C :: Vector{CandidateEdge} = [CandidateEdge(Inf, 0, 0) for i=1:length(Q.subset)]

        @time find_component_neighbors(Q, Q, forest, C)

        # and now add the edges..
        for nn::CandidateEdge in C
            if nn.a != 0 && !in_same_set(forest, nn.a, nn.b)
                union!(forest, nn.a, nn.b) ## add check if union occured
                push!(edges, Edge(nn.a, nn.b))
            end
        end

        update_node_roots(Q, forest)

        if plot
            p = plot_debug(Q.data, forest, edges)
            gui(p)
            sleep(1)
        end
    end
    return edges
end

function update_node_roots(Q::DTBNode, forest)
    if is_leaf(Q)
       Q.roots = unique([find_root!(forest, x) for x in Q.subset])
    else
        update_node_roots(Q.left, forest)
        update_node_roots(Q.right, forest)
        Q.roots = unique([Q.left.roots; Q.right.roots])
    end
end

Base.copy(s::IntDisjointSets) = IntDisjointSets(copy(s.parents), copy(s.ranks), s.ngroups)

function find_component_neighbors(Q::DTBNode, R::DTBNode, forest::IntDisjointSets, C::Vector{CandidateEdge})
    if length(Q.roots) == 1 && length(R.roots) == 1 && R.roots[1] == Q.roots[1]
        return
    end

    if distance(Q, R) > Q.dQ
        return
    end

    if is_leaf(Q) && is_leaf(R)
        Q.dQ = 0.0

        dims = size(Q.data)[1]
        for iq=1:length(Q.subset)
            qq :: Int64 = Q.subset[iq]
            cq :: Int64 = find_root!(forest, qq) # tree of q
            nn ::CandidateEdge = C[cq]
            for ir=1:length(R.subset)
                rr :: Int64 = R.subset[ir]

                if in_same_set(forest, qq, rr)
                    continue
                end

                dist_qr :: Float64 = 0
                @simd for i=1:dims
                    dist_qr += (Q.data[i, iq] - R.data[i, ir])^2
                end
                dist_qr = sqrt(dist_qr)

                if dist_qr < nn.distance
                    nn.distance = dist_qr
                    nn.a = qq
                    nn.b = rr
                end
            end
            Q.dQ = max(Q.dQ, nn.distance)
        end
        return
    end

    if is_leaf(Q)
        find_component_neighbors(Q, R.left, forest,  C)
        find_component_neighbors(Q, R.right, forest,  C)
    elseif is_leaf(R)
        find_component_neighbors(Q.left, R, forest,  C)
        find_component_neighbors(Q.right, R, forest,  C)
    else
        find_component_neighbors(Q.left, R.left, forest,  C)
        find_component_neighbors(Q.right, R.left, forest,  C)
        find_component_neighbors(Q.left, R.right, forest,  C)
        find_component_neighbors(Q.right, R.right, forest,  C)
    end
    Q.dQ = max(Q.left.dQ, Q.right.dQ)
end


"""
compute min dist. between bounding boxes, i.e. between rectangular boxes Q/R with
"""
function distance(q::DTBNode, r::DTBNode)
    dim ::Int64 = length(q.lb)
    dist :: Float64 = 0
    for i=1:dim
        delta = max(q.lb[i], r.lb[i]) - min(q.ub[i], r.ub[i])
        if delta> 0
            dist += delta^2
        end
    end
    dist = sqrt(dist)

    return dist
end

function plot_boxes!(Q::DTBNode, p)
    if !is_leaf(Q)
        plot_boxes!(Q.left, p)
        plot_boxes!(Q.right, p)
        return
    end
    plot!(p, [Q.lb[1], Q.ub[1], Q.ub[1], Q.lb[1], Q.lb[1]],[Q.lb[2], Q.lb[2], Q.ub[2], Q.ub[2], Q.lb[2]], lc=:black)
end

function plot_debug(x, forest, edges::Set{Edge})
    p = plot(legend=:none, aspect_ratio=:equal)

    # nodes
    roots = unique(forest.parents)
    d = Dict(a => Vector() for a in roots)

    for (i, xx) in enumerate(eachcol(x))
        r = find_root!(forest, i)
        push!(d[r], xx)
    end

    for component in values(d)
        if !(component isa Vector)
            continue
        end
        component = hcat(component...)
        scatter!(component[1,:], component[2,:], ms=3.0)
    end

    # edges
    m_edges = to_matrix(collect(edges))
    for zr in 1:size(m_edges, 1)
        plot!(x[1,m_edges[zr,:]], x[2,m_edges[zr,:]], linecolor=:gray)
    end

    return p
end

function make_tree(N :: Int64)
    points :: Vector{Vector{Float64}} = []

    for layer=0:(N-1)
        for i=1:(2^layer)
            push!(points, [Float64(i - 2^layer / 2 - 0.5) * 2^(N-layer), Float64(2^(N-layer) / 2)])
        end
    end

    return hcat(points...)
end