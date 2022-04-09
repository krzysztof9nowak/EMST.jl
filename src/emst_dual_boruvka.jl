
using DataStructures
using Distances
using Statistics

import Base.isequal
import Base.hash



function plot_debug(x, forest, edges::Set{Edge})
    p = plot(legend=:none)

    # points

    roots = unique(forest.parents)
    d = Dict(a => Vector() for a in roots)

    for (i, xx) in enumerate(eachcol(x))
        r = find_root!(forest, i)
        push!(d[r], xx)
    end

    for component in values(d)
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
    for zi = 1:length(ee); m[zi,1] = ee[zi].a; m[zi,2] = ee[zi].b; end
return m
end

"""
  DTBNode

implements node of a KD tree
"""
mutable struct DTBNode
    id::Int64 # identifies the node ; root has 00..00 , then
              # levels (right/left) are encoded as 01 or 10
              # i.e. node root->"right"->"right"->"left" ends with:
              # 10 01 01 00
              #  l r  r  R

    # data / indidces of this kd tree node
    data::Array{Float64,2}
    subset::Array{Int64}

    # kd tree node bounds
    lb::Array{Float64,1}
    ub::Array{Float64,1}

    # max_{v in this node} ( length of candidate edge of v)
    dQ::Float64

    left::DTBNode
    right::DTBNode

    # list of forest roots
    roots::Vector{Int64}

    function DTBNode(id::Int64, data::Array{Float64,2}, subset::Array{Int64}, box_lb::Array{Float64,1}, box_ub::Array{Float64,1}, dQ::Float64, roots::Vector{Int64})
        n = new(id, data, subset, box_lb, box_ub, dQ)
        n.left = n
        n.right = n
        n.roots = roots
        return n
    end
end
function isequal(a::DTBNode, b::DTBNode)
    return a.id == b.id
end
function hash(a::DTBNode)
    # a.id
    hash(a.id) # we could also drop this hash.. should not matter..
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
    root = DTBNode(Int64(0), xx, collect(Int64(1):Int64(size(xx, 2))), fill!(ones(size(xx, 1)), -Inf), fill!(ones(size(xx, 1)), Inf), Inf, collect(Int64(1):Int64(size(xx, 2))))
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

    id::Int64 = node.id
    id_depth  = Int(ceil((64 - leading_zeros(id)) / 2)) + 1 # in this cell we are, i.e. we have to shift id_depth times left by two bits..
    id_l = id   |  (1) << (2 * id_depth)
    id_r = id   |  (2) << (2 * id_depth)


    node.left = DTBNode(id_l, data_a, range_a, box_lb_a, box_ub_a, Inf, range_a)
    node.right = DTBNode(id_r, data_b, range_b, box_lb_b, box_ub_b, Inf, range_b)

    kdtree_split!(node.left, nmin)
    kdtree_split!(node.right, nmin)
end

"""
  compute_emst(data::Array{Float64,2};nmin::Int64=64)

Computes EMST for the given data (where columns are samples).
nmin is the max number of elements in kd-tree node.
"""
function compute_emst(data::Array{Float64,2};nmin::Int=64)
    nmin64 = Int64(nmin)
    root = kdtree(data)
    kdtree_split!(root, nmin64)
    edges = dtb(root, IntDisjointSets(size(data, 2)))
    return to_matrix(collect(edges))
end


"""
  dtb(q::DTBNode,e::IntDisjointSets)

implements the dual-tree Boruvka algorithm which computes the EMST for the
given KD-tree.
"""
function dtb(Q::DTBNode, forest::IntDisjointSets)
    edges = Set{Edge}()

    while forest.ngroups > 1
        ngroups = forest.ngroups;
        println("--> ngroups: $ngroups")

        # prepare dicts for candidate edges..
        C_dcq = Dict{Int64,Float64}()
        C_e   = Dict{Int64,Edge}()
        # init themz
        roots = unique(forest.parents)
        for ri in roots
            C_dcq[ri] = Inf
        end

        @time find_component_neighbors(Q, Q, forest, C_dcq, C_e)
        sleep(3)

        # and now add the edges..
        for ne::Edge in values(C_e)
            union!(forest, ne.a, ne.b)
            push!(edges, ne)
        end

        @time update_node_roots(Q, forest)

        p = plot_debug(Q.data, forest, edges)
        gui(p)
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


"""
Find Component Neighbors

Components are identified by the root of the component
C_dcq : component distances to candidate edges
C_e   : component candidate edges (comp i candidate edge is C_e[:,i])
"""
function find_component_neighbors(Q::DTBNode, R::DTBNode, forest::IntDisjointSets, C_dcq::Dict{Int64,Float64}, C_e::Dict{Int64,Edge})
    # check R and Q contain only one tree
    if length(Q.roots) == 1 && length(R.roots) == 1 && R.roots[1] == Q.roots[1]
        return
    end

    if distance(Q, R) > Q.dQ
        return
    end

    if is_leaf(Q) && is_leaf(R)
        n_dQ::Float64 = Q.dQ

        pairwise_d = Distances.pairwise(Euclidean(), Q.data, R.data, dims=2) # O(nmin^2)
        i=0
        @inbounds for (iq, qq) in enumerate(Q.subset), (ir, rr) in enumerate(R.subset) # O(nmin^2) * O(log n)
            if in_same_set(forest, qq, rr)
                continue
            end

            cq = find_root!(forest, qq) # tree of q, O(log n)

            dist_qr = pairwise_d[iq, ir]
            if dist_qr < C_dcq[ cq ]
                C_dcq[ cq ] = dist_qr
                C_e[ cq ]   = Edge(qq, rr)
                if dist_qr > n_dQ
                    n_dQ = dist_qr
                end
                i += 1
            end
        end
        Q.dQ = n_dQ

        return
    end

    find_component_neighbors(Q.left, R.left, forest,  C_dcq, C_e)
    find_component_neighbors(Q.right, R.left, forest,  C_dcq, C_e)
    find_component_neighbors(Q.left, R.right, forest,  C_dcq, C_e)
    find_component_neighbors(Q.right, R.right, forest,  C_dcq, C_e)
    Q.dQ = max(Q.left.dQ, Q.right.dQ)
end



function intervals_distance(xl::Float64, xu::Float64, yl::Float64, yu::Float64) :: Float64
    if xl <= yl <= xu || xl <= yu <= xu
        return 0
    end
    return min(abs(xl - yu), abs(yl - yu))
end

"""
compute min dist. between bounding boxes, i.e. between rectangular boxes Q/R with
"""
function distance(q::DTBNode, r::DTBNode)
    rdists = [intervals_distance(xl, xu, yl, yu) for (xl, xu, yl, yu) in zip(q.lb, q.ub, r.lb, r.ub)]
    return sqrt(sum(rdists.^2))
end
