
using DataStructures
using Distances
using Statistics

import Base.isequal
import Base.hash

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
    # splits::Array{DTBSplit,1}

    # kd tree node bounds
    lb::Array{Float64,1}
    ub::Array{Float64,1}

    # max_{v in this node} ( length of candidate edge of v)
    dQ::Float64

    left::DTBNode
    right::DTBNode

    # DTBNode(data::Array{Float64,2}) = (x = new(data); x.left=x; x.right=x)
    # function DTBNode(data::Array{Float64,2},subset::Array{Int64},splits::Array{DTBSplit,1},dQ::Float64)
    function DTBNode(id::Int64, data::Array{Float64,2}, subset::Array{Int64}, box_lb::Array{Float64,1}, box_ub::Array{Float64,1}, dQ::Float64)
        n = new(id, data, subset, box_lb, box_ub, dQ)
        n.left = n
        n.right = n
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
    root = DTBNode(Int64(0), xx, collect(Int64(1):Int64(size(xx, 2))), fill!(ones(size(xx, 1)), -Inf), fill!(ones(size(xx, 1)), Inf), Inf)
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
    range_b = node.subset[.~bx]

    data_a = node.data[:,bx]
    data_b = node.data[:,.~bx]

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


    node_left  = kdtree_split!(DTBNode(id_l, data_a, range_a, box_lb_a, box_ub_a, Inf), nmin)
    node_right = kdtree_split!(DTBNode(id_r, data_b, range_b, box_lb_b, box_ub_b, Inf), nmin)

    node.left  = node_left
    node.right = node_right

    return node
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
    return EMST.to_matrix(collect(edges))
end


"""
  dtb(q::DTBNode,e::IntDisjointSets)

implements the dual-tree Boruvka algorithm which computes the EMST for the
given KD-tree.
"""
function dtb(q::DTBNode, e::IntDisjointSets)
    edges = Set{Edge}()

    while (e.ngroups > 1)
        ngroups = e.ngroups;
        println("--> ngroups: $ngroups")

        # prepare dicts for candidate edges..
        C_dcq = Dict{Int64,Float64}()
        C_e   = Dict{Int64,Edge}()
        # init themz
        roots = unique(e.parents)
        for ri in roots
            C_dcq[ri] = Inf
        end

        # compute..
        # fcn(q,q,e,C_dcq,C_e , hashmap_QvsR )
        fcn(q, q, e, C_dcq, C_e)

        # and now add the edges..
        for ne::Edge in values(C_e)
            union!(e, ne.a, ne.b)
            push!(edges, ne)
    end
    end
    return edges
end


"""
  fcn(q::DTBNode,r::DTBNode,e::IntDisjointSets , C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,MyE} )
Find Component Neighbors

Components are identified by the root of the component
C_dcq : component distances to candidate edges
C_e   : component candidate edges (comp i candidate edge is C_e[:,i])
"""
function fcn(q::DTBNode, r::DTBNode, e::IntDisjointSets, C_dcq::Dict{Int64,Float64}, C_e::Dict{Int64,Edge})
    # check all in same component
    onecomp::Bool  = true
    joined         = [q.subset;r.subset]
    # println(joined)
    for ji in joined
        if (~in_same_set(e, joined[1], ji))
            onecomp = false
            break
    end
    end
    if (onecomp)
        return
    end

    # check d(Q,R) > d(Q)
    dqr = computeDQR(q, r)
    if ( dqr > q.dQ)
        return
    end

    # check if R and Q in a leaf node
    if ( is_leaf(q) && is_leaf(r) )
        # n_dQ = -Inf
        n_dQ::Float64 = q.dQ

        all_d_qr = Distances.pairwise(Euclidean(), q.data, r.data, dims=2)

        # --> hashing these values does not really speed up things..
        #     and it uses tons of memory, such that fancier handling of
        #     hash map entries would be required
        # if(!haskey(hashmap_QvsR,(q,r)))
        #    if(!haskey(hashmap_QvsR,(r,q)))
        #        all_d_qr = Distances.pairwise(Euclidean(),q.data,r.data)
        #        #hashmap_QvsR[(q,r)] = all_d_qr
        #    else
        #        all_d_qr = transpose(hashmap_QvsR[(r,q)])
        #    end
        #    hashmap_QvsR[(r,q)] = transpose(all_d_qr)
        # else
        #    all_d_qr = hashmap_QvsR[(q,r)]
        #    #print('.')
        # end

        for iq = 1:size(q.subset, 1)
            for ir = 1:size(r.subset, 1)
                qq = q.subset[iq]
                rr = r.subset[ir]
                if ( in_same_set(e, qq, rr) )
                    continue
                end

                cq = find_root(e, rr) # compoment of q

                # check distance:
                # dist_qr::Float64 = sqrt( sum( (q.data[:,iq]-r.data[:,ir]).^2 ) )
                dist_qr = all_d_qr[iq,ir]

                if ( dist_qr < C_dcq[ cq ] )
                    C_dcq[ cq ] = dist_qr
                    C_e[ cq ]   = Edge(qq, rr) # (qq,rr)
                    # and dQ !
                    n_dQ = max(n_dQ, dist_qr)
                    # n_dQ = dist_qr
            end
            end
        end
        q.dQ = n_dQ

        return
    end

    fcn(q.left, r.left, e,  C_dcq, C_e)
    fcn(q.right, r.left, e,  C_dcq, C_e)
    fcn(q.left, r.right, e,  C_dcq, C_e)
    fcn(q.right, r.right, e,  C_dcq, C_e)
    q.dQ = max(q.left.dQ, q.right.dQ)
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
function computeDQR(q::DTBNode, r::DTBNode)
    rdists = [intervals_distance(xl, xu, yl, yu) for (xl, xu, yl, yu) in zip(q.lb, q.ub, r.lb, r.ub)]
    return sqrt(sum(rdists.^2))
end
