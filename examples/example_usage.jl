using Base: Float64, Int64
# using EMST
include("../src/emst_dual_boruvka.jl")
using Plots

function plot_emst_2d(x, edges)
    p = plot(legend=:none)
    scatter!(x[1,:], x[2,:], ms=3.0, marker=(stroke(0, :gray)))
    for zr in 1:size(edges, 1)
        plot!(x[1,edges[zr,:]], x[2,edges[zr,:]], linecolor=:gray)
    end
    return p
end
function plot_emst_3d(x, edges)
    p = plot(legend=:none)
    scatter!(x[1,:], x[2,:], x[3,:], ms=3.0, marker=(stroke(0, :gray)))
    for zr in 1:size(edges, 1)
        plot!(x[1,edges[zr,:]], x[2,edges[zr,:]], x[3,edges[zr,:]], linecolor=:gray)
    end
    return p
end


x_test  = rand(2, 400)
@time e_test  = compute_emst(x_test;nmin=64);
# verify_emst(x_test,e_test,200)
plot_emst_2d(x_test,e_test)

# x_test  = rand(3,400)
# @time e_test  = EMST.compute_emst(x_test;nmin=64);
# verify_emst(x_test,e_test,200)
# plot_emst_3d(x_test,e_test)
