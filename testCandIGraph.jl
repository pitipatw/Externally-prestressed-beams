include("Interpolate.jl")

A = 0.:230.
list_c = zeros(length(A))
list_I = zeros(length(A))

for i in eachindex(A)
    Ai = A[i]
    c = get_C(Ai)
    I = get_Icrack(c)

    list_c[i] = c 
    list_I[i] = I 
end

#plot A vs list_c and list_I using Makie
using Makie , GLMakie
scene = Scene()
lines!(scene, A, list_c, color = :red)
# lines!(scene, A, list_I, color = :blue)
display(scene)

