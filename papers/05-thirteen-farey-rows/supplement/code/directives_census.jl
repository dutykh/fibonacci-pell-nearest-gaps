# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
#
# Second, independent verification of the finite facts on directives.
#
# Every continued fraction [2; B_1, ..., B_j, 2] with blocks B_i = (2) or (1, 1) and numerator at most
# N = 10^12 is enumerated: the blocks are the matrices T = [2 1; 1 0] and A = [2 1; 1 1], a product W of
# blocks gives the numerator m = (2, 1) W (2, 1)' and the denominator n = first entry of W (2, 1)', and
# the numerator grows with every block, so the depth-first enumeration stops when it exceeds N.
# The Markov numbers up to N and their roots are computed independently, from the Cohn matrices
# mu(a) = [2 1; 1 1], mu(b) = [5 2; 2 1] multiplied along the Farey tree: the label is the upper right
# entry of the product, and the root is its upper left entry reduced modulo the label below m/2.
# For each Markov number m >= 5, the denominators n found for the numerator m are listed; the article
# asserts that the root is the only one.  The integers 31 and 1130 and the label of 9/10 are recorded
# as well.  All arithmetic is exact (Int128).
#
# Usage, from the supplement directory (about ten seconds and 0.3 GB on one thread):
#     julia --threads=1 code/directives_census.jl
# It writes data/directives-census.txt, which checks/check_directives.py compares with its own results.

const N = Int128(10)^12
const OUT = joinpath(@__DIR__, "..", "data", "directives-census.txt")

mul(X, Y) = (X[1]*Y[1] + X[2]*Y[3], X[1]*Y[2] + X[2]*Y[4], X[3]*Y[1] + X[4]*Y[3], X[3]*Y[2] + X[4]*Y[4])

function markov_roots(bound::Int128)
    roots = Dict{Int128,Vector{Int128}}()
    stack = Tuple{NTuple{4,Int128},NTuple{4,Int128}}[((Int128(2), Int128(1), Int128(1), Int128(1)),
                                                       (Int128(5), Int128(2), Int128(2), Int128(1)))]
    while !isempty(stack)
        L, R = pop!(stack)
        M = mul(L, R)
        m = M[2]
        m > bound && continue
        a = mod(M[1], m)
        push!(get!(roots, m, Int128[]), min(a, m - a))
        push!(stack, (L, M))
        push!(stack, (M, R))
    end
    roots
end

function census(bound::Int128, wanted::Set{Int128})
    found = Dict{Int128,Vector{Int128}}()
    count = 0
    stack = NTuple{4,Int128}[(Int128(1), Int128(0), Int128(0), Int128(1))]
    while !isempty(stack)
        w11, w12, w21, w22 = pop!(stack)
        m = 4w11 + 2w12 + 2w21 + w22
        m > bound && continue
        count += 1
        if m in wanted
            push!(get!(found, m, Int128[]), 2w11 + w12)
        end
        push!(stack, (2w11 + w12, w11, 2w21 + w22, w21))
        push!(stack, (2w11 + w12, w11 + w12, 2w21 + w22, w21 + w22))
    end
    found, count
end

function cfsum(p::Int128, q::Int128)
    s = Int128(0)
    while q != 0
        s += div(p, q)
        p, q = q, mod(p, q)
    end
    s
end

function label_and_root(P::Int, Q::Int)
    X = (Int128(1), Int128(0), Int128(0), Int128(1))
    for i in 1:Q
        X = mul(X, div(i * P, Q) - div((i - 1) * P, Q) == 1 ? (Int128(5), Int128(2), Int128(2), Int128(1)) :
                                                            (Int128(2), Int128(1), Int128(1), Int128(1)))
    end
    m = X[2]
    a = mod(X[1], m)
    m, min(a, m - a)
end

function main()
    roots = markov_roots(N)
    labels = sort([m for m in keys(roots) if m >= 5])
    wanted = Set{Int128}(labels)
    push!(wanted, Int128(31))
    push!(wanted, Int128(1130))
    found, count = census(N, wanted)
    m910, u910 = label_and_root(9, 10)
    n910 = Int128(2561343)
    open(OUT, "w") do io
        println(io, "# bound: ", N)
        println(io, "# continued fractions [2; B, 2]: ", count)
        println(io, "# Markov numbers from 5: ", length(labels))
        println(io, "# 31: ", join(sort(get(found, Int128(31), Int128[])), ","))
        println(io, "# 1130: ", join(sort(get(found, Int128(1130), Int128[])), ","))
        println(io, "# 9/10: ", m910, " ", u910, " ", cfsum(m910, u910) - 1, " ", n910, " ", cfsum(m910, n910) - 1)
        for m in labels
            length(roots[m]) == 1 || error("the label $m occurs more than once")
            println(io, m, " ", roots[m][1], " ", join(sort(get(found, m, Int128[])), ","))
        end
    end
    exceptions = count_exceptions(labels, roots, found)
    println("bound $N: $(length(labels)) Markov numbers from 5, $count continued fractions, $exceptions exceptions")
    exceptions == 0 || exit(1)
end

count_exceptions(labels, roots, found) = sum(m -> sort(get(found, m, Int128[])) == roots[m] ? 0 : 1, labels)

main()
