# Autor Adam Jam 
# index 268423

module blocksys
export MatrixContainer, read_regular_matrix, read_matrix, insert, at, read_vector,generate_right_sides,square_line_multiply, triangulate!, triangulate_and_decompose!, triangulate_with_partial_choice!, gaussian_elimination, gaussian_elimination!, gaussian_elimination_with_choice, save_right_sides, save_right_sides,LU_decomposition_solve

struct MatrixContainer
    n::Int  # size of the matrix  
    l::Int  # size of submatrices
    v::Int  # number of submatrices
    C_matrices::Vector{Matrix{Float64}}     
    A_matrices::Vector{Matrix{Float64}}
    B_lists::Vector{Vector{Vector{Float64}}}
    D_lists::Vector{Vector{Vector{Float64}}}
    function MatrixContainer(n::Int, l::Int)
        v = trunc(Int, n/l);
        if v != n/l
            throw(DomainError("INVALID BLOCK MATRIX SIZE"))
        end


        # [ A C D 0 0 0 ... 0 ]
        # [ B A C D 0 0 ... 0 ]
        # [ 0 B A C D 0 ... 0 ]
        # [ 0 0 B A C D ... 0 ]
        # [ .               . ]
        # [ .            .  0 ]
        # [ .             . D ]
        # [ .              .C ]
        # [ .      ...  0 B A ]

        C_matrices = [zeros(Float64, l, l) for _ in 1:(v-1)]
        A_matrices = [zeros(Float64, l, l) for _ in 1:v]
        # B_lists = [[zeros(Float64, l), zeros(Float64, l)] for _ in 1:(v-1)]
        B_lists = [[zeros(Float64, 2) for _ in 1:l] for _ in 1:(v-1)]
        D_lists = [[zeros(Float64, l) for _ in 1:2] for _ in 1:(v-2)]
        new(n, l, v, C_matrices, A_matrices, B_lists, D_lists)
    end

    function MatrixContainer(filename::String)
        return read_matrix(filename)        
    end

end


function read_regular_matrix(filename::String)::Matrix{Float64}
    new_matrix = nothing
    open(filename) do file
        first_line = split(readline(file))
        if (size(first_line)[1] != 2)
            throw(DomainError("INVALID BLOCK MATRIX SIZE"))
        end
        # println("$(first_ line[1]) $(first_line[2])")

        n = parse(Int, first_line[1])
        new_matrix = zeros(n,n)

        for line in eachline(file)
            args = split(line)
            if (size(args)[1] == 3) 
                new_matrix[parse(Int, args[1]), parse(Int, args[2])] = parse(Float64, args[3])
            else
                throw(DomainError("invalid input line $line"))
            end
            # println("Read line: $line");  # Example processing
        end
    end

    return new_matrix
end


function read_matrix(filename::String)::MatrixContainer
    
    new_matrix = nothing
    open(filename) do file
        first_line = split(readline(file))
        if (size(first_line)[1] != 2)
            throw(DomainError("INVALID BLOCK MATRIX SIZE"))
        end
        # println("$(first_ line[1]) $(first_line[2])")

        new_matrix = MatrixContainer(parse(Int, first_line[1]), parse(Int, first_line[2]))

        for line in eachline(file)
            args = split(line)
            if (size(args)[1] == 3) 
                insert(new_matrix, parse(Int, args[1]), parse(Int, args[2]), parse(Float64, args[3]))
            else
                throw(DomainError("invalid input line $line"))
            end
            # println("Read line: $line");  # Example processing
        end
    end

    return new_matrix
end

function insert(mc::MatrixContainer, i::Int, j::Int, value::Float64)
    big_i = Int(ceil(i / mc.l))
    little_i = i % mc.l
    if (little_i == 0)
        little_i = mc.l
    end
    big_j = Int(ceil(j / mc.l))
    little_j = j % mc.l
    if (little_j == 0)
        little_j = mc.l
    end
    
    if (i > mc.n || j > mc.n || i < 1 || j < 1)
        throw(BoundsError("INVALID INSERTION OUT OF BOUNDS i = $i, j = $j"))
    end

    diagonal = big_i - big_j
    if (diagonal == 0) # member of A_matrix
        mc.A_matrices[big_i][little_i, little_j] = value
    elseif (diagonal == 1) # member of B_list
        if (little_j == mc.l)
            mc.B_lists[big_j][little_i][2] = value
        elseif (little_j == mc.l-1)
            mc.B_lists[big_j][little_i][1] = value
        else
            throw(DomainError("INVALID INSERTION $i $j"))
        end
    elseif (diagonal == -1) # member of C_matrix
        mc.C_matrices[big_i][little_i, little_j] = value
    elseif (diagonal == -2) # memeber of D_list
        if (little_i == mc.l)
            mc.D_lists[big_i][2][little_j] = value
        elseif (little_i == mc.l-1)
            mc.D_lists[big_i][1][little_j] = value
        else
            throw(DomainError("INVALID INSERTION $i $j"))
        end
    else
        throw(DomainError("INVALID INSERTION $i $j"))
    end
end

function at(mc::MatrixContainer, i::Int, j::Int)
    big_i = Int(ceil(i / mc.l))
    little_i = i % mc.l
    if (little_i == 0)
        little_i = mc.l
    end
    big_j = Int(ceil(j / mc.l))
    little_j = j % mc.l
    if (little_j == 0)
        little_j = mc.l
    end
    
    if (i > mc.n || j > mc.n || i < 1 || j < 1)
        throw(BoundsError("POSITION ($i, $j) OUT OF BOUNDS"))
    end

    diagonal = big_i - big_j
    if (diagonal == 0) # member of A_matrix
        return mc.A_matrices[big_i][little_i, little_j]
    elseif (diagonal == 1) # member of B_list
        if (little_j == mc.l)
            return mc.B_lists[big_j][little_i][2]
        elseif (little_j == mc.l-1)
            return mc.B_lists[big_j][little_i][1]
        else
            return 0.0
        end
    elseif (diagonal == -1) # member of C_matrices
        return mc.C_matrices[big_i][little_i, little_j]
    elseif (diagonal == -2) # memeber of D_list
        if (little_i == mc.l)
            return mc.D_lists[big_i][2][little_j]
        elseif (little_i == mc.l-1)
            return mc.D_lists[big_i][1][little_j]
        else
            return 0.0
        end
    else
        return 0.0
    end
end

function read_vector(filename::String)::Vector{Float64}
    new_vector = nothing
    open(filename) do file
        first_line = split(readline(file))
        if (size(first_line)[1] != 1)
            throw(DomainError("INVALID VECTOR SIZE"))
        end
        # println("$(first_ line[1]) $(first_line[2])")

        vec_size = parse(Int, first_line[1])

        vals = readlines(file)
        new_vector = [parse(Float64, val) for val in vals]
        resize!(new_vector, vec_size)
    end

    return new_vector
end

function generate_right_sides(mc::MatrixContainer)
    size = mc.n
    vec = ones(size)
    t_line_multiply(mc, vec)
end

function square_line_multiply(mc::MatrixContainer, v::Vector{Float64})::Vector{Float64}
    result = zeros(mc.n)
    for real_i in 1:(mc.n)
        # calculate result[real_i] value
        big_i = Int(ceil(real_i/mc.l))
        little_i = real_i % mc.l == 0 ? mc.l : real_i % mc.l
        
        # B block
        if (big_i > 1)
            for real_j in ((big_i - 1)*mc.l - 1):((big_i - 1)*mc.l)
                big_j = Int(ceil(real_j/mc.l))
                little_j = real_j % mc.l == 0 ? 2 : 1
                result[real_i] += mc.B_lists[big_j][little_i][little_j] * v[real_j]
            end
        end
        # A block
        for real_j in ((big_i - 1)*mc.l + 1):(big_i*mc.l)
            big_j = Int(ceil(real_j/mc.l))
            little_j = real_j % mc.l == 0 ? mc.l : real_j % mc.l 
            result[real_i] += mc.A_matrices[big_j][little_i, little_j] * v[real_j]
        end
        # C block
        if (big_i < mc.v)
            for real_j in ((big_i)*mc.l + 1):((big_i + 1)*mc.l)
                # big_j = Int(ceil(real_j/mc.l))
                little_j = real_j % mc.l == 0 ? mc.l : real_j % mc.l 
                result[real_i] += mc.C_matrices[big_i][little_i, little_j] * v[real_j]
            end         
        end
        
        # deprecated C block
        # if (big_i < mc.v)
            # real_j = real_i+mc.l
            # little_j = real_j % mc.l == 0 ? mc.l : real_j % mc.l
            # result[real_i] += mc.C_matrices[big_i][little_j] * v[real_j]
        # end

        # D block
        if big_i < mc.v - 1 && little_i >= mc.n - 1
            for real_j in ((big_i + 1)*mc.l + 1):((big_i + 2)*mc.l)
                big_j = Int(ceil(real_j/mc.l))
                little_j = real_j % mc.l == 0 ? mc.l : real_j % mc.l 
                result[real_i] += mc.D_lists[big_j][little_i - n + 2][little_j] * v[real_j]
            end
        end
    end
        
    return result
end

function triangulate_with_partial_choice!(mc::MatrixContainer, b::Vector{Float64})
    for real_j in 1:(mc.n - 1)
        big_j = Int(ceil(real_j/mc.l))
        little_j = real_j % mc.l == 0 ? mc.l : (real_j % mc.l)

        # THE (PARTIAL) CHOICE 
        max_row_arg = real_j;
        diagonal = mc.A_matrices[big_j][little_j, little_j]
        current_max = diagonal
        for real_i in (real_j + 1):(mc.n)
            big_i = Int(ceil(real_i/mc.l))
            little_i = real_i % mc.l == 0 ? mc.l : (real_i % mc.l)
            
            if ((big_i - big_j == 1 && little_j < mc.l - 1) || big_i - big_j > 1) # we hit first zero below diagonal (j,j)
                break; # we hit first set zero
            end
            if current_max < at(mc, real_i, real_j)
                current_max = at(mc, real_i, real_j)
                max_row_arg = real_i
            end
        end
        if max_row_arg != real_j
            # swap rows 
            # row_j <-> row_{max_row_arg}

            # tmp_copy = [at(mc, real_j,p) for p in real_j:min(mc.n, real_j + mc.l * 3)]

            # swap with something within the same block
            for real_k in real_j:(min(mc.n, (big_j + 1) * mc.l)) 
                #stop iterating at the last column of C block
                tmp_jk = at(mc, real_j, real_k)
                insert(mc, real_j, real_k, at(mc, max_row_arg, real_k))
                insert(mc, max_row_arg, real_k , tmp_jk)
            end
            if (little_j > mc.l - 2) # swap with presumably something outside 
                # iterate through a D block
                for real_k in ((big_j + 1) * mc.l + 1):(min(mc.n, (big_j + 2) * mc.l))
                    tmp_jk = at(mc, real_j, real_k)
                    insert(mc, real_j, real_k, at(mc, max_row_arg, real_k))
                    insert(mc, max_row_arg, real_k , tmp_jk)
                end
            end

            b[max_row_arg], b[real_j] = b[real_j], b[max_row_arg]

            # if (tmp_copy != [at(mc, max_row_arg,p) for p in real_j:min(mc.n, real_j + mc.l * 3)])
            #     println("at row $real_j and $max_row_arg")
            #     println(tmp_copy)
            #     println([at(mc, max_row_arg, p) for p in real_j:min(mc.n, real_j + mc.l * 3)])
            #     throw(ErrorException("swap operantion failed"))
            # end
        end
        # END OF PARTIAL CHOICE

        if (mc.A_matrices[big_j][little_j, little_j] == 0)
            throw(DomainError("triangulation failed"))
        end

        for real_i in (real_j + 1):(mc.n)
            big_i = Int(ceil(real_i/mc.l))
            little_i = real_i % mc.l == 0 ? mc.l : (real_i % mc.l)
            
            if ((big_i - big_j == 1 && little_j < mc.l - 1) || big_i - big_j > 1) # we hit first zero below diagonal (j,j)
                break; # we hit first set zero
            end

            t = at(mc, real_i, real_j) / mc.A_matrices[big_j][little_j, little_j]
            for real_k in real_j:min((big_j + 1) * mc.l, mc.n)
                value = at(mc, real_i, real_k) - t * at(mc, real_j, real_k)
                insert(mc, real_i, real_k, value) 
            end
            if (little_j > mc.l - 2) # swap with presumably something outside 
                # iterate through a D block
                for real_k in (min((big_j + 1) * mc.l, mc.n) + 1):(min((big_j + 2) * mc.l, mc.n))
                    value = at(mc, real_i, real_k) - t * at(mc, real_j, real_k)
                    insert(mc, real_i, real_k, value) 
                end
            end

            b[real_i] = b[real_i] - b[real_j] * t
        end
    end        
end

function triangulate!(mc::MatrixContainer, b::Vector{Float64})
    for real_j in 1:(mc.n - 1)
        big_j = Int(ceil(real_j/mc.l))
        little_j = real_j % mc.l == 0 ? mc.l : (real_j % mc.l)

        if (mc.A_matrices[big_j][little_j, little_j] == 0)
            throw(DomainError("triangulation failed"))
        end

        for real_i in (real_j + 1):(mc.n)
            big_i = Int(ceil(real_i/mc.l))
            little_i = real_i % mc.l == 0 ? mc.l : (real_i % mc.l)
            
            if ((big_i - big_j == 1 && little_j < mc.l - 1) || big_i - big_j > 1) # we hit first zero below diagonal (j,j)
                break; # we hit first set zero
            end

            t = at(mc, real_i, real_j) / mc.A_matrices[big_j][little_j, little_j]
            for real_k in real_j:(min(real_j + mc.l, mc.n))
                value = at(mc, real_i, real_k) - t * at(mc, real_j, real_k)
                # println("$real_i $real_k $(t) $(at(mc, real_i, real_k)) $(t * at(mc, real_j, real_k)) $value")
                insert(mc, real_i, real_k, value) 
            end

            b[real_i] = b[real_i] - b[real_j] * t
        end
    end
end


function gaussian_elimination_with_choice(mc::MatrixContainer, b::Vector{Float64})
    result = zeros(mc.n)
    U = deepcopy(mc)
    b = deepcopy(b)
    # MAKE'EM TRIANGULAR
    triangulate_with_partial_choice!(U, b)
    
    result[mc.n] = b[mc.n] / at(U, mc.n, mc.n) # calculate first result value
    for i in (U.n-1):-1:1
        sum = 0
        for j in (i+1):min(i + 1 + mc.l * 2, mc.n) # missed optimization with bounds
            sum += result[j] * at(U, i, j)
        end
        result[i] = (b[i] - sum) / at(U, i, i) 
    end

    return result
end

function gaussian_elimination(mc::MatrixContainer, b::Vector{Float64})
    result = zeros(mc.n)
    U = deepcopy(mc)
    b = deepcopy(b)

    # MAKE'EM TRIANGULAR
    triangulate!(U, b)
    
    result[mc.n] = b[mc.n] / at(U, mc.n, mc.n)
    for i in (U.n-1):-1:1
        sum = 0
        for j in (i+1):min(mc.n, i + 1 + mc.l)
            sum += result[j] * at(U, i, j)
        end
        result[i] = (b[i] - sum) / at(U, i, i) 
    end

    return result
end

function gaussian_elimination!(mc::MatrixContainer, b::Vector{Float64})
    result = zeros(mc.n)

    # MAKE'EM TRIANGULAR
    triangulate!(mc, b)
    
    result[mc.n] = b[mc.n] / at(mc, mc.n, mc.n)
    for i in (mc.n-1):-1:1
        sum = 0
        for j in (i+1):min(mc.n, i + 1 + mc.l)
            sum += result[j] * at(mc, i, j)
        end
        result[i] = (b[i] - sum) / at(mc, i, i) 
    end

    return result
end



function save_right_sides(vec::Vector{Float64}, outputfile::String)
    open(outputfile, "w") do file
        for val in vec
            write(file, "$val\n")
        end
    end
end

function save_right_sides(vec::Vector{Float64}, outputfile::String, error::Float64)
    open(outputfile, "w") do file
        write(file, "$error\n")
        for val in vec
            write(file, "$val\n")
        end
    end
end





# A = LU decomposition
function LU_decomposition_solve(mc::MatrixContainer, b::Vector{Float64})::Vector{Float64}
    mc = deepcopy(mc)
    y = deepcopy(b)
    g = deepcopy(b)

    # to get U = A^(n) we triangulate it as before
    # but we keep track of eliminations and store them into L_inv
    # as confusing as it might be L_inv is actually
    # a series of L^(k) k = 1 : n-1 operations
    # (as used in simple triangulation) 
    # and L is a reversed order series of L^(k)^(1) operations 
    # so we can obtain pretty looking equation "A = LU" at the end

    # triangulation (U):
    # println("before: $mc")
    (U, L_inv, g) = triangulate_and_decompose!(mc, g)
    # println(L_inv)
    # println("after: $mc")

    # y = square_line_multiply(L_inv, b)
    # compute L_inv * b 
    # iterate diagonaly then top-down
    for i in 1:(mc.n)
        for j in 1:(mc.l + 1)
            if (i + j > mc.n)
                break
            end
            y[i + j] += b[i] * L_inv[i][j] # missed optimization to omit repeated multiplication by b[i]  
        end
    end
    # println("original b = $b")
    # println("y = $y")
    # println("what y should be = $g")

    result = zeros(mc.n)
    
    result[U.n] = g[U.n] / at(U, U.n, U.n)
    for i in (U.n-1):-1:1 
        sum = 0
        for j in (i+1):(min(i + U.l, U.n)) # at most l values to the right might be nonzero
            sum += result[j] * at(U, i, j)
        end
        result[i] = (g[i] - sum) / at(mc, i, i) 
    end

    return result
end

function triangulate_and_decompose!(mc::MatrixContainer, b::Vector{Float64})
    # Returns the LU decomposition
    # turns mc into U matrix
    # keeps track of eliminations in L_inv
    L_inv = [ zeros(1 + mc.l) for _ in 1:(mc.n) ]
    # one could toil some more and add necessary E_lists to the matrix container so that partial choice could be implemented to the LU decomposition
    # as well as L_inv to be lists of size 2l for every last two columns of a block
    

    for real_j in 1:(mc.n - 1)
        big_j = Int(ceil(real_j/mc.l))
        little_j = real_j % mc.l == 0 ? mc.l : (real_j % mc.l)

        if (mc.A_matrices[big_j][little_j, little_j] == 0)
            throw(DomainError("triangulation failed"))
        end

        for real_i in (real_j + 1):(mc.n)
            big_i = Int(ceil(real_i/mc.l))
            little_i = real_i % mc.l == 0 ? mc.l : (real_i % mc.l)
            
            if ((big_i - big_j == 1 && little_j < mc.l - 1) || big_i - big_j > 1) # we hit first zero below diagonal (j,j)
                break; # we hit first set zero
            end

            t = at(mc, real_i, real_j) / mc.A_matrices[big_j][little_j, little_j]
            L_inv[real_j][real_i - real_j] = -t
            
            for real_k in real_j:(min(real_j + mc.l, mc.n))
                value = at(mc, real_i, real_k) - t * at(mc, real_j, real_k)
                # println("$real_i $real_k $(t) $(at(mc, real_i, real_k)) $(t * at(mc, real_j, real_k)) $value")
                insert(mc, real_i, real_k, value) 
            end

            b[real_i] = b[real_i] - b[real_j] * t
        end

    end

    # copy L_inv onto the zeros of U and negate it
    # since L_inv is L with swapped signs by doing so we obtain LU
    
    #traverse diagonally then top-down omit diagonals
    # for j in 1:(mc.n)
        # for i in (j+1):(mc.n)
            # to do LU in one matrixContainer we would've needed to add full B_matrices instead of B_lists
            # insert(mc, i, j, -at(L_inv, i, j))
            # insert(mc, i, j, -L_inv[i-j,j])
        # end
    # end

    return mc, L_inv, b # = U, L_inv 
end


end #module end



# TESTS 
A = MatrixContainer("ON_LAB5/Dane/Dane0/A.txt")
triangulate!(A)
for i in 1:A.n
    for j in i:A.n
        println(at(A, i, j))
    end
end
for i in 1:A.n
    for j in 1:(i-1)
        if (at(A, i, j) != 0)
            println("non-zero val at ($i, $j) = $(at(A, i, j))")
        end
    end
end


# TEST ACCURACY WITH gaussian_elimination
# A = MatrixContainer("ON_LAB5/Dane/Dane1/A.txt")
# b = read_vector("ON_LAB5/Dane/Dane1/b.txt")
# generated_b = generate_right_sides(A)
# for _ in 1:100
        



# TEST square_line_multiply multiplying MatrixContainer by a vector
A = MatrixContainer("ON_LAB5/Dane/Dane1/A.txt")
A.n
A.v
reg_mat = read_regular_matrix("ON_LAB5/Dane/Dane1/A.txt")
accumulator = 0
for _ in 1:100
    random_vec = rand(A.n) * 100
    difs = abs.(reg_mat * random_vec - square_line_multiply(A,random_vec))
    difs = difs./(square_line_multiply(A,random_vec))
    accumulator += sum(difs)
end
println("error: $accumulator")



# TRST I/O
save_right_sides([0.0,1.0,2.0], "mock.txt", 12.3)
# function read_vector(filename::String)

# end


b = read_vector("ON_LAB5/Dane/Dane1/b.txt")
A = MatrixContainer("ON_LAB5/Dane/Dane1/A.txt")
A = MatrixContainer("ON_LAB5/Dane/Dane0/A.txt")
A = MatrixContainer("ON_LAB5/Dane/Dane10000/A.txt")
square_line_multiply(A,b)
reg_mat = read_regular_matrix("ON_LAB5/Dane/Dane1/A.txt")
reg_mat = read_regular_matrix("ON_LAB5/Dane/Dane0/A.txt")
reg_mat * b - square_line_multiply(A,b)
random_vector = rand(A.n)
dif = reg_mat * (ones(A.n) + random_vector) - square_line_multiply(A,ones(A.n) + random_vector)
println(dif)
println(reg_mat)
println(A.A_matrices)
println(A.B_lists)
println(A.C_matrices)



# GAUSSIAN ELIMINATION TESTS

A = MatrixContainer("ON_LAB5/Dane/Dane0/A.txt")
U = deepcopy(A)
A.C_matrices[1][1,1] = 2
A.C_matrices[1][1,1]
U
at(U, 1,1)
at(A, 1,1)
insert(U, 1, 1, 2.0)




println(A_mat \ generated_b)
println(res)
println(res2)
println(A.A_matrices)
println(A.B_lists)
println(A.C_matrices)

A = MatrixContainer("ON_LAB5/Dane/Dane1/A.txt")
b = read_vector("ON_LAB5/Dane/Dane1/b.txt")
println(triangulate!(A, b))
println(A)
println(b)
A = MatrixContainer("ON_LAB5/Dane/Dane1/A.txt")
b = read_vector("ON_LAB5/Dane/Dane1/b.txt")
triangulate_and_decompose!(A)
println(A)



b = read_vector("ON_LAB5/Dane/Dane1/b.txt")
A_mat = read_regular_matrix("ON_LAB5/Dane/Dane1/A.txt")


test_matrices_data = [
    "ON_LAB5/Dane/Dane1/A.txt",
    "ON_LAB5/Dane/Dane10000/A.txt",
    "ON_LAB5/Dane/Dane50000/A.txt",
    "ON_LAB5/Dane/Dane100000/A.txt",
    "ON_LAB5/Dane/Dane300000/A.txt",
    "ON_LAB5/Dane/Dane500000/A.txt",
]

test_vector_data = [
    "ON_LAB5/Dane/Dane1/b.txt",
    "ON_LAB5/Dane/Dane10000/b.txt",
    "ON_LAB5/Dane/Dane50000/b.txt",
    "ON_LAB5/Dane/Dane100000/b.txt",
    "ON_LAB5/Dane/Dane300000/b.txt",
    "ON_LAB5/Dane/Dane500000/b.txt",
]

# GENERATED RIGHT SIDES ACCURACY TESTS

for test_matrix in test_matrices_data
    println("accuracy test: $test_matrix")
    A = MatrixContainer(test_matrix);
    b = generated_b = generate_right_sides(A);
    
    res  = gaussian_elimination(A, generated_b);
    res2 = gaussian_elimination_with_choice(A, generated_b);

    sum = 0
    for i in 1:res.size[1]
        sum += abs(res[i] - 1)
    end
    println("no choice: $(sum / res.size[1])")

    sum = 0
    for i in 1:res2.size[1]
        sum += abs(res2[i] - 1)
    end
    println("partial choice: $(sum / res2.size[1])")
end


# time TESTS
using BenchmarkTools
using Plots
for test_matrix in test_matrices_data
    println("test: $test_matrix")
    A = MatrixContainer(test_matrix);
    b = generated_b = generate_right_sides(A);
    
    res  = @time gaussian_elimination(A, generated_b);
    res2 = @time gaussian_elimination_with_choice(A, generated_b);
end