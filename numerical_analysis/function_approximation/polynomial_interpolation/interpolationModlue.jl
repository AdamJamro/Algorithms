# author
# Adam Jamrozinski
# 268423

module InterpolationModlue

export ilorazyRoznicowe, warNewton, naturalna, rysujNnfx

function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    if (x.size != f.size)
        return throw(DomainError(x, "x and f.(x) vectors size does not match!"));
    end
    n = x.size[1] - 1;
    quotients = zeros(n+1);
    quotients[1] = f[1];

    for i in 1:n
        for j in 1:(n+1-i)
            f[j] = (f[j] - f[j+1])/(x[j] - x[j+i]);
        end
        quotients[i+1] = f[1];
    end

    return quotients;
end


function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64)
    if (x.size != fx.size)
        return throw(DomainError(x, "x and fx vectors size does not match!"));
    end
    n = x.size[1] - 1;

    result = fx[n+1];
    for k in n:-1:1
        result = result * (t-x[k]) + fx[k];
    end
    
    return result;
end

function naturalna(x::Vector{Float64}, fx::Vector{Float64})
    if (x.size != fx.size)
        return throw(DomainError(x, "x and fx vectors size does not match!"));
    end

    n = length(fx) - 1; # polynomial degree
    coeffs = zeros(n + 1);

    for i in (n+1):-1:1
        for j in i:n
            coeffs[j] = coeffs[j+1] - x[i] * coeffs[j] 
        end
        coeffs[n+1] = fx[i] - coeffs[n+1] * x[i];
    end

    return coeffs
end

using Plots
function rysujNnfx(f, a::Float64, b::Float64, n::Int)
    if (n < 1)
        throw(DomainError(n, "n must be positive"))
    end

    h = (b-a)/n;
    x = [i for i in a:h:b];

    if (any(isinf, f.(x)))
        throw(DomainError(a, "unable to create $n points on [$a,$b]"));
    end

    quotients = ilorazyRoznicowe(x, f.(x));
    args = [i for i in a:0.1:b];
    vals = [warNewton(x, quotients, i) for i in args];


    plot1 = plot(
    args, vals,
    label="wielomian interpolacyjny", 
    xlabel="x", ylabel="N_$n(x)", 
    color=:"red",
    lw=2
    )

    plot!(plot1,
        args, f.(args),
        label="funkcja interpolowana",
    )

    scatter!(plot1,
        x, f.(x), 
        label = "węzły", 
        # Optional: customize point appearance
        marker = :circle,  # point shape
        markercolor = :blue,
        markersize = 5
    )
    
    display(plot1);
end

end # Interpolation module

