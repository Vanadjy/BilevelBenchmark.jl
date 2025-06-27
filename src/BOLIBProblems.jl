export BilevelProblem, get_bilevel_problem

using LinearAlgebra

struct BilevelProblem
    name::String
    dim::Vector{Int}  # [n_x, n_y, n_G, n_g]
    xy0::Vector{Float64}
    Ff::Vector{Float64}  # [F, f, Status]
    F_func::Function
    f_func::Function
    G_func::Function
    g_func::Function
    sol::Union{Vector{Float64}, Float64}  # Solution or optimal value
end

function get_bilevel_problem(prob_no::Union{Int,String})
    ## ---------------- NONLINEAR BILEVEL PROBLEMS ---------------- ##
    if prob_no == 1 || prob_no == "AiyoshiShimizu1984Ex2"
        return BilevelProblem(
            "AiyoshiShimizu1984Ex2",
            [2, 2, 5, 6],
            [10.0, 10.0, 20.0, 20.0],
            [5.0, 0.0, 1.0],
            (x, y) -> 2*x[1] + 2*x[2] - 3*y[1] - 3*y[2] - 60,
            (x, y) -> (y[1] - x[1] + 20)^2 + (y[2] - x[2] + 20)^2,
            (x, y) -> [
                x[1] + x[2] + y[1] - 2*y[2] - 40,
                x[1] - 50,
                x[2] - 50,
                -x[1],
                -x[2]
            ],
            (x, y) -> [
                2*y[1] - x[1] + 10,
                2*y[2] - x[2] + 10,
                -y[1] - 10,
                -y[2] - 10,
                y[1] - 20,
                y[2] - 20
            ],
            [25.0, 30.0, 5.0, 10.0]
        )

    elseif prob_no == 2 || prob_no == "AllendeStill2013"
        return(BilevelProblem(
            "AllendeStill2013",
            [2, 2, 5, 2],
            [0.0, 0.0, 0.0, 0.0],
            [1.0, -0.5, 1.0],
            (x,y) -> (x[1])^2 - 2*x[1] + (x[2])^2 - 2*x[2] + y[1]^2 + y[2]^2,
            (x,y) -> (y[1])^2 - 2*x[1]*y[1] + y[2]^2 - 2*x[2]*y[2],
            (x,y) -> [
                - x[1],
                - y[1],
                - x[2],
                - y[2],
                x[1]-2
            ],
            (x,y) -> [
                (y[1] - 1)^2 - 0.25,
                (y[2] - 1)^2 - 0.25
            ],
            [0.5, 0.5, 0.5, 0.5])
            
        )
    elseif prob_no == 3 || prob_no == "AnEtal2009"
        return(BilevelProblem(
            "AnEtal2009",
            [2, 2, 6, 4],
            [1.0, 1.0, 1.0, 1.0],
            [2251.55, 565.78, 1.0],
            (x,y) -> (x[1]-5)^2 + (2x[2]+1)^2,
            (x,y) -> (y[1]-4)^2 + (y[2]-2)^4,
            (x,y) -> [
                4x[1]+5x[2] <= 20,
                -4x[1]-5x[2] <= -20,
                x[1]-2y[1]+4y[2] <= 10,
                -x[1]+2y[1]-4y[2] <= -10,
                x[2]-3y[1]+y[2] <= 10,
                -x[2]+3y[1]-y[2] <= -10
            ],
            (x,y) -> [
                y[1] <= 10,
                y[2] <= 10,
                -y[1] <= 0,
                -y[2] <= 0
            ],
            [0.200001, 1.999997, 3.999998, 4.600005])
        )
    
    elseif prob_no == 4 || prob_no == "Bard1988Ex1"
        return BilevelProblem(
            "Bard1988Ex1",
            [1, 1, 1, 4],
            [4.0, 0.0],
            [17.0, 1.0, 1.0],
            (x, y) -> (x[1] - 5)^2 + (2*y[1] + 1)^2,
            (x, y) -> (y[1] - 1)^2 - 1.5*x[1]*y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -3*x[1] + y[1] + 3,
                x[1] - 0.5*y[1] - 4,
                x[1] + y[1] - 7,
                -y[1]
            ]
        )

    # TODO corriger ce probleme
    elseif prob_no == 5 || prob_no == "Bard1988Ex2"
        return BilevelProblem(
            "Bard1988Ex2",
            [4, 4, 9, 12],
            [5.0, 5.0, 15.0, 15.0, 0.0, 0.0, 0.0, 0.0],
            [-6600.0, 54.0, 1.0],
            (x, y) -> -100*x[1] - 100*x[2] - x[3] - x[4],
            (x, y) -> y[1]^2 + y[2]^2 + y[3]^2 + y[4]^2,
            (x, y) -> [
                x[1] + x[2] + x[3] + x[4] - 40,
                x[1] + x[2] - 20,
                x[3] + x[4] - 20,
                -x[1],
                -x[2],
                -x[3],
                -x[4],
                x[1] - 15,
                x[2] - 15
            ],
            (x, y) -> [
                y[1] + y[2] + y[3] + y[4] - 40,
                y[1] + y[2] - 20,
                y[3] + y[4] - 20,
                -y[1],
                -y[2],
                -y[3],
                -y[4],
                y[1] - 15,
                y[2] - 15,
                y[1] - x[1],
                y[2] - x[2],
                y[3] - x[3]
            ],
            −6600.00
        )
    elseif prob_no == 6 || prob_no == "Bard1988Ex3"
        return BilevelProblem(
            "Bard1988Ex3",
            [2, 2, 3, 4],
            [0.0, 2.0, 4.0, 1.0],
            [-12.68, -1.02, 1.0],
            # Upper-level objective
            (x, y) -> -x[1]^2 - 3*x[2] - 4*y[1] + y[2]^2,
            # Lower-level objective
            (x, y) -> 2*x[1]^2 + y[1]^2 - 5*y[2],
            # Upper-level constraints
            (x, y) -> [
                x[1]^2 + 2*x[2] - 4,
                -x[1],
                -x[2]
            ],
            # Lower-level constraints
            (x, y) -> [
                -x[1]^2 + 2*x[1] - x[2]^2 + 2*y[1] - y[2] - 3,
                -x[2] - 3*y[1] + 4*y[2] + 4,
                -y[1],
                -y[2]
            ],
            -12.68  # Optimal value
        )
    elseif prob_no == 7 || prob_no == "Bard1991Ex1"
        return BilevelProblem(
            "Bard1991Ex1",
            [1, 2, 2, 3],
            [1.0, 1.0, 1.0],
            [2.0, 12.0, 1.0],
            (x, y) -> x[1] + y[2],
            (x, y) -> 2*y[1] + x[1]*y[2],
            (x, y) -> [
                -x[1] + 2,
                x[1] - 4
            ],
            (x, y) -> [
                x[1] - y[1] - y[2] + 4,
                -y[1],
                -y[2]
            ],
            [2.0, 6.0, 0.0]
        )
    elseif prob_no == 8 || prob_no == "BardBook1998"
        return BilevelProblem(
            "BardBook1998",
            [2, 2, 4, 7],
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 5.0, 1.0],
            (x, y) -> sum((y .- x .+ 20).^2),
            (x, y) -> 2*x[1] + 2*x[2] - 3*y[1] - 3*y[2] - 60,
            (x, y) -> [
                x[1] - 50,
                x[2] - 50,
                -x[1],
                -x[2]
            ],
            (x, y) -> [
                x[1] + x[2] + y[1] - 2*y[2] - 40,
                2*y[1] - x[1] + 10,
                y[1] - 20,
                -y[1] - 10,
                2*y[2] - x[2] + 10,
                y[2] - 20,
                -y[2] - 10
            ],
            [25.0, 30.0, 5.0, 10.0]
        )
    elseif prob_no == 9 || prob_no == "CalamaiVicente1994a"
        return BilevelProblem(
            "CalamaiVicente1994a",
            [1, 1, 0, 3],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> 0.5*(x[1] - 1)^2 + 0.5*y[1]^2,
            (x, y) -> 0.5*y[1] - x[1]*y[1],
            (x, y) -> Float64[],
            (x, y) -> [
                x[1] - y[1] - 1,
                -x[1] - y[1] + 1,
                x[1] + y[1] - 1
            ],
            [1.0, 0.0] # Disclaimer: this problem should depend on a ρ parameter not included here but not included in BOLIB in Matlab anyway.
        )
    elseif prob_no == 10 || prob_no == "CalamaiVicente1994b"
        return BilevelProblem(
            "CalamaiVicente1994b",
            [4, 2, 0, 6],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.3125, -0.4063, 1.0],
            (x, y) -> 0.5*sum((x .- 1).^2) + 0.5*sum(y.^2),
            (x, y) -> 0.5*sum(y.^2) - x[1]*y[1] - x[2]*y[2],
            (x, y) -> Float64[],
            (x, y) -> [
                x[1] - y[1] - 1,
                x[2] - y[2] - 1,
                -x[1] - y[1] - 1,
                -x[2] - y[2] - 1,
                x[1] + y[1] - 1.5,
                x[1] + y[2] - 3
            ],
            [1.25, 0.5, 1, 1, 0.25, 0.5]
        )
    elseif prob_no == 11 || prob_no == "CalamaiVicente1994c"
        # Data matrices as in the Matlab function
        A = [197.2  32.4  -129.6 -43.2;
            32.4   110.8 -43.2  -14.4;
            -129.6 -43.2   302.8 -32.4;
            -43.2  -14.4  -32.4   389.2]
        B = [100.0 0.0; 0.0 100.0]
        a = [-8.56; -9.52; -9.92; -16.64]
        C = [-132.4 -10.8;
            -10.8  -103.6;
            43.2   14.4;
            14.4   4.8]
        D = [13.24  1.08  -4.32 -1.44;
            1.08   10.36 -1.44 -0.48;
            13.24  1.08  -4.32 -1.44;
            1.08   10.36 -1.44 -0.48;
            -13.24 -1.08   4.32  1.44;
            -1.08  -10.36  1.44  0.48]
        E = [-10.0 0.0; 0.0 -10.0; 10.0 0.0; 0.0 10.0; -10.0 0.0; 0.0 -10.0]
        d = [-1.0; -1.0; -1.5; -3.0; 1.0; 1.0]
        return BilevelProblem(
            "CalamaiVicente1994c",
            [4, 2, 0, 6],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.3125, -0.4063, 1.0],
            (x, y) -> 0.5*(x' * A * x + y' * B * y) + dot(a, x) + 2,
            (x, y) -> x' * C * y + 0.5 * y' * B * y,
            (x, y) -> Float64[],
            (x, y) -> D * x + E * y + d,
            [0.13085, 0.05195, 0.1022, 0.0674, 0.025, 0.05]
        )
    # ...existing code...

    elseif prob_no == 12 || prob_no == "CalveteGale1999P1"
        return BilevelProblem(
            "CalveteGale1999P1",
            [2, 3, 2, 6],
            [0.0, 0.5, 0.0, 0.5, 0.0],
            [-29.2, 0.31, 1.0],
            (x, y) -> -8*x[1] - 4*x[2] + y[1] - 40*y[2] - 4*y[3],
            (x, y) -> (1 + x[1] + x[2] + 2*y[1] - y[2] + y[3]) / (6 + 2*x[1] + y[1] + y[2] - 3*y[3]),
            (x, y) -> [
                -x[1],
                -x[2]
            ],
            (x, y) -> [
                -y[1],
                -y[2],
                -y[3],
                -y[1] + y[2] + y[3] - 1,
                2*x[1] - y[1] + 2*y[2] - 0.5*y[3] - 1,
                2*x[2] + 2*y[1] - y[2] - 0.5*y[3] - 1
            ],
            [0.0, 0.9, 0.0, 0.6, 0.4]
        )

    elseif prob_no == 13 || prob_no == "ClarkWesterberg1990a"
        return BilevelProblem(
            "ClarkWesterberg1990a",
            [1, 1, 2, 3],
            [1.0, 1.0],
            [5.0, 4.0, 1.0],
            (x, y) -> (x[1] - 3)^2 + (y[1] - 2)^2,
            (x, y) -> (y[1] - 5)^2,
            (x, y) -> [
                x[1] - 8,
                -x[1]
            ],
            (x, y) -> [
                -2*x[1] + y[1] - 1,
                x[1] - 2*y[1] + 2,
                x[1] + 2*y[1] - 14
            ],
            [1.0, 3.0]
        )

    elseif prob_no == 14 || prob_no == "Colson2002BIPA1"
        return BilevelProblem(
            "Colson2002BIPA1",
            [1, 1, 3, 3],
            [10.0, 10.0],
            [250.0, 0.0, 1.0],
            (x, y) -> (10 - x[1])^3 + (10 - y[1])^3,
            (x, y) -> (x[1] + 2*y[1] - 15)^4,
            (x, y) -> [
                x[1] - 5,
                -x[1] + y[1],
                -x[1]
            ],
            (x, y) -> [
                x[1] + y[1] - 20,
                y[1] - 20,
                -y[1]
            ],
            [5.0, 5.0]
        )

    elseif prob_no == 15 || prob_no == "Colson2002BIPA2"
        return BilevelProblem(
            "Colson2002BIPA2",
            [1, 1, 1, 4],
            [4.0, 0.0],
            [17.0, 2.0, 2.0],
            (x, y) -> (x[1] - 5)^2 + (2*y[1] + 1)^2,
            (x, y) -> (y[1] - 1)^2 - 1.5*x[1]*y[1] + x[1]^3,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -3*x[1] + y[1] + 3,
                x[1] - 0.5*y[1] - 4,
                x[1] + y[1] - 7,
                -y[1]
            ],
            [1.0, 0.0]
        )

    elseif prob_no == 16 || prob_no == "Colson2002BIPA3"
        return BilevelProblem(
            "Colson2002BIPA3",
            [1, 1, 2, 2],
            [4.0, 0.0],
            [2.0, 24.02, 2.0],
            (x, y) -> (x[1] - 5)^4 + (2*y[1] + 1)^4,
            (x, y) -> exp(-x[1] + y[1]) + x[1]^2 + 2*x[1]*y[1] + y[1]^2 + 2*x[1] + 6*y[1],
            (x, y) -> [
                x[1] + y[1] - 4,
                -x[1]
            ],
            (x, y) -> [
                -x[1] + y[1] - 2,
                -y[1]
            ],
            [4.0, 0.0]
        )

    elseif prob_no == 17 || prob_no == "Colson2002BIPA4"
        return BilevelProblem(
            "Colson2002BIPA4",
            [1, 1, 2, 2],
            [1.5, 2.25],
            [88.79, -0.77, 2.0],
            (x, y) -> x[1]^2 + (y[1] - 10)^2,
            (x, y) -> x[1]^3 + 2*y[1]^3 + x[1] - 2*y[1] - x[1]^2,
            (x, y) -> [
                x[1] + 2*y[1] - 6,
                -x[1]
            ],
            (x, y) -> [
                -x[1] + 2*y[1] - 3,
                -y[1]
            ],
            [0.0, 0.6039]
        )

    elseif prob_no == 18 || prob_no == "Colson2002BIPA5"
    return BilevelProblem(
        "Colson2002BIPA5",
        [2, 2, 2, 2],
        [1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (x[2] - 1)^2 + (y[1] - 1)^2 + (y[2] - 1)^2,
        (x, y) -> (y[1] - x[1])^2 + (y[2] - x[2])^2,
        (x, y) -> [
            x[1] + x[2] + y[1] + y[2] - 6,
            -x[1]
        ],
        (x, y) -> [
            y[1] + y[2] - 3,
            -y[1]
        ]
    )
elseif prob_no == 19 || prob_no == "Dempe1992a"
    return BilevelProblem(
        "Dempe1992a",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 20 || prob_no == "Dempe1992b"
    return BilevelProblem(
        "Dempe1992b",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 21 || prob_no == "DeSilva1978"
    return BilevelProblem(
        "DeSilva1978",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 22 || prob_no == "FalkLiu1995"
    return BilevelProblem(
        "FalkLiu1995",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 23 || prob_no == "Vogel2012"
    return BilevelProblem(
        "Vogel2012",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 24 || prob_no == "Zlobec2001a"
    return BilevelProblem(
        "Zlobec2001a",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 25 || prob_no == "Zlobec2001b"
    return BilevelProblem(
        "Zlobec2001b",
        [1, 1, 1, 2],
        [1.0, 1.0],
        [0.0, 0.0, 1.0],
        (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
        (x, y) -> (y[1] - x[1])^2,
        (x, y) -> [-x[1]],
        (x, y) -> [
            y[1] - 1,
            -y[1]
        ]
    )
elseif prob_no == 26 || prob_no == "AnEtal2009"
    return BilevelProblem(
        "AnEtal2009",
        [2, 2, 2, 2],
        [1.0, 1.0, 1.0, 1.0],
        [2251.55, 565.78, 1.0],
        (x, y) -> (x[1] - 1)^2 + (x[2] - 1)^2 + (y[1] - 1)^2 + (y[2] - 1)^2,
        (x, y) -> (y[1] - x[1])^2 + (y[2] - x[2])^2,
        (x, y) -> [
            x[1] + x[2] + y[1] + y[2] - 6,
            -x[1]
        ],
        (x, y) -> [
            y[1] + y[2] - 3,
            -y[1]
        ]
    )
elseif prob_no == 27 || prob_no == "AllendeStill2013"
    return BilevelProblem(
        "AllendeStill2013",
        [2, 2, 5, 2],
        [0.0, 0.0, 0.0, 0.0],
        [1.0, -0.5, 1.0],
        (x, y) -> (x[1]-1)^2 + (x[2]-1)^2 + y[1]^2 + y[2]^2,
        (x, y) -> y[1]^2 + y[2]^2 - 2*x[1]*y[1] - 2*x[2]*y[2],
        (x, y) -> [
            -x[1],
            -x[2],
            x[1]-2,
            -y[1],
            -y[2]
        ],
        (x, y) -> [
            (y[1]-1)^2-0.25,
            (y[2]-1)^2-0.25
        ]
    )

    elseif prob_no == 25 || prob_no == "DempeFranke2011Ex42"
        return BilevelProblem(
            "DempeFranke2011Ex42",
            [2, 2, 4, 3],
            [1.0, 1.0, 1.0, 1.0],
            [2.13, -3.5, 1.0],
            (x, y) -> x[1] + sum((y .- [1.0, 0.0]).^2),
            (x, y) -> x[1]*y[1] + x[2]*y[2],
            (x, y) -> [
                -1 - x[1],
                -1 + x[1],
                -1 - x[2],
                1 + x[2]
            ],
            (x, y) -> [
                -y[1] + y[2] - 1,
                y[1] + y[2] - 3.5,
                y[2] - 2
            ],
            [1.0, −1.0, 0.0, 1.0]
        )

    elseif prob_no == 26 || prob_no == "DempeFranke2014Ex38"
        return BilevelProblem(
            "DempeFranke2014Ex38",
            [2, 2, 4, 4],
            [1.0, 1.0, 1.0, 1.0],
            [-1.0, -4.0, 1.0],
            (x, y) -> 2*x[1] + x[2] + 2*y[1] - y[2],
            (x, y) -> x[1]*y[1] + x[2]*y[2],
            (x, y) -> [
                -1 - x[1],
                -1 + x[1],
                -1 - x[2],
                0.75 + x[2]
            ],
            (x, y) -> [
                -2*y[1] + y[2],
                y[1] - 2,
                y[2] - 2,
                -y[2]
            ],
            [−1.0, −1.0, 2.0, 2.0]
        )

    elseif prob_no == 27 || prob_no == "DempeLohse2011Ex31a"
        return BilevelProblem(
            "DempeLohse2011Ex31a",
            [2, 2, 0, 4],
            [1.0, 1.0, 1.0, 1.0],
            [-5.5, 0.0, 1.0],
            (x, y) -> sum((x .- 0.5).^2) - 3*y[1] - 3*y[2],
            (x, y) -> x[1]*y[1] + x[2]*y[2],
            (x, y) -> Float64[],
            (x, y) -> [
                sum(y) - 2,
                -y[1] + y[2],
                -y[1],
                -y[2]
            ],
            [0.0, 0.0, 1.0, 1.0]
        )

    elseif prob_no == 28 || prob_no == "DempeLohse2011Ex31b"
        return BilevelProblem(
            "DempeLohse2011Ex31b",
            [3, 3, 0, 5],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-12.0, 0.0, 1.0],
            (x, y) -> sum((x .- [0.5, 0.5, 0.0]).^2) - 3*y[1] - 3*y[2] - 6*y[3],
            (x, y) -> sum(x .* y),
            (x, y) -> Float64[],
            (x, y) -> [
                sum(y) - 2,
                -y[1] + y[2],
                -y[1],
                -y[2],
                -y[3]
            ],
            [0.5, 0.5, 0.0, 1.0, 1.0, 0.0]
        )

    elseif prob_no == 29 || prob_no == "DeSilva1978"
        return BilevelProblem(
            "DeSilva1978",
            [2, 2, 0, 4],
            [0.0, 0.0, 1.0, 1.0],
            [-1.0, 0.0, 1.0],
            (x, y) -> sum((x .- 1.0).^2) + sum(y.^2) - 2,
            (x, y) -> sum((y .- x).^2),
            (x, y) -> Float64[],
            (x, y) -> [
                -y[1] + 0.5,
                -y[2] + 0.5,
                y[1] - 1.5,
                y[2] - 1.5
            ],
            [0.5, 0.5, 0.5, 0.5]
        )

    elseif prob_no == 30 || prob_no == "FalkLiu1995"
        return BilevelProblem(
            "FalkLiu1995",
            [2, 2, 0, 4],
            [0.0, 0.0, 1.0, 1.0],
            [-2.1962, 0.0, 1.0],
            (x, y) -> sum((x .- 1.5).^2) + sum(y.^2) - 4.5,
            (x, y) -> sum((y .- x).^2),
            (x, y) -> Float64[],
            (x, y) -> [
                -y[1] + 0.5,
                -y[2] + 0.5,
                y[1] - 1.5,
                y[2] - 1.5
            ],
            [sqrt(3)/2, sqrt(3)/2, sqrt(3)/2, sqrt(3)/2]
        )

    elseif prob_no == 31 || prob_no == "FloudasEtal2013"
        return BilevelProblem(
            "FloudasEtal2013",
            [2, 2, 4, 7],
            [10.0, 10.0, 20.0, 20.0],
            [0.0, 200.0, 1.0],
            (x, y) -> 2*x[1] + 2*x[2] - 3*y[1] - 3*y[2] - 60,
            (x, y) -> sum((y .- x .+ 20).^2),
            (x, y) -> [
                x[1] - 50,
                x[2] - 50,
                -x[1],
                -x[2]
            ],
            (x, y) -> [
                2*y[1] - x[1] + 10,
                x[1] + x[2] + y[1] - 2*y[2] - 40,
                -y[1] - 10,
                y[1] - 20,
                2*y[2] - x[2] + 10,
                -y[2] - 10,
                y[2] - 20
            ],
            [0.0, 0.0, −10.0, −10.0]
        )

    elseif prob_no == 32 || prob_no == "FloudasZlobec1998"
        return BilevelProblem(
            "FloudasZlobec1998",
            [1, 2, 2, 6],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            (x, y) -> x[1]^3 * y[1] + y[2],
            (x, y) -> -y[2],
            (x, y) -> [
                -x[1],
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                -y[2],
                y[1] - 1,
                y[2] - 100,
                x[1]*y[1] - 10,
                y[1]^2 + x[1]*y[2] - 1
            ],
            [1.0, 0.0, 1.0]
        )

    elseif prob_no == 33 || prob_no == "GumusFloudas2001Ex1"
        return BilevelProblem(
            "GumusFloudas2001Ex1",
            [1, 1, 3, 3],
            [1.0, 1.0],
            [2250.0, 197.75, 1.0],
            (x, y) -> 16*x[1]^2 + 9*y[1]^2,
            (x, y) -> (x[1] + y[1] - 20)^4,
            (x, y) -> [
                -x[1],
                x[1] - 12.5,
                -4*x[1] + y[1]
            ],
            (x, y) -> [
                -y[1],
                y[1] - 50,
                4*x[1] + y[1] - 50
            ],
            [11.25, 5.0]
        )

    elseif prob_no == 34 || prob_no == "GumusFloudas2001Ex3"
        return BilevelProblem(
            "GumusFloudas2001Ex3",
            [2, 3, 4, 9],
            [0.0, 0.5, 0.0, 0.5, 0.0],
            [-29.2, 0.31, 1.0],
            (x, y) -> -8*x[1] - 4*x[2] + y[1] - 40*y[2] - 4*y[3],
            (x, y) -> (1 + x[1] + x[2] + 2*y[1] - y[2] + y[3]) / (6 + 2*x[1] + y[1] + y[2] - 3*y[3]),
            (x, y) -> [
                -x[1],
                -x[2],
                x[1] - 2,
                x[2] - 2
            ],
            (x, y) -> [
                -y[1],
                -y[2],
                -y[3],
                y[1] - 2,
                y[2] - 2,
                y[3] - 2,
                -y[1] + y[2] + y[3] - 1,
                2*x[1] - y[1] + 2*y[2] - 0.5*y[3] - 1,
                2*x[2] + 2*y[1] - y[2] - 0.5*y[3] - 1
            ],
            [0.0, 0.9, 0.0, 0.6, 0.4]
        )

    elseif prob_no == 35 || prob_no == "GumusFloudas2001Ex4"
        return BilevelProblem(
            "GumusFloudas2001Ex4",
            [1, 1, 5, 2],
            [1.0, 1.0],
            [9.0, 0.0, 1.0],
            (x, y) -> (x[1] - 3)^2 + (y[1] - 2)^2,
            (x, y) -> (y[1] - 5)^2,
            (x, y) -> [
                -x[1],
                x[1] - 8,
                -2*x[1] + y[1] - 1,
                x[1] - 2*y[1] + 2,
                x[1] + 2*y[1] - 14
            ],
            (x, y) -> [
                -y[1],
                y[1] - 10
            ],
            [3.0, 5.0]
        )

    elseif prob_no == 36 || prob_no == "GumusFloudas2001Ex5"
        return BilevelProblem(
            "GumusFloudas2001Ex5",
            [1, 2, 2, 6],
            [1.0, 1.0, 1.0],
            [0.194, -7.23, 1.0],
            (x, y) -> x[1],
            (x, y) -> -y[1] + 0.5864*y[1]^0.67,
            (x, y) -> [
                -x[1] + 0.1,
                x[1] - 10
            ],
            (x, y) -> [
                -y[1] + 0.1,
                -y[2] + 0.1,
                y[1] - 10,
                y[2] - 10,
                0.0332333/y[2] + 0.1*y[1] - 1,
                (4*x[1] + 2*x[1]^(-0.71))/y[2] + 0.0332333*x[1]^(-1.3) - 1
            ],
            [0.193616, 9.9667667, 10.0]
        )

    elseif prob_no == 37 || prob_no == "HatzEtal2013"
        return BilevelProblem(
            "HatzEtal2013",
            [1, 2, 0, 2],
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> -x[1] + 2*y[1] + y[2],
            (x, y) -> (x[1] - y[1])^2 + y[2]^2,
            (x, y) -> Float64[],
            (x, y) -> [
                -y[1],
                -y[2]
            ],
            [0.0, 0.0, 0.0]
        )

    elseif prob_no == 38 || prob_no == "HendersonQuandt1958"
        return BilevelProblem(
            "HendersonQuandt1958",
            [1, 1, 2, 1],
            [0.0, 0.0],
            [-3266.67, -711.11, 2.0],
            (x, y) -> (0.5*(x[1] + y[1]) - 95)*x[1],
            (x, y) -> (y[1] + 0.5*x[1] - 100)*y[1],
            (x, y) -> [
                x[1] - 200,
                -x[1]
            ],
            (x, y) -> [-y[1]],
            [93.33333, 26.667]
        )

    elseif prob_no == 39 || prob_no == "HenrionSurowiec2011"
        return BilevelProblem(
            "HenrionSurowiec2011",
            [1, 1, 0, 0],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> x[1]^2,
            (x, y) -> (y[1]/2 - x[1])*y[1],
            (x, y) -> Float64[],
            (x, y) -> Float64[],
            [0.0, 0.0]
        )

    elseif prob_no == 40 || prob_no == "IshizukaAiyoshi1992a"
        return BilevelProblem(
            "IshizukaAiyoshi1992a",
            [1, 2, 1, 5],
            [1.0, 1.0, 1.0],
            [0.0, -1.5, 1.0],
            (x, y) -> x[1]*y[2]^2,
            (x, y) -> y[1],
            (x, y) -> [-x[1] - 1.5],
            (x, y) -> [
                -x[1] - y[1] - y[2],
                -x[1] - y[1] + y[2],
                -x[1] + y[1] - y[2],
                -y[1] - y[2] - 1.5,
                -y[1] + y[2] - 1.5
            ],
            []
        )

    elseif prob_no == 41 || prob_no == "KleniatiAdjiman2014Ex3"
        return BilevelProblem(
            "KleniatiAdjiman2014Ex3",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1.0, 0.0, 1.0],
            (x, y) -> x[1] - y[1],
            (x, y) -> x[1]*y[1]^2/2 - x[1]*y[1]^3,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.0, 1.0]
        )

    elseif prob_no == 42 || prob_no == "KleniatiAdjiman2014Ex4"
        return BilevelProblem(
            "KleniatiAdjiman2014Ex4",
            [5, 5, 13, 11],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-10.0, -3.1, 2.0],
            (x, y) -> -sum(x.^2 + y.^2),
            (x, y) -> y[1]^3 + (x[1] + x[2])*y[2]^2 + 0.1*y[3] + (y[4]^2 + y[5]^2)*x[3]*x[4]*x[5],
            (x, y) -> [
                -x[1] - 1,
                -x[2] - 1,
                -x[3] - 1,
                -x[4] - 1,
                -x[5] - 1,
                x[1] - 1,
                x[2] - 1,
                x[3] - 1,
                x[4] - 1,
                x[5] - 1,
                y[1]*y[2] - x[1],
                x[1] - exp(x[2]) + y[3],
                x[2]*y[1]^2
            ],
            (x, y) -> [
                -y[1] - 1,
                -y[2] - 1,
                -y[3] - 1,
                -y[4] - 1,
                -y[5] - 1,
                y[1] - 1,
                y[2] - 1,
                y[3] - 1,
                y[4] - 1,
                y[5] - 1,
                x[1] - 0.2 - y[3]^2
            ],
            [1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
        )

    elseif prob_no == 43 || prob_no == "LamparielloSagratella2017Ex23"
        return BilevelProblem(
            "LamparielloSagratella2017Ex23",
            [1, 2, 2, 2],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
            (x, y) -> x[1],
            (x, y) -> (x[1] - y[1])^2 + (y[2] + 1)^2,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                y[1]^3 - y[2],
                -y[2]
            ],
            [-1.0, -1.0, 0.0]
        )

    elseif prob_no == 44 || prob_no == "LamparielloSagratella2017Ex31"
        return BilevelProblem(
            "LamparielloSagratella2017Ex31",
            [1, 1, 1, 1],
            [1.0, 1.0],
            [1.0, 0.0, 1.0],
            (x, y) -> x[1]^2 + y[1]^2,
            (x, y) -> y[1],
            (x, y) -> [1 - x[1]],
            (x, y) -> [1 - x[1] - y[1]],
            [1.0, 0.0]
        )

    elseif prob_no == 45 || prob_no == "LamparielloSagratella2017Ex32"
        return BilevelProblem(
            "LamparielloSagratella2017Ex32",
            [1, 1, 0, 0],
            [1.0, 1.0],
            [0.5, 0.0, 1.0],
            (x, y) -> x[1]^2 + y[1]^2,
            (x, y) -> (x[1] + y[1] - 1)^2,
            (x, y) -> Float64[],
            (x, y) -> Float64[],
            [0.5, 0.5]
        )

    elseif prob_no == 46 || prob_no == "LamparielloSagratella2017Ex33"
        return BilevelProblem(
            "LamparielloSagratella2017Ex33",
            [1, 2, 1, 3],
            [1.0, 1.0, 1.0],
            [0.5, 0.0, 1.0],
            (x, y) -> x[1]^2 + (y[1] + y[2])^2,
            (x, y) -> y[1],
            (x, y) -> [0.5 - x[1]],
            (x, y) -> [
                1 - x[1] - y[1] - y[2],
                -y[1],
                -y[2]
            ],
            [0.5, 0.0, 0.5]
        )

    elseif prob_no == 47 || prob_no == "LamparielloSagratella2017Ex35"
        return BilevelProblem(
            "LamparielloSagratella2017Ex35",
            [1, 1, 2, 3],
            [1.0, 1.0],
            [0.8, -0.4, 1.0],
            (x, y) -> x[1]^2 + y[1]^2,
            (x, y) -> -y[1],
            (x, y) -> [
                -1 - x[1],
                x[1] - 1
            ],
            (x, y) -> [
                2*x[1] + y[1] - 2,
                -y[1],
                y[1] - 1
            ],
            [4/5, 2/5]
        )

    elseif prob_no == 48 || prob_no == "LucchettiEtal1987"
        return BilevelProblem(
            "LucchettiEtal1987",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> (1 - x[1])/2 + x[1]*y[1],
            (x, y) -> (x[1] - 1)*y[1],
            (x, y) -> [
                -x[1],
                x[1] - 1
            ],
            (x, y) -> [
                -y[1],
                y[1] - 1
            ],
            [1.0, 0.0]
        )

    elseif prob_no == 49 || prob_no == "LuDebSinha2016a"
        return BilevelProblem(
            "LuDebSinha2016a",
            [1, 1, 4, 0],
            [0.0, 2.0],
            [1.1360, 1.1838, 2.0],
            (x, y) -> begin
                a = (0.2*y[1] - x[1] + 0.6)/0.055
                2 - exp(-a^0.4) - 0.8*exp(-((0.15*y[1] + x[1] - 0.4)/0.3)^2)
            end,
            (x, y) -> begin
                a = (1.5*y[1] - x[1])/0.055
                2 - exp(-a^0.4) - 0.8*exp(-((2*y[1] + x[1] - 3)/0.5)^2)
            end,
            (x, y) -> [
                -x[1],
                x[1] - 1,
                -y[1],
                y[1] - 2
            ],
            (x, y) -> Float64[],
            [0.2, 1.4]
        )

    elseif prob_no == 50 || prob_no == "LuDebSinha2016b"
        return BilevelProblem(
            "LuDebSinha2016b",
            [1, 1, 4, 0],
            [1.0, 1.0],
            [0.0, 1.6645, 2.0],
            (x, y) -> (x[1] - 0.5)^2 + (y[1] - 1)^2,
            (x, y) -> begin
                a = (1.5*y[1] - x[1])/0.055
                2 - exp(-a^0.4) - 0.8*exp(-((2*y[1] + x[1] - 3)/0.5)^2)
            end,
            (x, y) -> [
                -x[1],
                x[1] - 1,
                -y[1],
                y[1] - 2
            ],
            (x, y) -> Float64[],
            [0.5, 1.0]
        )

    elseif prob_no == 51 || prob_no == "LuDebSinha2016c"
        return BilevelProblem(
            "LuDebSinha2016c",
            [1, 1, 4, 0],
            [1.0, 1.0],
            [1.12, 0.06, 2.0],
            (x, y) -> begin
                a = (0.2*y[1] - x[1] + 0.6)/0.055
                2 - exp(-a^0.4) - 0.8*exp(-((0.15*y[1] + x[1] - 0.4)/0.3)^2)
            end,
            (x, y) -> (x[1] - 0.5)^2 + (y[1] - 1)^2,
            (x, y) -> [
                -x[1],
                x[1] - 1,
                -y[1],
                y[1] - 2
            ],
            (x, y) -> Float64[],
            [0.26, 1.0]
        )

    elseif prob_no == 52 || prob_no == "LuDebSinha2016d"
        return BilevelProblem(
            "LuDebSinha2016d",
            [2, 2, 11, 3],
            [1.0, 1.0, 1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> -x[2],
            (x, y) -> -y[2],
            (x, y) -> [
                -4 - x[1],
                -100 - x[2],
                x[1] - 10,
                x[2] - 200,
                -4 - y[1],
                -100 - y[2],
                y[1] - 10,
                y[2] - 200,
                -(y[1]/14 + 16/7)*(x[1] - 2)^2 + x[2],
                12.5*(y[1]/14 + 16/7)*(x[1] - 5) - x[2],
                -5*(x[1] + 4 - (y[1]/14 + 16/7))*(x[1] + 8 - (y[1]/14 + 16/7)) + x[2]
            ],
            (x, y) -> [
                -(x[1]/14 + 16/7)*(y[1] - 2)^2 + y[2],
                12.5*(x[1]/14 + 16/7)*(y[1] - 5) - y[2],
                -5*(y[1] + 4 - (x[1]/14 + 16/7))*(y[1] + 8 - (x[1]/14 + 16/7)) + y[2]
            ],
            [10.0, 192.0, 10.0, 192.0]
        )

    elseif prob_no == 53 || prob_no == "LuDebSinha2016e"
        return BilevelProblem(
            "LuDebSinha2016e",
            [1, 2, 6, 3],
            [1.0, 1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> (x[1] - 2.5)^2 / 0.04 + (y[2] - 50)^2 / 900,
            (x, y) -> -y[2],
            (x, y) -> [
                2 - x[1],
                x[1] - 3,
                -4 - y[1],
                -100 - y[2],
                y[1] - 10,
                y[2] - 200
            ],
            (x, y) -> [
                -x[1]*(y[1] - 2)^2 + y[2],
                12.5*x[1]*(y[1] - 5) - y[2],
                -5*(y[1] + 4 - x[1])*(y[1] + 8 - x[1]) + y[2]
            ],
            []
        )

    elseif prob_no == 54 || prob_no == "LuDebSinha2016f"
        return BilevelProblem(
            "LuDebSinha2016f",
            [2, 1, 9, 0],
            [1.0, 1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> -x[2],
            (x, y) -> (x[1] - 50)^2 / 784 + (y[1] - 2.5)^2 / 0.04,
            (x, y) -> [
                2 - y[1],
                y[1] - 4,
                -80 - x[1],
                -100 - x[2],
                x[1] - 200,
                x[2] - 200,
                -y[1]*(x[1]/20 - 2)^2 + x[2],
                12.5*y[1]*(x[1]/20 - 5) - x[2],
                -5*(x[1]/20 + 4 - y[1])*(x[1]/20 + 8 - y[1]) + x[2]
            ],
            (x, y) -> Float64[],
            []
        )

    elseif prob_no == 55 || prob_no == "MacalHurter1997"
        return BilevelProblem(
            "MacalHurter1997",
            [1, 1, 0, 0],
            [1.0, 1.0],
            [81.327, -0.333, 1.0],
            (x, y) -> (x[1] - 1)^2 + (y[1] - 1)^2,
            (x, y) -> 0.5*y[1]^2 + 500*y[1] - 50*x[1]*y[1],
            (x, y) -> Float64[],
            (x, y) -> Float64[],
            [10.0163, 0.8197]
        )

    elseif prob_no == 56 || prob_no == "Mirrlees1999"
        return BilevelProblem(
            "Mirrlees1999",
            [1, 1, 0, 2],
            [1.0, 1.0],
            [1.002, -1.02, 1.0],
            (x, y) -> (x[1] - 2)^2 + (y[1] - 1)^2,
            (x, y) -> -x[1]*exp(-(y[1] + 1)^2) - exp(-(y[1] - 1)^2),
            (x, y) -> Float64[],
            (x, y) -> [
                y[1] - 2,
                -y[1] - 2
            ],
            [1.0, 0.95753]
        )

    elseif prob_no == 57 || prob_no == "MitsosBarton2006Ex38"
        return BilevelProblem(
            "MitsosBarton2006Ex38",
            [1, 1, 4, 2],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> y[1]^2,
            (x, y) -> x[1]*y[1] + exp(x[1])*y[1],
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1,
                -y[1] - 0.1,
                y[1] - 0.1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [−0.567, 0.0]
        )

    elseif prob_no == 58 || prob_no == "MitsosBarton2006Ex39"
        return BilevelProblem(
            "MitsosBarton2006Ex39",
            [1, 1, 3, 2],
            [1.0, 1.0],
            [-1.0, -1.0, 1.0],
            (x, y) -> x[1],
            (x, y) -> y[1]^3,
            (x, y) -> [
                - x[1] + y[1],
                - x[1] - 10,
                x[1] - 10
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
        [−1.0, −1.0]
        )

    elseif prob_no == 59 || prob_no == "MitsosBarton2006Ex310"
        return BilevelProblem(
            "MitsosBarton2006Ex310",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.5, -0.1, 1.0],
            (x, y) -> y[1],
            (x, y) -> x[1]*(16*y[1]^4 + 2*y[1]^3 - 8*y[1]^2 - 1.5*y[1] + 0.5),
            (x, y) -> [
                -x[1] + 0.1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [1.0, 0.5] # Careful, not the only global optimum
        )

    elseif prob_no == 60 || prob_no == "MitsosBarton2006Ex311"
        return BilevelProblem(
            "MitsosBarton2006Ex311",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-0.8, 0.0, 1.0],
            (x, y) -> y[1],
            (x, y) -> x[1]*(16*y[1]^4 + 2*y[1]^3 - 8*y[1]^2 - 1.5*y[1] + 0.5),
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 0.8,
                y[1] - 1
            ],
            [0.0, −0.8]
        )

    elseif prob_no == 61 || prob_no == "MitsosBarton2006Ex312"
        return BilevelProblem(
            "MitsosBarton2006Ex312",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> -x[1] + x[1]*y[1] + 10*y[1]^2,
            (x, y) -> -x[1]*y[1]^2 + 0.5*y[1]^4,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.0, 0.0]
        )

    elseif prob_no == 62 || prob_no == "MitsosBarton2006Ex313"
        return BilevelProblem(
            "MitsosBarton2006Ex313",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1.0, 0.0, 1.0],
            (x, y) -> x[1] - y[1],
            (x, y) -> x[1]*y[1]*(y[1]/2 - x[1]^2),
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.0, 1.0]
        )

    elseif prob_no == 63 || prob_no == "MitsosBarton2006Ex314"
        return BilevelProblem(
            "MitsosBarton2006Ex314",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.25, -1/12, 1.0],
            (x, y) -> (x[1] - 0.25)^2 + y[1]^2,
            (x, y) -> y[1]^3 / 3 - x[1]*y[1],
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.25, 0.5]
        )

    elseif prob_no == 64 || prob_no == "MitsosBarton2006Ex315"
        return BilevelProblem(
            "MitsosBarton2006Ex315",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.0, -0.83, 1.0],
            (x, y) -> x[1] + y[1],
            (x, y) -> x[1]*y[1]^2 / 2 - y[1]^3 / 3,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
        [-1.0, 1.0]
        )

    elseif prob_no == 65 || prob_no == "MitsosBarton2006Ex316"
        return BilevelProblem(
            "MitsosBarton2006Ex316",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-2.0, 0.0, 1.0],
            (x, y) -> 2*x[1] + y[1],
            (x, y) -> -x[1]*y[1]^2 / 2 - y[1]^4 / 4,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [−1.0, 0.0] # [-0.5, -1.0] is also a global optimum
        )

    elseif prob_no == 66 || prob_no == "MitsosBarton2006Ex317"
        return BilevelProblem(
            "MitsosBarton2006Ex317",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [3/16, -1/64, 1.0],
            (x, y) -> (x[1] + 0.5)^2 + y[1]^2 / 2,
            (x, y) -> x[1]*y[1]^2 / 2 + y[1]^4 / 4,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [−0.25, 0.5] # [-0.25, -0.5] is also a global optimum
        )

    elseif prob_no == 67 || prob_no == "MitsosBarton2006Ex318"
        return BilevelProblem(
            "MitsosBarton2006Ex318",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1/4, 0.0, 1.0],
            (x, y) -> -x[1]^2 + y[1]^2,
            (x, y) -> x[1]*y[1]^2 - y[1]^4 / 2,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.5, 0.0]
        )

    elseif prob_no == 68 || prob_no == "MitsosBarton2006Ex319"
        return BilevelProblem(
            "MitsosBarton2006Ex319",
            [1, 1, 2, 2],
            [-1.0, 1.0],
            [-0.258, -0.0178, 1.0],
            (x, y) -> (x[1] - 1 + y[1]/2)*y[1],
            (x, y) -> (-x[1] + y[1]^2 / 2)*y[1]^2,
            (x, y) -> [
                -x[1] - 1,
                x[1] - 1
            ],
            (x, y) -> [
                -y[1] - 1,
                y[1] - 1
            ],
            [0.189, 0.4343]
        )

    # ...existing code...

    elseif prob_no == 69 || prob_no == "MitsosBarton2006Ex320"
        return BilevelProblem(
            "MitsosBarton2006Ex320",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.3125, -0.0833, 1.0],
            (x, y) -> (x[1] - 0.25)^2 + y[1]^2,
            (x, y) -> y[1]^3/3 - x[1]^2*y[1],
            (x, y) -> [-x[1] - 1, x[1] - 1],
            (x, y) -> [-y[1] - 1, y[1] - 1],
            [0.5, 0.5]
        )

    elseif prob_no == 70 || prob_no == "MitsosBarton2006Ex321"
        return BilevelProblem(
            "MitsosBarton2006Ex321",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [0.2095, -0.0656, 1.0],
            (x, y) -> (x[1] + 0.6)^2 + y[1]^2,
            (x, y) -> y[1]^4 + (4/30)*(1-x[1])*y[1]^3 + (-0.02*x[1]^2 + 0.16*x[1] - 0.4)*y[1]^2 + (0.004*x[1]^3 - 0.036*x[1]^2 + 0.08*x[1])*y[1],
            (x, y) -> [-x[1] - 1, x[1] - 1],
            (x, y) -> [-y[1] - 1, y[1] - 1],
            [−0.5545, 0.4554]
        )

    elseif prob_no == 71 || prob_no == "MitsosBarton2006Ex322"
        return BilevelProblem(
            "MitsosBarton2006Ex322",
            [1, 1, 2, 3],
            [1.0, 1.0],
            [0.2095, -0.0656, 1.0],
            (x, y) -> (x[1] + 0.6)^2 + y[1]^2,
            (x, y) -> y[1]^4 + (4/30)*(1-x[1])*y[1]^3 + (-0.02*x[1]^2 + 0.16*x[1] - 0.4)*y[1]^2 + (0.004*x[1]^3 - 0.036*x[1]^2 + 0.08*x[1])*y[1],
            (x, y) -> [-x[1] - 1, x[1] - 1],
            (x, y) -> [-y[1] - 1, y[1] - 1, 0.01*(1 + x[1]^2) - y[1]^2],
            [−0.5545, 0.4554]
        )

    elseif prob_no == 72 || prob_no == "MitsosBarton2006Ex323"
        return BilevelProblem(
            "MitsosBarton2006Ex323",
            [1, 1, 3, 3],
            [0.0, 1.0],
            [0.176, -1.0, 1.0],
            (x, y) -> x[1]^2,
            (x, y) -> y[1],
            (x, y) -> [-x[1] - 1, x[1] - 1, 1 + x[1] - 9*x[1]^2 - y[1]],
            (x, y) -> [-y[1] - 1, y[1] - 1, y[1]^2*(x[1] - 0.5)],
            [−0.4191, −1.0]
        )

    elseif prob_no == 73 || prob_no == "MitsosBarton2006Ex324"
        return BilevelProblem(
            "MitsosBarton2006Ex324",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1.755, 0.0, 1.0],
            (x, y) -> x[1]^2 - y[1],
            (x, y) -> ((y[1] - 1 - 0.1*x[1])^2 - 0.5 - 0.5*x[1])^2,
            (x, y) -> [-x[1], x[1] - 1],
            (x, y) -> [-y[1], y[1] - 3],
            [0.2106, 1.799]
        )

    elseif prob_no == 74 || prob_no == "MitsosBarton2006Ex325"
        return BilevelProblem(
            "MitsosBarton2006Ex325",
            [2, 3, 6, 9],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [-1.0, -2.0, 2.0],
            (x, y) -> x[1]*y[1] + x[2]*y[1]^2 - x[1]*x[2]*y[3],
            (x, y) -> x[1]*y[1]^2 + x[2]*y[2]*y[3],
            (x, y) -> [-x[1] - 1, x[1] - 1, 0.1*y[1]*y[2] - x[1]^2, x[2]*y[1]^2],
            (x, y) -> [
                -y[1] - 1, -y[2] - 1, -y[3] - 1,
                y[1] - 1, y[2] - 1, y[3] - 1,
                y[1]^2 - y[2]*y[3],
                y[2]^2*y[3] - y[1]*x[1],
                -y[3]^2 + 0.1
            ],
            [−1.0, −1.0, −1.0, 1.0, 1.0]
        )

    elseif prob_no == 75 || prob_no == "MitsosBarton2006Ex326"
        return BilevelProblem(
            "MitsosBarton2006Ex326",
            [2, 3, 7, 6],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [-2.354, -2.0, 1.0],
            (x, y) -> x[1]*y[1] + x[2]*y[2]^2 + x[1]*x[2]*y[3]^3,
            (x, y) -> x[1]*y[1]^2 + x[2]*y[2]^2 + (x[1] - x[2])*y[3]^2,
            (x, y) -> [
                -x[1] - 1, -x[2] - 1, x[1] - 1, x[2] - 1,
                0.1 - x[1]^2, 1.5 - sum(y.^2), -2.5 + sum(y.^2)
            ],
            (x, y) -> [
                -y[1] - 1, -y[2] - 1, -y[3] - 1,
                y[1] - 1, y[2] - 1, y[3] - 1
            ],
            [−1.0, −1.0, 1.0, 1.0, −0.707]
        )

    elseif prob_no == 76 || prob_no == "MitsosBarton2006Ex327"
        return BilevelProblem(
            "MitsosBarton2006Ex327",
            [5, 5, 13, 13],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [2.0, -1.1, 2.0],
            (x, y) -> sum(x.^2 + y.^2),
            (x, y) -> y[1]^3 + (x[1] + x[2])*y[2]^2 + 0.1*y[3] + (y[4]^2 + y[5]^2)*x[3]*x[4]*x[5],
            (x, y) -> [
                -x[1] - 1, -x[2] - 1, -x[3] - 1, -x[4] - 1, -x[5] - 1,
                x[1] - 1, x[2] - 1, x[3] - 1, x[4] - 1, x[5] - 1,
                y[1]*y[2] - x[1], x[1] - exp(x[2]) + y[3], x[2]*y[1]^2
            ],
            (x, y) -> [
                -y[1] - 1, -y[2] - 1, -y[3] - 1, -y[4] - 1, -y[5] - 1,
                y[1] - 1, y[2] - 1, y[3] - 1, y[4] - 1, y[5] - 1,
                y[1]*y[2] - 0.3, x[1] - 0.2 - y[3]^2, -exp(y[3]) + y[4]*y[5] - 0.1
            ],
            [0.0, 0.0, 0.0, 0.0, 0.0, −1.0, 0.0, −1.0, 0.0, 0.0]
        )

    elseif prob_no == 77 || prob_no == "MitsosBarton2006Ex328"
        return BilevelProblem(
            "MitsosBarton2006Ex328",
            [5, 5, 13, 13],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-10.0, -3.1, 2.0],
            (x, y) -> -sum(x.^2 + y.^2),
            (x, y) -> y[1]^3 + (x[1] + x[2])*y[2]^2 + 0.1*y[3] + (y[4]^2 + y[5]^2)*x[3]*x[4]*x[5],
            (x, y) -> [
                -x[1] - 1, -x[2] - 1, -x[3] - 1, -x[4] - 1, -x[5] - 1,
                x[1] - 1, x[2] - 1, x[3] - 1, x[4] - 1, x[5] - 1,
                y[1]*y[2] - x[1], x[1] - exp(x[2]) + y[3], x[2]*y[1]^2
            ],
            (x, y) -> [
                -y[1] - 1, -y[2] - 1, -y[3] - 1, -y[4] - 1, -y[5] - 1,
                y[1] - 1, y[2] - 1, y[3] - 1, y[4] - 1, y[5] - 1,
                y[1]*y[2] - 0.3, x[1] - 0.2 - y[3]^2, -exp(y[3]) + y[4]*y[5] - 0.1
            ],
            [1.0, −1.0 , −1.0 , −1.0 , −1.0 , −1.0 , −1.0 , 1.0, −1.0, −1.0, 1.0]
        )

    elseif prob_no == 78 || prob_no == "MorganPatrone2006a"
        return BilevelProblem(
            "MorganPatrone2006a",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1.0, 0.0, 1.0],
            (x, y) -> -(x[1] + y[1]),
            (x, y) -> x[1]*y[1],
            (x, y) -> [-x[1] - 0.5, x[1] - 0.5],
            (x, y) -> [-y[1] - 1, y[1] - 1],
            [0.0, 1.0]
        )

    elseif prob_no == 79 || prob_no == "MorganPatrone2006b"
        return BilevelProblem(
            "MorganPatrone2006b",
            [1, 1, 0, 4],
            [1.0, 1.0],
            [-1.25, 0.0, 1.0],
            (x, y) -> -(x[1] + y[1]),
            (x, y) -> begin
                if x[1] >= -0.5 && x[1] <= -0.25
                    (x[1] + 0.25)*y[1]
                elseif x[1] > -0.25 && x[1] < 0.25
                    0.0
                elseif x[1] >= 0.25 && x[1] < 0.51
                    (x[1] - 0.25)*y[1]
                else
                    1e10
                end
            end,
            (x, y) -> Float64[],
            (x, y) -> [-x[1] - 0.5, x[1] - 0.5, -y[1] - 1, y[1] - 1],
            [0.25, 1.0]
        )

    elseif prob_no == 80 || prob_no == "MorganPatrone2006c"
        return BilevelProblem(
            "MorganPatrone2006c",
            [1, 1, 0, 4],
            [1.0, 1.0],
            [-1.0, -0.25, 0.0],
            (x, y) -> -(x[1] + y[1]),
            (x, y) -> begin
                if x[1] >= 7/4 && x[1] <= 2.0
                    (x[1] - 7/4)*y[1]
                elseif x[1] > -7/4 && x[1] < 7/4
                    0.0
                elseif x[1] > -2.0 && x[1] <= -7/4
                    (x[1] + 7/4)*y[1]
                else
                    1e10
                end
            end,
            (x, y) -> Float64[],
            (x, y) -> [-x[1] - 2, x[1] - 2, -y[1] - 1, y[1] - 1],
            [2.0, -1.0]
        )

    elseif prob_no == 81 || prob_no == "MuuQuy2003Ex1"
        return BilevelProblem(
            "MuuQuy2003Ex1",
            [1, 2, 2, 3],
            [2.0, 0.0, 0.0],
            [-2.0769, -0.5868, 2.0],
            (x, y) -> x[1]^2 - 4*x[1] + sum(y.^2),
            (x, y) -> y[1]^2 + y[2]^2/2 + y[1]*y[2] + [1 - 3*x[1], 1 + x[1]]*y,
            (x, y) -> [-x[1], x[1] - 2],
            (x, y) -> [2*y[1] + y[2] - 2*x[1] - 1, -y[1], -y[2]],
            [0.8438, 0.7657, 0.0]
        )

    elseif prob_no == 82 || prob_no == "MuuQuy2003Ex2"
        return BilevelProblem(
            "MuuQuy2003Ex2",
            [2, 3, 3, 4],
            [0.0, 1.0, 1.0, 1.0, 1.0],
            [0.6426, 1.6708, 2.0],
            (x, y) -> -7*x[1] + 4*x[2] + y[1]^2 + y[3]^2 - y[1]*y[3] - 4*y[2],
            (x, y) -> y[1]^2 + y[2]^2/2 + y[3]^2/2 + y[1]*y[2] + (1 - 3*x[1])*y[1] + (1 + x[2])*y[2],
            (x, y) -> [-x[1], -x[2], x[1] + x[2] - 1],
            (x, y) -> [
                2*y[1] + y[2] - y[3] + x[1] - 2*x[2] + 2,
                -y[1], -y[2], -y[3]
            ],
            [0.609, 0.391, 0.0, 0.0, 1.828]
        )

    elseif prob_no == 83 || prob_no == "NieWangYe2017Ex34"
        return BilevelProblem(
            "NieWangYe2017Ex34",
            [1, 2, 2, 2],
            [1.0, 1.0, 1.0],
            [2.0, 0.0, 1.0],
            (x, y) -> x[1] + sum(y),
            (x, y) -> x[1]*sum(y),
            (x, y) -> [-x[1] + 2, x[1] - 3],
            (x, y) -> [
                -y[1]^2 + y[2]^2 + (y[1]^2 + y[2]^2)^2,
                -y[1]
            ],
            [2.0, 0.0, 0.0]
        )

    elseif prob_no == 84 || prob_no == "NieWangYe2017Ex52"
        return BilevelProblem(
            "NieWangYe2017Ex52",
            [2, 3, 5, 2],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [-1.710, -2.232, 1.0],
            (x, y) -> x[1]*y[1] + x[2]*y[2] + prod(x)*prod(y),
            (x, y) -> x[1]*y[1]^2 + x[2]^2*y[2]*y[3] - y[1]*y[3]^2,
            (x, y) -> [-x[1] - 1, -x[2] - 1, x[1] - 1, x[2] - 1, y[1]*y[2] - x[1]^2],
            (x, y) -> [1 - sum(y.^2), sum(y.^2) - 2],
            [−1.0, −1.0, 1.1097, 0.3143, −0.8184]
        )

    elseif prob_no == 85 || prob_no == "NieWangYe2017Ex54"
        return BilevelProblem(
            "NieWangYe2017Ex54",
            [4, 4, 3, 2],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-0.437, -1.190, 1.0],
            (x, y) -> x[1]^2*y[1] + x[2]*y[2] + x[3]*y[3]^2 + x[4]*y[4]^2,
            (x, y) -> y[1]^2 - y[2]*(x[1] + x[2]) - (y[3] + y[4])*(x[3] + x[4]),
            (x, y) -> [
                sum(x.^2) - 1,
                y[1]*y[2] - x[1],
                y[3]*y[4] - x[3]^2
            ],
            (x, y) -> [sum(y.^2) - 1, sum(y.^2) - y[1]^2 - y[1]],
            [0.0, 0.0, −0.7071, −0.7071, 0.6180, 0.0, −0.5559, −0.5559]
        )

    elseif prob_no == 86 || prob_no == "NieWangYe2017Ex57"
        return BilevelProblem(
            "NieWangYe2017Ex57",
            [2, 3, 5, 2],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [-2.0, -1.0, 2.0],
            (x, y) -> x[1]^2*y[1]/2 + x[2]*y[2]^2 - (x[1] + x[2]^2)*y[3],
            (x, y) -> x[2]*(prod(y) + y[2]^2 - y[3]^3),
            (x, y) -> [
                x[1] - 1, x[2] - 1, -x[1] - 1, -x[2] - 1,
                -sum(x) + x[1]^2 + y[1]^2 + y[2]^2 + y[3]^2
            ],
            (x, y) -> [sum(y.^2) - x[1], 2*y[2]*y[3] - 1],
            [1.0, 1.0, 0.0, 0.0, 1.0]
        )

    elseif prob_no == 87 || prob_no == "NieWangYe2017Ex58"
        return BilevelProblem(
            "NieWangYe2017Ex58",
            [4, 4, 3, 2],
            [0.5442,  0.4682,  0.4904,  0.4942, -0.7792, -0.5034, -0.2871, -0.1855],
            [-3.488, -0.862, 2.0],
            (x, y) -> sum(x)*sum(y),
            (x, y) -> x[1]*y[1] + x[2]*y[2] + 0.1*y[3] + 0.5*y[4] - y[3]*y[4],
            (x, y) -> [
                sum(x.^2) - 1,
                y[3]^2 - x[4],
                y[2]*y[4] - x[1]
            ],
            (x, y) -> [
                [1, 2, 3, 4] * (y.^2) - x[1]^2 - x[3]^2 - x[2] - x[4],
                y[2]*y[3] - y[1]*y[4]
            ],
            [0.5135, 0.5050, 0.4882, 0.4929, −0.8346, −0.4104, −0.2106, −0.2887]
        )

    elseif prob_no == 88 || prob_no == "NieWangYe2017Ex61"
        return BilevelProblem(
            "NieWangYe2017Ex61",
            [2, 2, 5, 1],
            [1.0, -1.0, -1.0, 1.0],
            [-1.022, -1.084, 2.0],
            (x, y) -> y[1]^3*(x[1]^2 - 3*x[1]*x[2]) - y[1]^2*y[2] + y[2]*x[2]^3,
            (x, y) -> y[1]*y[2]^2 - y[2]^3 - y[1]^2*(x[2] - x[1]^2),
            (x, y) -> [
                x[1] - 1, x[2] - 1, -x[1] - 1, -x[2] - 1, -y[2] - y[1]*(1 - x[1]^2)
            ],
            (x, y) -> [sum(y.^2) - 1],
            [0.5708, −1.0, −0.1639, 0.9865]
        )

    elseif prob_no == 89 || prob_no == "Outrata1990Ex1a"
        return BilevelProblem(
            "Outrata1990Ex1a",
            [2, 2, 0, 4],
            [0.0, 0.0, 3.0, 3.0],
            [-8.92, -6.05, 2.0],
            (x, y) -> 0.1*sum(x.^2) + 0.5*sum((y .- [3.0, 4.0]).^2) - 12.5,
            (x, y) -> 0.5*(y'*[1.0 -2.0; -2.0 5.0]*y) - x'*y,
            (x, y) -> Float64[],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [0.97, 3.14, 2.6, 1.8]
        )

    elseif prob_no == 90 || prob_no == "Outrata1990Ex1b"
        return BilevelProblem(
            "Outrata1990Ex1b",
            [2, 2, 0, 4],
            [0.0, 0.0, 3.0, 3.0],
            [-7.56, -0.58, 2.0],
            (x, y) -> sum(x.^2) + 0.5*sum((y .- [3.0, 4.0]).^2) - 12.5,
            (x, y) -> 0.5*(y'*[1.0 -2.0; -2.0 5.0]*y) - x'*y,
            (x, y) -> Float64[],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [0.28, 0.48, 2.34, 1.03]
        )

    elseif prob_no == 91 || prob_no == "Outrata1990Ex1c"
        return BilevelProblem(
            "Outrata1990Ex1c",
            [2, 2, 0, 4],
            [0.0, 0.0, 3.0, 3.0],
            [-12.0, -112.71, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2) - 12.5,
            (x, y) -> 0.5*(y'*[1.0 3.0; 3.0 10.0]*y) - x'*y,
            (x, y) -> Float64[],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [20.26, 42.81, 3.0, 3.0]
        )

    elseif prob_no == 92 || prob_no == "Outrata1990Ex1d"
        return BilevelProblem(
            "Outrata1990Ex1d",
            [2, 2, 0, 4],
            [6.0, -3.0, 3.0, 3.0],
            [-3.60, -2.0, 2.0],
            (x, y) -> 0.1*sum(x.^2) + 0.5*sum((y .- [3.0, 4.0]).^2) - 12.5,
            (x, y) -> 0.5*(y'*[1.0 3.0; 3.0 10.0]*y) - x'*y,
            (x, y) -> Float64[],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [2, 0.06, 2.0, 0.0]
        )

    elseif prob_no == 93 || prob_no == "Outrata1990Ex1e"
        return BilevelProblem(
            "Outrata1990Ex1e",
            [2, 2, 0, 4],
            [6.0, -3.0, 3.0, 3.0],
            [-3.15, -16.29, 2],
            (x, y) -> 0.1*sum(x.^2) + 0.5*sum((y .- [3.0, 4.0]).^2) - 12.5,
            (x, y) -> 0.5*(y'*[1.0 3.0; 3.0 10.0]*y) - y'*[-1.0 2.0; 3.0 -3.0]*x,
            (x, y) -> Float64[],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [2.42, −3.65, 0, 1.58]
        )

    elseif prob_no == 94 || prob_no == "Outrata1990Ex2a"
        return BilevelProblem(
            "Outrata1990Ex2a",
            [1, 2, 1, 4],
            [1.0, 1.0, 1.0],
            [0.50, -14.53, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*sum(y.^2) - [3.0 + 1.333*x[1], x[1]]*y,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [2.07, 3.0, 3.0]
        )

    elseif prob_no == 95 || prob_no == "Outrata1990Ex2b"
        return BilevelProblem(
            "Outrata1990Ex2b",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [0.50, -4.50, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(y'*[1.0 + x[1] 0.0; 0.0 0.0]*y) - [3.0 + 1.333*x[1], x[1]]*y,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [0.0, 3.0, 3.0]
        )

    elseif prob_no == 96 || prob_no == "Outrata1990Ex2c"
        return BilevelProblem(
            "Outrata1990Ex2c",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [1.860,  -10.931, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(y'*[1.0 + x[1] 0.0; 0.0 1.0 + 0.1*x[1]]*y) - [3.0 + 1.333*x[1], x[1]]*y,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [3.456, 1.707, 2.569]
        )

    elseif prob_no == 97 || prob_no == "Outrata1990Ex2d"
        return BilevelProblem(
            "Outrata1990Ex2d",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [0.92,  -19.47, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*sum(y.^2) - [3.0 + 1.333*x[1], x[1]]*y,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [2.498, 3.632, 2.8]
        )

    elseif prob_no == 98 || prob_no == "Outrata1990Ex2e"
        return BilevelProblem(
            "Outrata1990Ex2e",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [0.90 -14.94 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(y'*[1.0 + x[1] 0.0; 0.0 1.0]*y) - [3.0 + 1.333*x[1], x[1]]*y,
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] - 2,
                y[1] - 0.333*y[2] - 2,
                -y[1], -y[2]
            ],
            [3.999, 1.665, 3.887]
        )

    elseif prob_no == 99 || prob_no == "Outrata1993Ex31"
        return BilevelProblem(
            "Outrata1993Ex31",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [3.208  -20.531 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(1 + 0.2*x[1])*y[1]^2 + 0.5*(1 + 0.1*x[1])*y[2]^2 - (3 + 1.333*x[1])*y[1] - x[1]*y[2],
            (x, y) -> [-x[1]],
            (x, y) -> [
                (-0.333*y[1] + 0.1*x[1])*y[1] + y[2] + 0.1*x[1] - 2,
                y[1] + (-0.333 - 0.1*x[1])*y[2] + 0.1*x[1] - 2,
                -y[1], -y[2]
            ],
            [1.90910, 2.97836, 2.23182]
        )

    elseif prob_no == 100 || prob_no == "Outrata1993Ex32"
        return BilevelProblem(
            "Outrata1993Ex32",
            [1, 2, 1, 4],
            [0.0, 0.0, 0.0],
            [1.56 -11.68 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(1 + 0.2*x[1])*y[1]^2 + 0.5*(1 + 0.1*x[1])*y[2]^2 - (3 + 1.333*x[1])*y[1] - x[1]*y[2],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] + 0.1*x[1] - 1,
                sum(y.^2) - 0.1*x[1] - 9,
                -y[1], -y[2]
            ],
            [4.06095, 2.68227, 1.48710]
        )

    elseif prob_no == 101 || prob_no == "Outrata1994Ex31"
        return BilevelProblem(
            "Outrata1994Ex31",
            [1, 2, 2, 4],
            [0.0, 0.0, 0.0],
            [3.208,  -20.531, 2.0],
            (x, y) -> 0.5*sum((y .- [3.0, 4.0]).^2),
            (x, y) -> 0.5*(1 + 0.2*x[1])*y[1]^2 + 0.5*(1 + 0.1*x[1])*y[2]^2 - (3 + 1.333*x[1])*y[1] - x[1]*y[2],
            (x, y) -> [x[1] - 10, -x[1]],
            (x, y) -> [
                -0.333*y[1] + y[2] + 0.1*x[1] - 1,
                sum(y.^2) - 0.1*x[1] - 9,
                -y[1],
                -y[2]
            ],
            [4.0604, 2.6822, 1.4871]
        )

    elseif prob_no == 102 || prob_no == "OutrataCervinka2009"
        return BilevelProblem(
            "OutrataCervinka2009",
            [2, 2, 1, 3],
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> -2*x[1] - 0.5*x[2] - y[2],
            (x, y) -> y[1] - y[2] + x[1]*y[1] + x[2]*y[2] + 0.5*y[1]^2 + 0.5*y[2]^2,
            (x, y) -> [x[1]],
            (x, y) -> [
                y[2],
                -y[1] + y[2],
                y[1] + y[2]
            ],
            zeros(Float64, 4)
        )

    elseif prob_no == 103 || prob_no == "PaulaviciusAdjiman2017a"
        return BilevelProblem(
            "PaulaviciusAdjiman2017a",
            [1, 1, 4, 2],
            [1.0, 1.0],
            [0.25, 0.0, 1.0],
            (x, y) -> x[1]^2 + y[1]^2,
            (x, y) -> x[1]*y[1]^2 - 0.5*y[1]^4,
            (x, y) -> [-x[1] - 1, x[1] - 1, -y[1] - 1, y[1] - 1],
            (x, y) -> [-y[1] - 1, y[1] - 1],
            [0.5, 0.0]
        )

    elseif prob_no == 104 || prob_no == "PaulaviciusAdjiman2017b"
        return BilevelProblem(
            "PaulaviciusAdjiman2017b",
            [1, 1, 4, 2],
            [0.0, 0.0],
            [0.25, 0.0, 1.0],
            (x, y) -> x[1] + y[1],
            (x, y) -> 0.5*x[1]*y[1]^2 - x[1]^3*y[1],
            (x, y) -> [-x[1] - 1, x[1] - 1, -y[1] - 1, y[1] - 1],
            (x, y) -> [-y[1] - 1, y[1] - 1],
            [−1.0, −1.0]
        )

    elseif prob_no == 105 || prob_no == "SahinCiric1998Ex2"
        return BilevelProblem(
            "SahinCiric1998Ex2",
            [1, 1, 2, 3],
            [1.0, 1.0],
            [-2.0, -1.5, 1.0],
            (x, y) -> (x[1] - 3)^2 + (y[1] - 2)^2,
            (x, y) -> (y[1] - 5)^2,
            (x, y) -> [-x[1], x[1] - 8],
            (x, y) -> [
                -2*x[1] + y[1] - 1,
                x[1] - 2*y[1] + 2,
                x[1] + 2*y[1] - 14
            ],
            [1.0, 3.0]
        )

    elseif prob_no == 106 || prob_no == "ShimizuAiyoshi1981Ex1"
        return BilevelProblem(
            "ShimizuAiyoshi1981Ex1",
            [1, 1, 3, 3],
            [1.0, 1.0],
            [100.0, 0.0, 1.0],
            (x, y) -> x[1]^2 + (y[1] - 10)^2,
            (x, y) -> (x[1] + 2*y[1] - 30)^2,
            (x, y) -> [x[1] - 15, -x[1] + y[1], -x[1]],
            (x, y) -> [x[1] + y[1] - 20, y[1] - 20, -y[1]],
            [10.0, 10.0]
        )

    elseif prob_no == 107 || prob_no == "ShimizuAiyoshi1981Ex2"
        return BilevelProblem(
            "ShimizuAiyoshi1981Ex2",
            [2, 2, 3, 4],
            [1.0, 1.0, 1.0, 1.0],
            [225.0, 100.0, 1.0],
            (x, y) -> sum((x .- [30.0, 20.0]).^2) - 20*y[1] + 20*y[2],
            (x, y) -> sum((x .- y).^2),
            (x, y) -> [-x[1] - 2*x[2] + 30, x[1] + x[2] - 25, x[2] - 15],
            (x, y) -> [y[1] - 10, y[2] - 10, -y[1], -y[2]],
            [20.0, 5.0, 10.0, 5.0]
        )

    elseif prob_no == 108 || prob_no == "ShimizuEtal1997a"
        return BilevelProblem(
            "ShimizuEtal1997a",
            [1, 1, 0, 3],
            [1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> (x[1] - 5)^2 + (2*y[1] + 1)^2,
            (x, y) -> (y[1] - 1)^2 - 1.5*x[1]*y[1],
            (x, y) -> Float64[],
            (x, y) -> [-3*x[1] + y[1] + 3, x[1] - 0.5*y[1] - 4, x[1] + y[1] - 7],
            [5.0, 2.0]
        )

    elseif prob_no == 109 || prob_no == "ShimizuEtal1997b"
        return BilevelProblem(
            "ShimizuEtal1997b",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [2250.0,  197.75, 1.0],
            (x, y) -> 16*x[1]^2 + 9*y[1]^2,
            (x, y) -> (x[1] + y[1] - 20)^4,
            (x, y) -> [-x[1], -4*x[1] + y[1]],
            (x, y) -> [-y[1], 4*x[1] + y[1] - 50],
            [11.25, 5.0] #Careful, [7.2, 12.8] is a local optimum
        )

    elseif prob_no == 110 || prob_no == "SinhaMaloDeb2014TP3"
        return BilevelProblem(
            "SinhaMaloDeb2014TP3",
            [2, 2, 3, 4],
            [1.0, 1.0, 1.0, 1.0],
            [-18.679, -1.016, 2.0],
            (x, y) -> x[1]^2*(-1) + x[2]^2*(-3) - 4*y[1] + y[2]^2,
            (x, y) -> 2*x[1]^2 + y[1]^2 - 5*y[2],
            (x, y) -> [-x[1], -x[2], x[1]^2 + 2*x[2] - 4],
            (x, y) -> [
                -y[1], -y[2],
                -x[2] + -y[2] + 4*y[1] + 4,
                -x[1]^2 + 2*x[1] - x[2]^2 + 2*y[1] - y[2] - 3
            ],
            -18.6787
        )

    elseif prob_no == 111 || prob_no == "SinhaMaloDeb2014TP6"
        return BilevelProblem(
            "SinhaMaloDeb2014TP6",
            [1, 2, 1, 6],
            [1.0, 1.0, 1.0],
            [-1.209, 7.615, 2.0],
            (x, y) -> (x[1] - 1)^2 - 2*x[1] + 2*y[1],
            (x, y) -> sum((2*y .- [4.0, 1.0]).^2) + x[1]*y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -y[1], -y[2],
                4*x[1] + 5*y[1] + 4*y[2] - 12,
                -4*x[1] - 5*y[1] + 4*y[2] + 4,
                4*x[1] - 4*y[1] + 5*y[2] - 4,
                -4*x[1] + 4*y[1] + 5*y[2] - 4
            ],
            -1.2091
        )

    elseif prob_no == 112 || prob_no == "SinhaMaloDeb2014TP7"
        return BilevelProblem(
            "SinhaMaloDeb2014TP7",
            [2, 2, 4, 4],
            [1.0, 1.0, 0.0, 1.0],
            [-1.961, 1.961, 2.0],
            (x, y) -> -(x[1] + y[1])*(x[2] + y[2])/(1 + x[1]*y[1] + x[2]*y[2]),
            (x, y) -> (x[1] + y[1])*(x[2] + y[2])/(1 + x[1]*y[1] + x[2]*y[2]),
            (x, y) -> [-x[1], -x[2], x[1] - x[2], x[1]^2 + x[2]^2 - 100],
            (x, y) -> [-y[1], -y[2], y[1] - x[2], y[2] - x[1]],
            -1.96
        )

    elseif prob_no == 113 || prob_no == "SinhaMaloDeb2014TP8"
        return BilevelProblem(
            "SinhaMaloDeb2014TP8",
            [2, 2, 5, 6],
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 100.0, 1.0],
            (x, y) -> abs(2*x[1] + 2*x[2] - 3*y[1] - 3*y[2] - 60),
            (x, y) -> sum((y .- x .+ 20).^2),
            (x, y) -> [
                -x[1], -x[2], x[1] - 50, x[2] - 50, x[1] + x[2] + y[1] - 2*y[2] - 40
            ],
            (x, y) -> [
                2*y[1] - x[1] + 10,
                y[1] - 20,
                -y[1] - 10,
                2*y[2] - x[2] + 10,
                y[2] - 20,
                -y[2] - 10
            ],
            0.0
        )

    elseif prob_no == 114 || prob_no == "SinhaMaloDeb2014TP9"
        return BilevelProblem(
            "SinhaMaloDeb2014TP9",
            [10, 10, 0, 20],
            zeros(Float64, 20),
            [0.0, 1.0, 2.0],
            (x, y) -> sum((x .- 1).^2 + y.^2),
            (x, y) -> begin
                t = sqrt.(1:10)
                exp((1 + sum(y.^2)/400 - prod(cos.(y ./ t)))*sum(x.^2))
            end,
            (x, y) -> Float64[],
            (x, y) -> [y .- pi; -y .- pi],
            0.0
        )

    elseif prob_no == 115 || prob_no == "SinhaMaloDeb2014TP10"
        return BilevelProblem(
            "SinhaMaloDeb2014TP10",
            [10, 10, 0, 20],
            zeros(Float64, 20),
            [0.0, 1.0, 2.0],
            (x, y) -> sum((x .- 1).^2 + y.^2),
            (x, y) -> begin
                t = sqrt.(1:10)
                xy = x .* y
                exp(1 + sum(xy.^2)/4000 - prod(cos.(xy ./ t)))
            end,
            (x, y) -> Float64[],
            (x, y) -> [y .- pi; -y .- pi],
            0.0
        )

    elseif prob_no == 116 || prob_no == "TuyEtal2007"
        return BilevelProblem(
            "TuyEtal2007",
            [1, 1, 2, 3],
            [5.0, 0.0],
            [22.5, -1.5, 1.0],
            (x, y) -> x[1]^2 + y[1]^2,
            (x, y) -> -y[1],
            (x, y) -> [-x[1], -y[1]],
            (x, y) -> [
                3*x[1] + y[1] - 15,
                x[1] + y[1] - 7,
                x[1] + 3*y[1] - 15
            ],
            [4.492188, 1.523438]
        )

    elseif prob_no == 117 || prob_no == "Vogel2012"
        return BilevelProblem(
            "Vogel2012",
            [1, 1, 2, 1],
            [1.0, 1.0],
            [1.0, -2.0, 1.0],
            (x, y) -> (y[1] + 1)^2,
            (x, y) -> y[1]^3 - 3*y[1],
            (x, y) -> [-3 - x[1], x[1] - 2],
            (x, y) -> [x[1] - y[1]],
            [-2.0, -2.0]
        )

    elseif prob_no == 118 || prob_no == "WanWangLv2011"
        return BilevelProblem(
            "WanWangLv2011",
            [2, 3, 0, 8],
            [0.0, 0.5, 0.0, 0.0, 0.0],
            [10.62, -0.50, 1.0],
            (x, y) -> (1 + [1.0, -1.0]'*x + 2*y[2])*(8 - x[1] + [-2.0, 1.0, 5.0]'*y),
            (x, y) -> [2.0, -1.0, 1.0]'*y,
            (x, y) -> Float64[],
            (x, y) -> [
                [-1.0, 1.0, 1.0]'*y - 1,
                2*x[1] + [-1.0, 2.0, -0.5]'*y - 1,
                2*x[2] + [2.0, -1.0, -0.5]'*y - 1,
                -x[1], -x[2], -y[1], -y[2], -y[3]
            ],
            [0.0, 0.75, 0.0, 0.5, 0.0]
        )

    elseif prob_no == 119 || prob_no == "YeZhu2010Ex42"
        return BilevelProblem(
            "YeZhu2010Ex42",
            [1, 1, 2, 1],
            [-1.0, -1.0],
            [1.0, -2.0, 1.0],
            (x, y) -> (x[1] - 1)^2 + y[1]^2,
            (x, y) -> y[1]^3 - 3*y[1],
            (x, y) -> [-3 - x[1], x[1] - 2],
            (x, y) -> [x[1] - y[1]],
            [1.0, 1.0]
        )

    elseif prob_no == 120 || prob_no == "YeZhu2010Ex43"
        return BilevelProblem(
            "YeZhu2010Ex43",
            [1, 1, 2, 1],
            [-1.0, -1.0],
            [1.0, -2.0, 1.0],
            (x, y) -> (x[1] - 0.5)^2 + (y[1] - 2)^2,
            (x, y) -> y[1]^3 - 3*y[1],
            (x, y) -> [-3 - x[1], x[1] - 2],
            (x, y) -> [x[1] - y[1]],
            [1.0, 1.0]
        )

    elseif prob_no == 121 || prob_no == "Yezza1996Ex31"
        return BilevelProblem(
            "Yezza1996Ex31",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [1.5, -2.5, 1.0],
            (x, y) -> -(4*x[1] - 3)*y[1] + 2*x[1] + 1,
            (x, y) -> -(1 - 4*x[1])*y[1] - 2*x[1] - 2,
            (x, y) -> [-x[1], x[1] - 1],
            (x, y) -> [-y[1], y[1] - 1],
            [0.25, 0.0]
        )

    elseif prob_no == 122 || prob_no == "Yezza1996Ex41"
        return BilevelProblem(
            "Yezza1996Ex41",
            [1, 1, 0, 2],
            [-1.0, -1.0],
            [0.5, 2.5, 1.0],
            (x, y) -> 0.5*(y[1] - 2)^2 + 0.5*(x[1] - y[1] - 2)^2,
            (x, y) -> 0.5*y[1]^2 + x[1] - y[1],
            (x, y) -> Float64[],
            (x, y) -> [-y[1], y[1] - x[1]],
            [3.0, 1.0]
        )

    elseif prob_no == 123 || prob_no == "Zlobec2001a"
        return BilevelProblem(
            "Zlobec2001a",
            [1, 2, 0, 4], # /!\ changed regarding the original one as the constraint dimensions does not fit
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            (x, y) -> -y[1]/x[1],
            (x, y) -> -y[1] - y[2],
            (x, y) -> Float64[],
            (x, y) -> [
                -1 + y[1],
                x[1] + y[1],
                -y[1],
                -y[2]
            ],
            [1.0, 1.0, 0.0]
        )

    elseif prob_no == 124 || prob_no == "Zlobec2001b"
        @warn "Zlobec2001b is not a closed problem, no optimal solution is available"
        return BilevelProblem(
            "Zlobec2001b",
            [1, 1, 2, 4],
            [1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> x[1] + y[1],
            (x, y) -> -y[1],
            (x, y) -> [-x[1], x[1] - 1],
            (x, y) -> [-y[1], y[1] - 1, -x[1]*y[1], x[1]*y[1]],
            [] # No optimal solution as problem not closed
        )

    elseif prob_no == 125 || prob_no == "DesignCentringP1"
        return BilevelProblem(
            "DesignCentringP1",
            [3, 6, 3, 3],
            ones(Float64, 9),
            [NaN, NaN, 0.0],
            (x, y) -> -π*x[3]^2,
            (x, y) -> y[1] + y[2]^2 - y[3]/4 - y[4] + y[6],
            (x, y) -> [-y[1] - y[2]^2, y[3]/4 + y[4] - 3/4, -y[6] - 1],
            (x, y) -> [
                (y[1] - x[1])^2 + (y[2] - x[2])^2 - x[3]^2,
                (y[3] - x[1])^2 + (y[4] - x[2])^2 - x[3]^2,
                (y[5] - x[1])^2 + (y[6] - x[2])^2 - x[3]^2
            ],
        [0.7486, −0.2304, 0.7696, −0.0084, −0.0917, 0.9352, 0.5162, 0.7486, −1.0]
        )

    elseif prob_no == 126 || prob_no == "DesignCentringP2"
        return BilevelProblem(
            "DesignCentringP2",
            [4, 6, 5, 3],
            ones(Float64, 10),
            [NaN, NaN, 0.0],
            (x, y) -> -π*x[3]*x[4],
            (x, y) -> y[1] + y[2]^2 - y[3]/4 - y[4] + y[6],
            (x, y) -> [
                -y[1] - y[2]^2,
                y[3]/4 + y[4] - 3/4,
                -y[6] - 1,
                x[3] - 1,
                x[4] - 1
            ],
            (x, y) -> [
                (x[1] - y[1])^2 / x[3]^2 + (x[2] - y[2])^2 / x[4]^2 - 1,
                (x[1] - y[3])^2 / x[3]^2 + (x[2] - y[4])^2 / x[4]^2 - 1,
                (x[1] - y[5])^2 / x[3]^2 + (x[2] - y[6])^2 / x[4]^2 - 1
            ],
            [3.0, 0.0, 1.0, 1.0, 3.0, 0.0, 3.0, 0.0, 3.0, 0.0]
        )

    elseif prob_no == 127 || prob_no == "DesignCentringP3"
        return BilevelProblem(
            "DesignCentringP3",
            [6, 6, 3, 3],
            [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> -π*abs(x[3]*x[6] - x[4]*x[5]),
            (x, y) -> y[1] + y[2]^2 - y[3]/4 - y[4] + y[6],
            (x, y) -> [
                -y[1] - y[2]^2,
                y[3]/4 + y[4] - 3/4,
                -y[6] - 1
            ],
            (x, y) -> begin
                z = [x[1]; x[2]]
                A = inv([x[3] x[4]; x[5] x[6]] * [x[3] x[5]; x[4] x[6]])
                [
                    ([y[1]; y[2]] - z)' * A * ([y[1]; y[2]] - z) - 1,
                    ([y[3]; y[4]] - z)' * A * ([y[3]; y[4]] - z) - 1,
                    ([y[5]; y[6]] - z)' * A * ([y[5]; y[6]] - z) - 1
                ]
            end,
            3.7234
        )

    elseif prob_no == 128 || prob_no == "DesignCentringP4"
        return BilevelProblem(
            "DesignCentringP4",
            [4, 6, 3, 12],
            ones(Float64, 10),
            [NaN, NaN, 0.0],
            (x, y) -> -(x[1] - x[3]) * (x[2] - x[4]),
            (x, y) -> y[1] + y[2]^2 - y[3]/4 - y[4] + y[6],
            (x, y) -> [
                -y[1] - y[2]^2,
                y[3]/4 + y[4] - 3/4,
                -y[6] - 1
            ],
            (x, y) -> [y .- [x[1], x[2], x[1], x[2], x[1], x[2]] ; 
                    y .- [x[3], x[4], x[3], x[4], x[3], x[4]]],
            [−1.0, 1.0, −1.0, 1.0, −1.0, 1.0, −1.0, 1.0, −1.0, 1.0]
        )
    # /!\ TODO: check NetworkDesign problems in details

    elseif prob_no == 129 || prob_no == "NetworkDesignP1"
        return BilevelProblem(
            "NetworkDesignP1",
            [5, 5, 5, 11],
            zeros(Float64, 10),
            [300.5, 419.8, 2.0],
            (x, y) -> begin
                xy = y ./ (1 .+ x)
                [50 + xy[1], 10*xy[2], 10 + xy[3], 10*xy[4], 50 + xy[5]] * y + 100*sum(x)
            end,
            (x, y) -> begin
                xy = y ./ (1 .+ x) ./ 2
                [50 + xy[1], 10*xy[2], 10 + xy[3], 10*xy[4], 50 + xy[5]] * y
            end,
            (x, y) -> -x .- 1,
            (x, y) -> begin
                A = [1 0 1 0 1; 0 1 -1 0 -1; -1 0 -1 1 0]
                vcat(A * y .+ [-6; 0; 0], -A * y .+ [6; 0; 0], -y)
            end,
            300.5
        )

    elseif prob_no == 130 || prob_no == "NetworkDesignP2"
        return BilevelProblem(
            "NetworkDesignP2",
            [5, 5, 5, 11],
            zeros(Float64, 10),
            [142.9, 81.95, 2.0],
            (x, y) -> begin
                xy = y ./ (1 .+ x)
                [50 + xy[1], 10*xy[2], 10 + xy[3], 10*xy[4], 50 + xy[5]] * y + sum(x)
            end,
            (x, y) -> begin
                xy = y ./ (1 .+ x) ./ 2
                [50 + xy[1], 10*xy[2], 10 + xy[3], 10*xy[4], 50 + xy[5]] * y
            end,
            (x, y) -> -x .- 1,
            (x, y) -> begin
                A = [1 0 1 0 1; 0 1 -1 0 -1; -1 0 -1 1 0]
                vcat(A * y .+ [-6; 0; 0], -A * y .+ [6; 0; 0], -y)
            end,
            142.9
        )

        # Example for N = 10
    elseif prob_no == 131 || prob_no == "RobustPortfolioP1"
        @warn "RobustPortfolioP1 is encoded as a special case where N=10 and δ=2."
        return BilevelProblem(
            "RobustPortfolioP1",
            [11, 10, 13, 11], # Initial point depends on N
            ones(Float64, 21),
            [1.15, 0.0, 2.0],
            (x, y) -> -x[end],
            (x, y) -> y' * x[1:end-1] - x[end],
            (x, y) -> begin
                N = length(y)
                [x[end] - y' * x[1:N]; -x[1:N]; sum(x[1:N]) - 1; 1 - sum(x[1:N])]
            end,
            (x, y) -> begin
                N = length(y)
                I = collect(1:N)
                d = 2
                si = ((0.05 / 3 / N) * sqrt(2 * N * (N + 1) .* I)).^d
                yi = 1.15 .+ (0.05 / N) .* I
                [sum(abs.(y .- yi).^d ./ si) - 1.5^d; -y]
            end,
            [fill(1/10, 10) ; fill(1.15, 11)]
        )

    elseif prob_no == 132 || prob_no == "RobustPortfolioP2"
        @warn "RobustPortfolioP2 is encoded as a special case where N=10"
        return BilevelProblem(
            "RobustPortfolioP2",
            [11, 10, 13, 11], # N is problem parameter
            ones(Float64, 21), # Initial point depends on N
            [1.15, 0.0, 2.0],
            (x, y) -> -x[end],
            (x, y) -> y' * x[1:end-1] - x[end],
            (x, y) -> begin
                N = length(y)
                [x[end] - y' * x[1:N]; -x[1:N]; sum(x[1:N]) - 1; 1 - sum(x[1:N])]
            end,
            (x, y) -> begin
                N = length(y)
                I = collect(1:N)
                si = ((0.05 / 3 / N)^2 * (2 * N * (N + 1) .* I))
                yi = 1.15 .+ (0.05 / N) .* I
                sx = 1.5 * (1 + sum((x[1:N] .- 1 / N).^2))
                [sum((y .- yi).^2 ./ si) - sx^2; -y]
            end,
            [fill(1/10, 10) ; fill(1.15, 11)]
        )

    elseif prob_no == 133 || prob_no == "TollSettingP1"
        return BilevelProblem(
            "TollSettingP1",
            [3, 8, 3, 18],
            [0.0, 5.0, 5.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [-7.0, 12.0, 2.0],
            (x, y) -> -[y[3], y[4], y[8]]' * x,
            (x, y) -> [2, 6, 5 + x[1], x[2], 4, 2, 6, x[3]]' * y,
            (x, y) -> -x,
            (x, y) -> begin
                A = [1 1 1 0 0 0 0 0; -1 0 0 1 1 0 0 0;
                    0 -1 0 -1 0 1 1 0; 0 0 0 0 -1 -1 0 1; 0 0 1 0 0 0 1 1]
                vcat(A * y .+ [-1; 0; 0; 0; -1], -A * y .+ [1; 0; 0; 0; 1], -y)
            end,
            [7.0, 4.0, 6.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        )

    elseif prob_no == 134 || prob_no == "TollSettingP2"
        return BilevelProblem(
            "TollSettingP2",
            [3, 18, 3, 38],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0],
            [-4.5, 32.0, 2.0],
            (x, y) -> -[y[1] + y[2], y[3] + y[4], y[5] + y[6]]' * x,
            (x, y) -> [2*x[1], 2*x[1], 2*x[2], 2*x[2], 2*x[3], 2*x[3], 5, 7, 14, 7, 2, 4, 29, 20, 12, 8, 5, 2]' * y,
            (x, y) -> -x,
            (x, y) -> begin
                rows1 = [fill(1,3); fill(2,3); fill(3,3); fill(4,3); fill(5,3); fill(6,3); fill(7,2); fill(8,2); fill(9,1); fill(10,1)]
                rows_minus1 = [5; 6; 7; 7; 8; 8; fill(9, 3); fill(10, 3)]
                cols1 = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 5, 13, 2, 6, 16, 3, 14, 4, 17, 15, 18]
                cols_minus1 = [7, 10, 1, 8, 2, 11, 3, 5, 9, 4, 6, 12]
                A = zeros(10, 18)
                for i in rows1
                    for j in cols1
                        A[i, j] = 1.0
                    end
                end
                for i in rows_minus1
                    for j in cols_minus1
                        A[i, j] = - 1.0
                    end
                end
                vcat(A * y .+ [fill(-1.0, 4) ; zeros(Float64, 6)], -A * y .+ [fill(-1.0, 4) ; zeros(Float64, 6)], -y)
            end,
            [0.5, 4.0, 4.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
        )

    elseif prob_no == 135 || prob_no == "TollSettingP3"
        return BilevelProblem(
            "TollSettingP3",
            [3, 18, 3, 38],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0],
            [-3.5, 32.0, 2.0],
            (x, y) -> -[y[1] + 10*y[2], y[3] + 10*y[4], y[5] + 10*y[6]]' * x,
            (x, y) -> [2*x[1], 20*x[1], 2*x[2], 20*x[2], 2*x[3], 20*x[3], 5, 7, 14, 7, 2, 4, 29, 20, 12, 8, 5, 2]' * y,
            (x, y) -> -x,
            (x, y) -> begin
                rows1 = [fill(1,3); fill(2,3); fill(3,3); fill(4,3); fill(5,3); 6 ; 7 ; 7 ; 8 ; 9 ; 10]
                rows_minus1 = [5 ; 6 ; 7 ; 7 ; 8 ; 9 ; 9 ; 9 ; 10]
                rows10 = [6 ; 6 ; 7]
                rows_minus10 = [8 ; 10 ; 10]

                cols1 = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 5, 13, 16, 3, 14, 17, 15, 18]
                cols_minus1 = [7, 10, 1, 8, 11, 3, 5, 9, 12]
                cols10 = [2, 6, 4]
                cols_minus10 = [2, 4, 6]

                # building matrix A
                A = zeros(10, 18)
                for i in rows1
                    for j in cols1
                        A[i, j] = 1.0
                    end
                end
                for i in rows_minus1
                    for j in cols_minus1
                        A[i, j] = - 1.0
                    end
                end
                for i in rows10
                    for j in cols10
                        A[i, j] = 10.0
                    end
                end
                for i in rows_minus10
                    for j in cols_minus10
                        A[i, j] = -10.0
                    end
                end
                vcat(A * y .+ [fill(-1.0, 4) ; zeros(Float64, 6)], -A * y .+ [fill(-1.0, 4) ; zeros(Float64, 6)], -y)
            end,
            -24.0
        )

    elseif prob_no == 136 || prob_no == "TollSettingP4"
        return BilevelProblem(
            "TollSettingP4",
            [2, 4, 0, 8],
            ones(Float64, 6),
            [-4.0, 14.0, 2.0],
            (x, y) -> -(y[2] + y[3]) * x[1] - y[3] * x[2],
            (x, y) -> [8, 3 + 2*x[1], 3 + 2*x[1] + 2*x[2], 6]' * y,
            (x, y) -> Float64[],
            (x, y) -> begin
                A = [1 1 0 0; 0 0 1 1]
                vcat(A * y .+ [-1; -1], -A * y .+ [1; 1], -y)
            end,
            -8.0
        )

    elseif prob_no == 137 || prob_no == "TollSettingP5"
        return BilevelProblem(
            "TollSettingP5",
            [1, 4, 0, 8],
            ones(Float64, 5),
            [-2.5, 14.0, 2.0],
            (x, y) -> -(y[2] + y[3]) * x[1],
            (x, y) -> [8, 3 + 2*x[1], 4 + 2*x[1], 6]' * y,
            (x, y) -> Float64[],
            (x, y) -> begin
                A = [1 1 0 0; 0 0 1 1]
                vcat(A * y .+ [-1; -1], -A * y .+ [1; 1], -y)
            end,
            -2.5
        )

    # Special cases that need additional data
    elseif prob_no == 138 || prob_no == "OptimalControl"
        error("OptimalControl requires additional parameters - implement separately")

    ## ---------------- LINEAR BILEVEL PROBLEMS ---------------- ##
        # ...existing code...

    elseif prob_no == 139 || prob_no == "AnandalinghamWhite1990"
        return BilevelProblem(
            "AnandalinghamWhite1990",
            [1, 1, 1, 6],
            [1.0, 1.0],
            [-49.0, 15.0, 1.0],
            (x, y) -> -x[1] - 3*y[1],
            (x, y) -> -x[1] + 3*y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -x[1] - 2*y[1] + 10,
                x[1] - 2*y[1] - 6,
                2*x[1] - y[1] - 21,
                x[1] + 2*y[1] - 38,
                -x[1] + 2*y[1] - 18,
                -y[1]
            ],
            [16.0, 11.0]
        )

    elseif prob_no == 140 || prob_no == "Bard1984a"
        return BilevelProblem(
            "Bard1984a",
            [1, 1, 1, 5],
            [1.0, 1.0],
            [28/9, -60/9, 1.0],
            (x, y) -> x[1] + y[1],
            (x, y) -> -5*x[1] - y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -x[1] - 0.5*y[1] + 2,
                -0.25*x[1] + y[1] - 2,
                x[1] + 0.5*y[1] - 8,
                x[1] - 2*y[1] - 4,
                -y[1]
            ],
            [8/9, 20/9]
        )

    elseif prob_no == 141 || prob_no == "Bard1984b"
        return BilevelProblem(
            "Bard1984b",
            [1, 1, 1, 5],
            [1.0, 1.0],
            [-37.6, 1.6, 1.0],
            (x, y) -> -5*x[1] - y[1],
            (x, y) -> y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -x[1] - 0.5*y[1] + 2,
                -0.25*x[1] + y[1] - 2,
                x[1] + 0.5*y[1] - 8,
                x[1] - 2*y[1] - 4,
                -y[1]
            ],
            [7.2, 1.6]
        )

    elseif prob_no == 142 || prob_no == "Bard1991Ex2"
        return BilevelProblem(
            "Bard1991Ex2",
            [1, 2, 1, 5],
            [0.0, 0.0, 0.0],
            [-1.0, -1.0, 1.0],
            (x, y) -> -x[1] + 10*y[1] - y[2],
            (x, y) -> -y[1] - y[2],
            (x, y) -> [-x[1]],
            (x, y) -> [
                x[1] + y[1] - 1,
                x[1] + y[2] - 1,
                y[1] + y[2] - 1,
                -y[1],
                -y[2]
            ],
            [0.0, 0.0, 1.0]
        )

    elseif prob_no == 143 || prob_no == "BardFalk1982Ex2"
        return BilevelProblem(
            "BardFalk1982Ex2",
            [2, 2, 2, 5],
            [1.0, 1.0, 1.0, 1.0],
            [-3.25, -4.0, 1.0],
            (x, y) -> -2*x[1] + x[2] + 0.5*y[1],
            (x, y) -> x[1] + x[2] - 4*y[1] + y[2],
            (x, y) -> -x,
            (x, y) -> [
                -2*x[1] + y[1] - y[2] + 2.5,
                x[1] - 3*x[2] + y[2] - 2,
                x[1] + x[2] - 2,
                -y[1],
                -y[2]
            ],
            [2.0, 0.0, 1.5, 0.0]
        )

    elseif prob_no == 144 || prob_no == "BenAyedBlair1990a"
        return BilevelProblem(
            "BenAyedBlair1990a",
            [1, 2, 2, 4],
            [1.0, 1.0, 1.0],
            [-2.5, -5.0, 1.0],
            (x, y) -> -1.5*x[1] - 6*y[1] - y[2],
            (x, y) -> -y[1] - 5*y[2],
            (x, y) -> [-x[1], x[1] - 1],
            (x, y) -> [
                x[1] + 3*y[1] + y[2] - 5,
                2*x[1] + y[1] + 3*y[2] - 5,
                -y[1],
                -y[2]
            ],
            [1.0, 0.0, 1.0]
        )

    elseif prob_no == 145 || prob_no == "BenAyedBlair1990b"
        return BilevelProblem(
            "BenAyedBlair1990b",
            [1, 1, 1, 4],
            [3.0, 4.0],
            [-6.0, 5.0, 1.0],
            (x, y) -> -x[1] - y[1],
            (x, y) -> y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -4*x[1] - 3*y[1] + 19,
                x[1] + 2*y[1] - 11,
                3*x[1] + y[1] - 13,
                -y[1]
            ],
            [1.0, 5.0]
        )

    elseif prob_no == 146 || prob_no == "BialasKarwan1984a"
        return BilevelProblem(
            "BialasKarwan1984a",
            [1, 2, 1, 7],
            [1.0, 1.0, 1.0],
            [-2.0, -0.5, 1.0],
            (x, y) -> -x[1] - y[2],
            (x, y) -> -y[2],
            (x, y) -> [-x[1]],
            (x, y) -> [
                x[1] + y[1] + y[2] - 3,
                -x[1] - y[1] + y[2] + 1,
                -x[1] + y[1] + y[2] - 1,
                x[1] - y[1] + y[2] - 1,
                y[2] - 0.5,
                -y[1],
                -y[2]
            ],
            [1.5, 1.0, 0.5]
        )

    elseif prob_no == 147 || prob_no == "BialasKarwan1984b"
        return BilevelProblem(
            "BialasKarwan1984b",
            [1, 1, 1, 6],
            [10.0, 14.0],
            [-11.0, 11.0, 1.0],
            (x, y) -> -y[1],
            (x, y) -> y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -x[1] - 2*y[1] + 10,
                x[1] + 2*y[1] - 38,
                -x[1] + 2*y[1] - 18,
                x[1] - 2*y[1] - 6,
                2*x[1] - y[1] - 21,
                -y[1]
            ],
            [16.0, 11.0]
        )

    elseif prob_no == 148 || prob_no == "CandlerTownsley1982"
        return BilevelProblem(
            "CandlerTownsley1982",
            [2, 3, 2, 6],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [-29.2, 3.2, 1.0],
            (x, y) -> -[8.0, 4.0]'*x + [4.0, -40.0, -4.0]'*y,
            (x, y) -> [1.0, 2.0]'*x + [1.0, 1.0, 2.0]'*y,
            (x, y) -> -x,
            (x, y) -> [
                -y[1] + y[2] + y[3] + 1,
                2*x[1] - y[1] + 2*y[2] - 0.5*y[3] - 1,
                2*x[2] + 2*y[1] - y[2] - 0.5*y[3] - 1,
                -y[1],
                -y[2],
                -y[3]
            ],
            [0.0, 0.9, 0.0, 0.6, 0.4]
        )

    elseif prob_no == 149 || prob_no == "ClarkWesterberg1988"
        return BilevelProblem(
            "ClarkWesterberg1988",
            [1, 1, 0, 3],
            [20.0, 0.0],
            [-37.0, 14.0, 1.0],
            (x, y) -> x[1] - 4*y[1],
            (x, y) -> y[1],
            (x, y) -> Float64[],
            (x, y) -> [
                -2*x[1] + y[1],
                2*x[1] + 5*y[1] - 108,
                2*x[1] - 3*y[1] + 4
            ],
            [19.0, 14.0]
        )

    elseif prob_no == 150 || prob_no == "ClarkWesterberg1990b"
        return BilevelProblem(
            "ClarkWesterberg1990b",
            [1, 2, 2, 5],
            [1.0, 1.0, 1.0],
            [-13.0, -4.0, 1.0],
            (x, y) -> -x[1] - 3*y[1] + 2*y[2],
            (x, y) -> -y[1],
            (x, y) -> [-x[1], x[1] - 8],
            (x, y) -> [
                0*x[1] - y[1] + 0*y[2] + 0,
                0*x[1] + y[1] + 0*y[2] - 4,
                -2*x[1] + y[1] + 4*y[2] - 16,
                8*x[1] + 3*y[1] - 2*y[2] - 48,
                -2*x[1] + y[1] - 3*y[2] + 12
            ],
            [5.0, 4.0, 2.0]
        )

    elseif prob_no == 151 || prob_no == "GlackinEtal2009"
        return BilevelProblem(
            "GlackinEtal2009",
            [2, 1, 3, 3],
            [0.0, 1.0, 0.0],
            [6.0, 0.0, 1.0],
            (x, y) -> -2*x[1] + 4*x[2] + 3*y[1],
            (x, y) -> -y[1],
            (x, y) -> [x[1] - x[2] + 1, -x[1], -x[2]],
            (x, y) -> [
                x[1] + x[2] + y[1] - 4,
                2*x[1] + 2*x[2] + y[1] - 6,
                -y[1]
            ],
            [1.0, 2.0, 0.0]
        )

    elseif prob_no == 152 || prob_no == "HaurieSavardWhite1990"
        return BilevelProblem(
            "HaurieSavardWhite1990",
            [1, 1, 0, 4],
            [1.0, 1.0],
            [27.0, -3.0, 1.0],
            (x, y) -> x[1] + 5*y[1],
            (x, y) -> -y[1],
            (x, y) -> Float64[],
            (x, y) -> [
                -3*x[1] + 2*y[1] - 6,
                3*x[1] + 4*y[1] - 48,
                2*x[1] - 5*y[1] - 9,
                -x[1] - y[1] + 8
            ],
            [12.0, 3.0]
        )

    elseif prob_no == 153 || prob_no == "HuHuangZhang2009"
        return BilevelProblem(
            "HuHuangZhang2009",
            [1, 2, 1, 5],
            [1.0, 1.0, 1.0],
            [-76/9, -41/9, 1.0],
            (x, y) -> -4*x[1] - y[1] - y[2],
            (x, y) -> -x[1] - 3*y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                x[1] + y[1] + y[2] - 25/9,
                x[1] + y[2] - 2,
                y[1] + y[2] - 8/9,
                -y[1],
                -y[2]
            ],
            [17/9, 8/9, 0.0]
        )

    elseif prob_no == 154 || prob_no == "LanWenShihLee2007"
        return BilevelProblem(
            "LanWenShihLee2007",
            [1, 1, 1, 7],
            [1.0, 1.0],
            [-85.0855, 50.174, 2.0],
            (x, y) -> 2*x[1] - 11*y[1],
            (x, y) -> x[1] + 3*y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                1*x[1] - 2*y[1] - 4,
                2*x[1] - y[1] - 24,
                3*x[1] + 4*y[1] - 96,
                1*x[1] + 7*y[1] - 126,
                -4*x[1] + 5*y[1] - 65,
                -1*x[1] - 4*y[1] + 8,
                -y[1]
            ],
            [17.45, 10.908]
        )

    elseif prob_no == 155 || prob_no == "LiuHart1994"
        return BilevelProblem(
            "LiuHart1994",
            [1, 1, 1, 4],
            [1.0, 1.0],
            [-16.0, 4.0, 1.0],
            (x, y) -> -x[1] - 3*y[1],
            (x, y) -> y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                -x[1] + y[1] - 3,
                x[1] + 2*y[1] - 12,
                4*x[1] - y[1] - 12,
                -y[1]
            ],
            [4.0, 4.0]
        )

    elseif prob_no == 156 || prob_no == "MershaDempe2006Ex1"
        @warn "MershaDempe2006Ex1 has no optimal solution."
        return BilevelProblem(
            "MershaDempe2006Ex1",
            [1, 1, 1, 5],
            [1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> x[1] - 8*y[1],
            (x, y) -> y[1],
            (x, y) -> [-x[1]],
            (x, y) -> [
                5*x[1] - 2*y[1] - 33,
                -x[1] - 2*y[1] + 9,
                -7*x[1] + 3*y[1] - 5,
                x[1] + y[1] - 15,
                -y[1]
            ],
            []
        )

    elseif prob_no == 157 || prob_no == "MershaDempe2006Ex2"
        return BilevelProblem(
            "MershaDempe2006Ex2",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-20.0, -6.0, 1.0],
            (x, y) -> -x[1] - 2*y[1],
            (x, y) -> -y[1],
            (x, y) -> [-2*x[1] + 3*y[1] - 12, x[1] + y[1] - 14],
            (x, y) -> [
                -3*x[1] + y[1] + 3,
                3*x[1] + y[1] - 30
            ],
            [8.0, 6.0]
        )

    elseif prob_no == 158 || prob_no == "TuyEtal1993"
        return BilevelProblem(
            "TuyEtal1993",
            [2, 2, 3, 4],
            [0.0, 0.0, 0.0, 0.0],
            [-3.25, -6.0, 1.0],
            (x, y) -> -2*x[1] + x[2] + 0.5*y[1],
            (x, y) -> -4*y[1] + y[2],
            (x, y) -> [x[1] + x[2] - 2, -x[1], -x[2]],
            (x, y) -> [
                -2*x[1] + y[1] - y[2] + 2.5,
                x[1] - 3*x[2] + y[2] - 2,
                -y[1],
                -y[2]
            ],
            [2.0, 0.0, 1.5, 0.0]
        )

    elseif prob_no == 159 || prob_no == "TuyEtal1994"
        return BilevelProblem(
        "GumusFloudas2001Ex1",
        [1, 1, 3, 3],
        [1.0, 1.0],
        [2250.0, 197.75, 1.0],
        (x,y) -> (x[1]-30)^2 + (y[1]-20)^2,
        (x,y) -> (y[1]-x[1])^2,
        (x,y) -> [
            x[1] + y[1] - 50,
            -x[1] - y[1] + 50,
            x[1]^2 + y[1]^2 - 100
        ],
        (x,y) -> [
            y[1] - 50,
            -y[1],
            x[1]^2 - y[1]^2
        ]
    )
    end

    # Special cases that need additional data
    if prob_no == 138 || prob_no == "OptimalControl"
        error("OptimalControl requires additional parameters - implement separately")
    elseif prob_no == 173 || prob_no == "ShehuEtal2019Ex42"
        error("ShehuEtal2019Ex42 requires additional parameters - implement separately")
    else
        error("Problem $prob_no not implemented yet")
    end
end