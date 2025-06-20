#using LinearAlgebra

#export get_bilevel_problem

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

# Problem: AiyoshiShimizu1984Ex2
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
            ],
            [1.0, 0.0]
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
            [1, 2, 1, 6],
            [1.0, 1.0, 1.0],
            [2.75 0.57 2],
            (x, y) -> (x[1] - y[2])^4 + (y[1] - 1)^2 + (y[1] - y[2])^2,
            (x, y) -> 2*x[1] + exp(y[1]) + y[1]^2 + 4*y[1] + 2*y[2]^2 - 6*y[2],
            (x, y) -> [-x[1]],
            (x, y) -> [
                6*x[1] + y[1]^2 + exp(y[2]) - 15,
                5*x[1] + y[1]^4 - y[2] - 25,
                y[1] - 4,
                y[2] - 2,
                -y[1],
                -y[2]
            ],
            [1.94, 0.0, 1.21]
        )

    elseif prob_no == 19 || prob_no == "Dempe1992a"
        return BilevelProblem(
            "Dempe1992a",
            [2, 2, 1, 2],
            [1.0, 1.0, 1.0, 1.0],
            [NaN, NaN, 0.0],
            (x, y) -> y[2],
            (x, y) -> 0.5*sum((y .- [1.0, 0.0]).^2),
            (x, y) -> [sum((x .+ [0.0, 1.0]).^2) - 1],
            (x, y) -> [
                y[1] + y[2]*x[1] + x[2],
                y[1]
            ],
            [0.0, 0.0, 0.0, −0.5]
        )

    elseif prob_no == 20 || prob_no == "Dempe1992b"
        return BilevelProblem(
            "Dempe1992b",
            [1, 1, 0, 1],
            [-1.0, -1.0],
            [31.25, 4.0, 1.0],
            (x, y) -> (x[1] - 3.5)^2 + (y[1] + 4)^2,
            (x, y) -> (y[1] - 3.0)^2,
            (x, y) -> Float64[],
            (x, y) -> [y[1]^2 - x[1]],
            31.25
        )
    # ...existing code...

    elseif prob_no == 21 || prob_no == "DempeDutta2012Ex24"
        return BilevelProblem(
            "DempeDutta2012Ex24",
            [1, 1, 0, 1],
            [1.0, 1.0],
            [0.0, 0.0, 1.0],
            (x, y) -> (x[1] - 1)^2 + y[1]^2,
            (x, y) -> x[1]^2 * y[1],
            (x, y) -> Float64[],
            (x, y) -> [y[1]^2],
            [1.0, 0.0]
        )

    elseif prob_no == 22 || prob_no == "DempeDutta2012Ex31"
        return BilevelProblem(
            "DempeDutta2012Ex31",
            [2, 2, 4, 2],
            [1.0, 1.0, 1.0, 1.0],
            [-1.0, 4.0, 1.0],
            (x, y) -> -y[2],
            (x, y) -> sum((y .+ [0.0, 1.0]).^2),
            (x, y) -> [
                -x[1],
                -x[2],
                y[1]*y[2],
                -y[1]*y[2]
            ],
            (x, y) -> [
                (y[1] - x[1])^2 + (y[2] - x[1] - 1)^2 - 1,
                (y[1] + x[2])^2 + (y[2] - x[2] - 1)^2 - 1
            ],
            [0.71, 0.71, 0.0, 1.0]
        )

    elseif prob_no == 23 || prob_no == "DempeEtal2012"
        return BilevelProblem(
            "DempeEtal2012",
            [1, 1, 2, 2],
            [1.0, 1.0],
            [-1.0, -1.0, 1.0],
            (x, y) -> x[1],
            (x, y) -> x[1]*y[1],
            (x, y) -> [
                -1 - x[1],
                x[1] - 1
            ],
            (x, y) -> [
                -y[1],
                y[1] - 1
            ],
            [-1.0, 1.0]
        )

    elseif prob_no == 24 || prob_no == "DempeFranke2011Ex41"
        return BilevelProblem(
            "DempeFranke2011Ex41",
            [2, 2, 4, 4],
            [1.0, 1.0, 1.0, 1.0],
            [5.0, -2.0, 1.0],
            (x, y) -> x[1] + y[1]^2 + y[2]^2,
            (x, y) -> x[1]*y[1] + x[2]*y[2],
            (x, y) -> [
                -1 - x[1],
                -1 + x[1],
                -1 - x[2],
                1 + x[2]
            ],
            (x, y) -> [
                -2*y[1] + y[2],
                y[1] - 2,
                y[2] - 2,
                -y[2]
            ],
            [0.0, -1.0, 1.0, 2.0]
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