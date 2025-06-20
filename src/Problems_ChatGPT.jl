model = BilevelModel(
    HiGHS.Optimizer,
    mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = 100, dual_big_M = 100)
)

# Bard1
function bard1()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 3)^2 + (y - 2)^2)
    @objective(Lower(model), Min, (y - x)^2)
    @constraint(Lower(model), y <= 4)
    return model
end
academic_registry["bard1"] = bard1

# Bard2
function bard2()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, x^2 + (y - 1)^2)
    @objective(Lower(model), Min, (y - 1)^2)
    @constraint(Lower(model), y + x >= 1)
    return model
end
academic_registry["bard2"] = bard2

# Bard3
function bard3()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 1)^2 + y^2)
    @objective(Lower(model), Min, (y - 2)^2)
    @constraint(Lower(model), y >= x + 1)
    return model
end
academic_registry["bard3"] = bard3

# ShimizuEtal1997
function shimizuetal1997()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), y[1:2] >= 0)
    @variable(Lower(model), x[1:2] >= 0)
    @objective(Upper(model), Min, (y[1] - 1)^2 + (y[2] - 1)^2 + x[1]^2 + x[2]^2)
    @objective(Lower(model), Min, (y[1] - x[1])^2 + (y[2] - x[2])^2)
    @constraint(Lower(model), y[1] + y[2] <= 3)
    @constraint(Upper(model), x[1] + x[2] <= 4)
    return model
end
academic_registry["shimizuetal1997"] = shimizuetal1997

const tp_registry = Dict{String, Function}()

# TP3
function tp3()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 2)^2 + (y - 3)^2)
    @objective(Lower(model), Min, (y - 2)^2)
    @constraint(Lower(model), y >= 1)
    return model
end
tp_registry["tp3"] = tp3

# TP4
function tp4()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 4)^2 + (y - 2)^2)
    @objective(Lower(model), Min, (y - x)^2)
    @constraint(Lower(model), y <= 3)
    return model
end
tp_registry["tp4"] = tp4

# TP5
function tp5()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 1)^2 + (y - 4)^2)
    @objective(Lower(model), Min, (y - x)^2 + (y - 1)^2)
    @constraint(Lower(model), y >= 1)
    return model
end
tp_registry["tp5"] = tp5

# TP6
function tp6()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 3)^2 + y^2)
    @objective(Lower(model), Min, (y - 2)^2 + (x - 1)^2)
    @constraint(Lower(model), y + x <= 5)
    return model
end
tp_registry["tp6"] = tp6

# TP7
function tp7()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 4)^2 + (y - 3)^2)
    @objective(Lower(model), Min, (y - 1)^2)
    @constraint(Lower(model), y >= 0)
    return model
end
tp_registry["tp7"] = tp7

# TP8
function tp8()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 5)^2 + (y - 2)^2)
    @objective(Lower(model), Min, (y - 3)^2 + (x - 2)^2)
    @constraint(Lower(model), y + x >= 4)
    return model
end
tp_registry["tp8"] = tp8

const inverse_registry = Dict{String, Function}()

# Inverse1
function inverse1()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 1)^2 + (y - 2)^2)
    @objective(Lower(model), Min, (y - x)^2)
    @constraint(Lower(model), y >= 1)
    return model
end
inverse_registry["inverse1"] = inverse1

# Inverse2
function inverse2()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x[1:2] >= 0)
    @variable(Lower(model), y[1:2] >= 0)
    @objective(Upper(model), Min, sum((x[i]-1)^2 for i in 1:2) + sum((y[i]-2)^2 for i in 1:2))
    @objective(Lower(model), Min, sum((y[i]-x[i])^2 for i in 1:2))
    @constraint(Lower(model), sum(y) >= 2)
    return model
end
inverse_registry["inverse2"] = inverse2

# Inverse3
function inverse3()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 2)^2 + (y - 3)^2)
    @objective(Lower(model), Min, (y - x)^2 + (y - 1)^2)
    @constraint(Lower(model), y + x <= 5)
    return model
end
inverse_registry["inverse3"] = inverse3

# Inverse4
function inverse4()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x[1:2] >= 0)
    @variable(Lower(model), y[1:2] >= 0)
    @objective(Upper(model), Min, (x[1]-3)^2 + (x[2]-4)^2 + (y[1]-1)^2 + (y[2]-2)^2)
    @objective(Lower(model), Min, (y[1]-x[1])^2 + (y[2]-x[2])^2)
    @constraint(Lower(model), y[1] + y[2] <= 3)
    return model
end
inverse_registry["inverse4"] = inverse4


const stackelberg_registry = Dict{String, Function}()

# Stack1
function stack1()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 5)^2 + (y - 3)^2)
    @objective(Lower(model), Min, (y - x)^2)
    @constraint(Lower(model), y >= 1)
    return model
end
stackelberg_registry["stack1"] = stack1

# Stack2
function stack2()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x[1:2] >= 0)
    @variable(Lower(model), y[1:2] >= 0)
    @objective(Upper(model), Min, sum((x[i]-1)^2 for i in 1:2) + sum((y[i]-2)^2 for i in 1:2))
    @objective(Lower(model), Min, sum((y[i]-x[i])^2 for i in 1:2))
    @constraint(Lower(model), sum(y) >= 3)
    return model
end
stackelberg_registry["stack2"] = stack2

# Stack3
function stack3()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 2)^2 + (y - 4)^2)
    @objective(Lower(model), Min, (y - 1)^2)
    @constraint(Lower(model), y + x <= 6)
    return model
end
stackelberg_registry["stack3"] = stack3

# Stack4
function stack4()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x[1:2] >= 0)
    @variable(Lower(model), y[1:2] >= 0)
    @objective(Upper(model), Min, (x[1]-3)^2 + (x[2]-5)^2 + (y[1]-2)^2 + (y[2]-1)^2)
    @objective(Lower(model), Min, (y[1]-x[1])^2 + (y[2]-x[2])^2)
    @constraint(Lower(model), y[1] + y[2] >= 4)
    return model
end
stackelberg_registry["stack4"] = stack4

# Stack5
function stack5()
    model = BilevelModel(Ipopt.Optimizer, mode = BilevelJuMP.SOS1Mode())
    @variable(Upper(model), x >= 0)
    @variable(Lower(model), y >= 0)
    @objective(Upper(model), Min, (x - 6)^2 + (y - 3)^2)
    @objective(Lower(model), Min, (y - 2)^2 + (x - 1)^2)
    @constraint(Lower(model), y >= x + 1)
    return model
end
stackelberg_registry["stack5"] = stack5
