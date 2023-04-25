include("MealModel.jl")
using DifferentialEquations, BenchmarkTools

model = MixedMealModel()

inputs = [
    75000., # Meal Glucose
    75.,    # Body Mass
    60000., # Meal TG
]

model.simulator.initialvalues = [0., 5., 0., 18., 0., 18., 18., 18., 0., 0., 0., 1.3, 0.33]
model.simulator.timespan = (0., 2880.)
#print(model)
solution = @benchmark simulate($model, $inputs)