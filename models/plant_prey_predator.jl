using DifferentialEquations, Plots

function plant_prey_predator(du,u,p,t)
    #x, y, z = u
   # dx, dy, dz = du
    x, y = u

    du[1] = -p[3]*x     #plant spontaneous death rate
    du[1] += -p[4]*y    #plants eaten by prey
    if x >= p[2]        #plant reproduction threshold
        du[1] += p[8]*x*(1 - x/p[1])     #plant reproduction with saturation
        du[1] += p[5]*(1 + p[6]*cos(2*pi/p[7]*t)^2)   #stagional growth factor
    end

    du[2] = -p[9]*y     #prey spontaneous death rate
    if y >= p[12]
        du[2] += p[10]*y    #prey reproduction
        du[2] += p[11]*x*y  #prey feeding growth factor
    end
        




end

#=
PARAMETERS

1   Saturation value for plant population
2   Plant reproduction threshold
3   Spontaneous death rate parameter
4   Plant eaten parameter
5   Stagional growth factor
6   Stagional variation strength
7   Stagional variation time scale (period)
8   Plant reproduction speed

9   Prey spontaneous death rate
10  Prey reproduction speed
11  Prey feeding efficiency
12  Prey reproduction threshold

=#

u0 = [100.,14.]
p_plant = [1000., 5., 0.1, 0.1, 1.5, 1.2, 4., 1.8]
p_prey = [0.2, 0.01, 0.01, 10.]

p = vcat(p_plant, p_prey)

tspan = (0.,10.)

prob = ODEProblem(plant_prey_predator,u0,tspan,p)

@time sol = solve(prob)

plot(sol)