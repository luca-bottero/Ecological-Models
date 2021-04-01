using DifferentialEquations, Plots

function plant_prey_predator(du,u,p,t)
    #x, y, z = u
   # dx, dy, dz = du
    x, y, z = u

    du[1] = -p[3]*x     #plant spontaneous death rate
    du[1] += -p[4]*y    #plants eaten by prey
    if x >= p[2]        #plant reproduction threshold
        du[1] += p[8]*x*(1 - x/p[1])     #plant reproduction with saturation
        du[1] += p[5]*x*(1 + p[6]*cos(2*pi/p[7]*t)^2)   #stagional growth factor
    end

    du[2] = -p[9]*y     #prey spontaneous death rate
    if y >= p[10]
        du[2] += p[11]/(1 + exp(y))    #prey reproduction
        du[2] += p[12]*x*y  #prey feeding growth factor
    end
    du[2] += -p[16]*y*z     #prey death due to predation
        
    du[3] = -p[13]*z    #predator spontaneous death rate
    if z >= p[14]
        du[3] += p[15]/(1 + exp(z))   #predator reproduction
        du[3] += p[16]*y*z  #predator feeding growth factor
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
10  Prey reproduction threshold
11  Prey reproduction speed
12  Prey feeding efficiency


13  Predator spontaneous death rate
14  Predator reproduction threshold
15  Predator reproduction speed
16  Predator feeding efficiency

=#

u0 = [100.,10., 4.]

p_plant = [10000., 0., 0.1, 0.1, 1.5, 1.2, 10., 3.8]
p_prey = [0.2, 0., 0., 1e-4]
p_predator = [0.5, 0., 0.0005, 1e-5]

p = vcat(vcat(p_plant, p_prey), p_predator)

#print(p)

tspan = (0.,10.)

prob = ODEProblem(plant_prey_predator,u0,tspan,p)

@time sol = solve(prob)

time_plot = plot(sol)
xyz_plot = plot(sol, vars = (1,2,3), xlabel = "Plant", ylabel="Prey", zlabel = "Predator")

plot(time_plot, xyz_plot)
