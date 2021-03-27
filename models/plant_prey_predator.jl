using DifferentialEquations

function plant_prey_predator(du,u,p,t)
    #x, y, z = u
   # dx, dy, dz = du
    x = u[1]

    du[1] = -0.9*x

    if x >= p[2]
        du[1] += x*(1 - x/p[1])
    end
end

u0 = [6.]
p = [100.,5.]
tspan = (0.,10.)

prob = ODEProblem(plant_prey_predator,u0,tspan,p)
sol = solve(prob)

plot(sol)
