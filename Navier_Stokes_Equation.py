
# ASSUMPTIONS
# - Incompressible flow
# - Isothermal
# - No buoyancy
# - Constant fluid properties
# - Absolute pressure is irrelevant (only pressure gradients matter)

# FLUID PROPERTIES (choose based on environment)
# have altitude &temperature effects implicitly
    # kg/m^3  (density)
    # m^2/s   (kinematic viscosity)


#u[:, :] = 0.0      # initially at rest
#v[:, :] = 0.0      # initially at rest
#p[:, :] = 0.0      # reference pressure (relative)



import numpy as np
import ezdxf as DXF
import matplotlib.pyplot as plt


# UPSCALED GRID
small_testcase = np.array([[1,1,1,1,1,1,1,1],
                           [1,0,0,0,0,0,0,1],
                           [1,0,0,0,0,0,0,1],
                           [1,0,0,0,0,0,0,1],
                           [1,1,1,1,1,1,1,1]])

scale_y = 10 # 50 / 5
scale_x = 10 # 80 / 8

# repeat each row/column to make bigger grid
grid = np.kron(small_testcase, np.ones((scale_y, scale_x))) #kronecker
ny, nx = grid.shape

# initialize fields
u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.ones((ny,nx)) *1e5

fx = 1.0 # m/s^2 external driving force
fy = 2.0

density = 1000.0
viscocity  = 1e-6  #slightly high(kind of) for smooth sims
g = 9.8

dx = 0.01
dy = 0.01
dt = 0.00001
nt = 200   #timestamps


# DIVERGENCE AND LAPLACIAN FUNCTIONS

def divergence(u,v,dx,dy): #This computes ∇.u(used in poisson pressure)
    div = np.zeros_like(u)
    div[1:-1,1:-1] = ((u[1:-1,2:]-u[1:-1,:-2])/(2*dx)) + (v[2:, 1:-1] - v[:-2, 1:-1]) / (2*dy) # 2*dx and 2*dy central difference
    return div

def laplacian(f,dx,dy): #viscous diffusion
    lap = np.zeros_like(f)
    lap[1:-1, 1:-1] = ((f[1:-1, 2:] - 2*f[1:-1, 1:-1] + f[1:-1, :-2]) / dx**2 +
                        (f[2:, 1:-1] - 2*f[1:-1, 1:-1] + f[:-2, 1:-1]) / dy**2)
    return lap


# PRESSURE POISSON SOLVER

def pressure_poisson(p,u_star,v_star,dx,dy,density,dt,iters = 2000):
    #gauss-seidel
    for _ in range(iters):
        p[1:-1,1:-1] = (((p[1:-1, 2:] + p[1:-1, :-2]) * dy**2 +
                          (p[2:, 1:-1] + p[:-2, 1:-1]) * dx**2 -
                          density * dx**2 * dy**2 / dt * divergence(u_star, v_star, dx, dy)[1:-1, 1:-1])
                         / (2*(dx**2 + dy**2)))
    return p 


# TIME STEPPING LOOP

for n in range(nt):
    #velocity prediction(without pressure)
    u_star = u + dt * (-u*np.gradient(u,dx,axis = 1)- v*np.gradient(u,dy,axis =0) +
                       viscocity * laplacian(u, dx, dy) + fx)
    v_star = v + dt * (- u * np.gradient(v, dx, axis=1)- v * np.gradient(v, dy, axis=0)  +
                       viscocity* laplacian(v, dx, dy)  + fy)
    
    #solve pressure poisson
    p = pressure_poisson(p, u_star, v_star, dx, dy, density, dt, iters=200)
    
    #velocity correction (removes divergence with poisson equation)
    u[1:-1, 1:-1] = u_star[1:-1, 1:-1] - dt/density * ((p[1:-1, 2:] - p[1:-1, :-2]) / (2*dx))
    v[1:-1, 1:-1] = v_star[1:-1, 1:-1] - dt/density * ((p[2:, 1:-1] - p[:-2, 1:-1]) / (2*dy))
    
    #boundary walls (enforce no-slip)
    u[grid == 1] = 0 
    v[grid ==1] = 0


# MATPLOTLIB CHECKING
#defining axis
x = np.linspace(0, (nx-1)*dx, nx)
y = np.linspace(0, (ny-1)*dy, ny)
X, Y = np.meshgrid(x, y)

#defining velocity vectors(quiver plot)
plt.figure()
plt.quiver(X[::2, ::2], Y[::2, ::2],
           u[::2, ::2], v[::2, ::2])
plt.xlabel("x")
plt.ylabel("y")
plt.title("Velocity field")
plt.axis("equal")

#defining velocity magnitude
speed = np.sqrt(u**2 + v**2)
plt.figure()
plt.imshow(speed, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()]) #velocity gradient
plt.colorbar(label="Speed")
plt.title("Velocity magnitude")
plt.xlabel("x")
plt.ylabel("y")

#defining pressure field
plt.figure()
plt.imshow(p, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()]) #pressure gradient
plt.colorbar(label="Pressure")
plt.title("Pressure field")
plt.xlabel("x")
plt.ylabel("y")


 
plt.show()
