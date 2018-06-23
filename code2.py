import fenics as fe
import matplotlib.pyplot as plt
import numpy as np
Mesh_size=100
L=1.01/0.045
mesh = fe.UnitIntervalMesh(Mesh_size)
P1 = fe.FiniteElement('P', mesh.ufl_cell(), 1)
V = fe.FunctionSpace(mesh, P1)
v = fe.TestFunction(V)
h = fe.Function(V)
ste=fe.Constant(0.045)
prandl_number=fe.Constant(1)
rho=1
K=1
c=1
def phi(T):
    r=0.025
    T_r=0
    return 0.5*(1. + fe.tanh((T_r - T)/r))
def phis(h):
	r=0.025
	return 0.5*(1. + fe.tanh((-h)/r))
def T(h):
	return c*(h**3)
T=T(h)
phi_h=phis(h)
hot_wall_enthalpy = 0
h_h = fe.Constant(hot_wall_enthalpy)
cold_wall_enthalpy = c
h_c = fe.Constant(cold_wall_enthalpy)
initial_melt_thickness = 10./float(Mesh_size)
h_n = fe.interpolate(
    fe.Expression(
        "(h_h - h_c)*(x[0] < x_m0) + h_c",
        h_h = hot_wall_enthalpy, 
        h_c = cold_wall_enthalpy,
        x_m0 = initial_melt_thickness,
        element = P1),
    V)

timestep_size = 1.e-2
Delta_t = fe.Constant(timestep_size)
T_t = (T )/Delta_t
phi_t = (phis(h) - phis(h))/Delta_t
h_t=(h-h_n)/Delta_t
diffth=3*c*(h**2)
F=(diffth*(fe.dot(fe.grad(v),fe.grad(h)))+(rho/K)*(v*h_t))*fe.dx
JF = fe.derivative(F, h, fe.TrialFunction(V))
hot_wall = "near(x[0],  0)"
cold_wall = "near(x[0],  1.0)"
hot_enthalpy=0
cold_enthalpy=c
boundary_conditions = [
fe.DirichletBC(V, hot_enthalpy, hot_wall),
fe.DirichletBC(V, cold_enthalpy, cold_wall)]
problem = fe.NonlinearVariationalProblem(F, h, boundary_conditions, JF)
solver = fe.NonlinearVariationalSolver(problem)
solver.solve()
print(h)
fe.plot(h)
plt.xlabel('x')
plt.ylabel('Enthalpy')
plt.title('Enthaply after solving using newton''s method')
plt.show()

