import fenics 
import matplotlib.pyplot as plt
def Enthalpyvariationalform(bc):
	rho=1
	k=1
	mesh = UnitSquareMesh(20,20)
	V = FunctionSpace(mesh,"Lagrange",1)
	u0 = Expression("2*pow(x[0],2)",degree=2)
	def u0_boundary(x,on_boundary):
		return on_boundary
	bc = DirichletBC(V,u0,u0_boundary)
	h = TrialFunction(V)
	v = TestFunction(V)
	a =inner(dt*del(h),v)
	L =k/rho*inner(nabla_grad(h),nabla_grad(v))*ds;
	h = Function(V)
	solve(a==L,bc)
	plot(h)
	plt.show()

