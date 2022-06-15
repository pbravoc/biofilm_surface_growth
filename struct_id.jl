using StructuralIdentifiability

ode = @ODEmodel(
	x1'(t) = K-k*x1(t),
	y(t) = x1(t)
)
##
global_id = assess_identifiability(ode, 0.99)
##