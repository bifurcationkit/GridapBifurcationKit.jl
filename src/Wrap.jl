# structure to help casting the functional in a way BifurcationKit can use
struct GridapProblem{Tres,Tjac,Tt_Ω,Td2res,Td3res,Ttri,Tquad,TV,TU, Tls}
	res::Tres		# res(u,p,v), residual
	jac::Tjac		# jac(u,p,du,v), jacobian
	t_Ω::Tt_Ω
	d2res::Td2res	# d2res(u,du1,du2,v)
	d3res::Td3res	# d3res(u,du1,du2,du3,v)
	trian::Ttri
	quad::Tquad
	V::TV
	U::TU
	ls::Tls
end

# constructors
function GridapProblem(res, jac, trian, quad, V, U; autodiff = false, linsolver = nothing)
	return GridapProblem(res, jac, nothing, nothing, nothing, trian, quad, V, U, linsolver)
end

function GridapProblem(res, jac, d2F, d3F, trian, quad, V, U; autodiff = false, linsolver = nothing)
	return GridapProblem(res, jac, nothing, d2F, d3F, trian, quad, V, U, linsolver)
end
################################################################################
function op_from_param(gp::GridapProblem, p)
	res(u,v) = gp.res(u, p, v)
	jac(u,du,v) = gp.jac(u, p, du, v)
	t_Ω = FETerm(res,jac, gp.trian, gp.quad)
	return FEOperator(gp.U, gp.V, t_Ω)
end

# residual
function (gp::GridapProblem)(::Val{:Res}, u, p)
	op = op_from_param(gp, p)
	algop = Gridap.FESpaces.get_algebraic_operator(op)
	return Gridap.FESpaces.residual(algop, u)
end

# (sparse) jacobian matrix
function (gp::GridapProblem)(::Val{:Jac}, u, p)
	op = op_from_param(gp, p)
	algop = Gridap.FESpaces.get_algebraic_operator(op)
	return Gridap.FESpaces.jacobian(algop, u)
end

# # matrix free jacobian
# function (gp::GridapProblem)(u, p, du)
# 	@error "WIP"
# end

# second derivative
function (gp::GridapProblem)(u, p, du1, du2)
	du1h = FEFunction(gp.U, du1)
	du2h = FEFunction(gp.U, du2)
	a(u,v) = gp.d2res(u, p, du1h, du2h, v)
	t_Ω = LinearFETerm(a, gp.trian, gp.quad)
	feop = FEOperator(gp.U, gp.V, t_Ω)
	alop = Gridap.FESpaces.get_algebraic_operator(feop)
	Gridap.FESpaces.residual(alop, u)
end

# third derivative
function (gp::GridapProblem)(u, p, du1, du2, du3)
	du1h = FEFunction(gp.U, du1)
	du2h = FEFunction(gp.U, du2)
	du3h = FEFunction(gp.U, du3)
	a(u,v) = gp.d3res(u, p, du1h, du2h, du3h, v)
	t_Ω = LinearFETerm(a, gp.trian, gp.quad)
	feop = FEOperator(gp.U, gp.V, t_Ω)
	alop = Gridap.FESpaces.get_algebraic_operator(feop)
	Gridap.FESpaces.residual(alop, u)
end

# # user_uh_to_cell_residual
# function (gp::GridapProblem)(::Val{:CellRes}, u, p)
# 	op = op_from_param(gp, p)
# 	algop = Gridap.FESpaces.get_algebraic_operator(op)
# 	return Gridap.FESpaces.residual(algop, u)
# end
