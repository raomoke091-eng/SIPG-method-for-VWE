function r=fe_solution(x,y,solution,derivative_order_x,derivative_order_y,basis_type)


r=triangular_local_basis(x,y,basis_type,derivative_order_x,derivative_order_y)*solution;
