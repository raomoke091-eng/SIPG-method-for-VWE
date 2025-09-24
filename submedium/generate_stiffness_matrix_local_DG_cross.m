function r=generate_stiffness_matrix_local_DG_cross(func_right_term,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,derivative_degree_x,derivative_degree_y)


r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

for n=1:size(DGT,2)
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    for k=1:length(Gauss_coefficient_local_triangle)
        a=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type...
                                ,derivative_degree_x,derivative_degree_y);
        c=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type...
                                ,derivative_degree_y,derivative_degree_x);                    
        b=b+Gauss_coefficient_local_triangle(k)*feval(func_right_term,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2))...
            *a'*c;
    end
    for i=1:(basis_type+1)*(basis_type+2)/2
        for j=1:(basis_type+1)*(basis_type+2)/2
            r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b(i,j);
        end
    end
end


        