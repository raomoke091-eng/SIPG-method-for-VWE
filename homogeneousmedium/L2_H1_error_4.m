function temp=L2_H1_error_4(DGM,DGT,X,func_name_1,func_name_2,func_name_3,func_name_4,t,basis_type, ...
                            derivative_order_x_1,derivative_order_y_1,derivative_order_x_2,derivative_order_y_2,derivative_order_x_3,derivative_order_y_3,derivative_order_x_4,derivative_order_y_4, ...
                            Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)
%参数t是在误差中要保留的部分（积分外的部分）
temp=0; %返回值 初始化temp用于存储误差的平方和
for n=1:size(DGT,2) %遍历每个三角单元
    temp1=0;temp2=0;temp3=0;temp4=0;
    uh_local=X((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1); %X=r_2
    %(n-1)*(basis_type+1)*(basis_type+2)/2+1为当前三角形数值解在 X 向量中的起始索引，n*(basis_type+1)*(basis_type+2)/2为当前三角形数值解在 X 向量中的结束索引
    %数值解是以基函数的系数向量的形式保存的，因此下面计算误差的公式中要把系数和基函数放在一起作线性组合展开
    vertices=DGM(:,DGT(:,n)); %2*3
    h_n=element_diameter(vertices);
    %DGT(:,n)当前三角单元的三个顶点编码（以列向量呈现），vertices提取并存储了当前三角单元的顶点坐标 2*3矩阵
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    %生成当前三角单元上的高斯积分点和系数
  
    for k=1:length(Gauss_coefficient_local_triangle)
        temp1=temp1+Gauss_coefficient_local_triangle(k)*(feval(func_name_1,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),t)...
             -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),uh_local,derivative_order_x_1,derivative_order_y_1,basis_type))^2;
        temp2=temp2+Gauss_coefficient_local_triangle(k)*(feval(func_name_2,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),t)...
             -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),uh_local,derivative_order_x_2,derivative_order_y_2,basis_type))^2;
        temp3=temp1+Gauss_coefficient_local_triangle(k)*(feval(func_name_3,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),t)...
             -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),uh_local,derivative_order_x_3,derivative_order_y_3,basis_type))^2;
        temp4=temp1+Gauss_coefficient_local_triangle(k)*(feval(func_name_4,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),t)...
             -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),uh_local,derivative_order_x_4,derivative_order_y_4,basis_type))^2;
    end  
    temp=temp+h_n^2*(temp1+temp2+temp3+temp4);
end
temp=sqrt(temp);