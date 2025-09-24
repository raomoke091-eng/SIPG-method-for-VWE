function temp=L2_H1_error_5(DGT,DGE,DG_edge_flag,X,func_term_1,func_term_2,t,basis_type,derivative_order_x,derivative_order_y,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

temp=0;
for n=1:size(DGT,2)
    if DGE(6,n)~=0 %interior edges
        ele=DGE(5,n);
        %ele_neighbor=DGE(5,DGE(6,n));
        neighbor_index=DGE(6,n);
        if mod(neighbor_index,3)~=0
            ele_neighbor=(neighbor_index-mod(neighbor_index,3))/3+1;
        else
            ele_neighbor=neighbor_index/3;
        end
        uh_local_1=X((ele-1)*(basis_type+1)*(basis_type+2)/2+1:ele*(basis_type+1)*(basis_type+2)/2,1);
        uh_local_2=X((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+1:ele_neighbor*(basis_type+1)*(basis_type+2)/2,1);
        begin_point=DGE(1:2,n);%提取内部边 n 的起点坐标
        end_point=DGE(3:4,n);%提取内部边 n 的终点坐标
        edge=sqrt((end_point(1)-begin_point(1))^2+(end_point(2)-begin_point(2))^2);
        edge_1=1/edge;
        normal_vector=DG_edge_flag(2:3,n);

        if begin_point(1)==end_point(1) 
            lower_bound=min(begin_point(2),end_point(2));%定义积分下界
            upper_bound=max(begin_point(2),end_point(2));%定义积分上界
            [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
            for k=1:length(Gauss_coefficient_local_1D)
                temp=temp+edge_1*Gauss_coefficient_local_1D(k)*(feval(func_term_1,begin_point(1),Gauss_point_local_1D(k),t))^2* ...
                    (fe_solution(begin_point(1),Gauss_point_local_1D(k),uh_local_2,derivative_order_x,derivative_order_y,basis_type)- ...
                    fe_solution(begin_point(1),Gauss_point_local_1D(k),uh_local_1,derivative_order_x,derivative_order_y,basis_type))^2 ...
                    *dot(normal_vector,normal_vector);
            end
           
        elseif begin_point(2)==end_point(2) 
                lower_bound=min(begin_point(1),end_point(1));%定义积分下界
                upper_bound=max(begin_point(1),end_point(1));%定义积分上界
                [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);%由标准元上的高斯积分点和积分系数导出一般元上的高斯积分点和积分系数
                for k=1:length(Gauss_coefficient_local_1D)
                    temp=temp+edge_1*Gauss_coefficient_local_1D(k)*(feval(func_term_1,Gauss_point_local_1D(k),begin_point(2),t))^2* ...
                    (fe_solution(Gauss_point_local_1D(k),begin_point(2),uh_local_2,derivative_order_x,derivative_order_y,basis_type)- ...
                    fe_solution(Gauss_point_local_1D(k),begin_point(2),uh_local_1,derivative_order_x,derivative_order_y,basis_type))^2 ...
                    *dot(normal_vector,normal_vector);
                end
               
        else
            lower_bound=min(begin_point(1),end_point(1));
            upper_bound=max(begin_point(1),end_point(1));
            [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
            slope=(end_point(2)-begin_point(2))/(end_point(1)-begin_point(1));
            Jacobi=sqrt(1+slope^2);
            for k=1:length(Gauss_coefficient_local_1D)
                x=Gauss_point_local_1D(k);
                y=slope*(x-begin_point(1))+begin_point(2);
                temp=temp+edge_1*Gauss_coefficient_local_1D(k)*Jacobi*(feval(func_term_1,x,y,t))^2* ...
                    (fe_solution(x,y,uh_local_2,derivative_order_x,derivative_order_y,basis_type)- ...
                    fe_solution(x,y,uh_local_1,derivative_order_x,derivative_order_y,basis_type))^2 ...
                    *dot(normal_vector,normal_vector);
            end
         
        end
    else
        ele=DGE(5,n);
        uh_local=X((ele-1)*(basis_type+1)*(basis_type+2)/2+1:ele*(basis_type+1)*(basis_type+2)/2,1);
        begin_point=DGE(1:2,n);
        end_point=DGE(3:4,n);
        edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
        edge_1=1/edge;
        normal_vector=DG_edge_flag(2:3,n);
        if begin_point(1)==end_point(1) %vertical edge
            lower_bound=min(begin_point(2),end_point(2));
            upper_bound=max(begin_point(2),end_point(2));
            [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
            for k=1:length(Gauss_coefficient_local_1D)
                temp=temp+edge_1*Gauss_coefficient_local_1D(k)*(feval(func_term_1,begin_point(1),Gauss_point_local_1D(k),t))^2* ...
                        (feval(func_term_2,begin_point(1),Gauss_point_local_1D(k),t)- ...
                        fe_solution(begin_point(1),Gauss_point_local_1D(k),uh_local,derivative_order_x,derivative_order_y,basis_type))^2 ...
                        *dot(normal_vector,normal_vector);
            end
          
        else  % horizontal edge
            lower_bound=min(begin_point(1),end_point(1));
            upper_bound=max(begin_point(1),end_point(1));
            [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
            for k=1:length(Gauss_coefficient_local_1D)
                temp=temp+edge_1*Gauss_coefficient_local_1D(k)*(feval(func_term_1,Gauss_point_local_1D(k),begin_point(2),t))^2* ...
                        (feval(func_term_2,Gauss_point_local_1D(k),begin_point(2),t)- ...
                        fe_solution(Gauss_point_local_1D(k),begin_point(2),uh_local,derivative_order_x,derivative_order_y,basis_type))^2 ...
                        *dot(normal_vector,normal_vector);
            end
           
        end
        
    end
    
end
temp=sqrt(temp);