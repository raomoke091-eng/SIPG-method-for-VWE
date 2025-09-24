function r=generate_local_matrix_edge(func_term,DGE,DGT,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,beta,sigma,interior_edges,boundary_edges)
%integrand_func_name_1,integrand_func_name_2分别用于内部边和边界边

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

%% Interior edges
for i=1:length(interior_edges)
    n=interior_edges(i);%从数组 interior_edges 中获取索引 i 处的元素，并将其赋值给变量 n，n 代表当前处理的内部边的索引(进行这一步操作是因为interior_edges中内部边的排列可能并不是按顺序排列的)
    ele=DGE(5,n);%从 DGE 中获取当前边所属的三角形元素索引记为 ele
    %ele_neighbor=DGE(5,DGE(6,n));%从 DGE 中获取公用该边的三角形元素的索引记为 ele_neighbor
    neighbor_index=DGE(6,n);
    if mod(neighbor_index,3)~=0
        ele_neighbor=(neighbor_index-mod(neighbor_index,3))/3+1;
    else
        ele_neighbor=neighbor_index/3;
    end
    % 输出调试信息
    %fprintf('n = %d, neighbor_index = %d\n', n, neighbor_index);
    
    begin_point=DGE(1:2,n);%提取内部边 n 的起点坐标
    end_point=DGE(3:4,n);%提取内部边 n 的终点坐标
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);%计算内部边 n 的长度 即公式中的|e|，用于加罚项系数
    normal_vector=DG_edge_flag(2:3,n);
    

    b11=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b12=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b21=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b22=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    %b_2=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    %以上两行初始化了两个单元刚度矩阵 b 和 b_2，它们的大小由形函数的阶数 basis_type 决定，即单个三角形单元上的形函数总数
    if begin_point(1)==end_point(1) %vertical edge竖直边(检查边的起点和终点的 x 坐标是否相等，如果相等，则说明这是一条竖直边)
       %对于竖直边，可以将其上的积分视作y维度上的一元积分
        lower_bound=min(begin_point(2),end_point(2));%定义积分下界
        upper_bound=max(begin_point(2),end_point(2));%定义积分上界
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        %initialize the quadrature points and weights
        
        for k=1:length(Gauss_coefficient_local_1D)
            a1=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0);  %compute values of basis functions on the first face neighbor at Gaussian points
            a2=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0);  %compute values of basis functions on the second face neighbor at Gaussian points
            c1=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,1,0)';
            c2=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,1,0)';
            d1=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,1)';
            d2=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,1)';
            e1=[c1 d1];%compute derivatives of basis functions on the first face neighbor at Gaussian points
            e2=[c2 d2];%compute derivatives of basis functions on the second face neighbor at Gaussian points
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b11(alpha,gamma)=b11(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b12(alpha,gamma)=b12(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)-1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b21(alpha,gamma)=b21(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a1(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b22(alpha,gamma)=b22(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a2(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
        end

    elseif begin_point(2)==end_point(2) % horizontal edge水平边(检查边的起点和终点的 y 坐标是否相等，如果相等，则说明这是一条水平边)
           %对于水平边，可以将其上的积分视作x维度上的一元积分
        lower_bound=min(begin_point(1),end_point(1));%定义积分下界
        upper_bound=max(begin_point(1),end_point(1));%定义积分上界
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);%由标准元上的高斯积分点和积分系数导出一般元上的高斯积分点和积分系数
        for k=1:length(Gauss_coefficient_local_1D)
            a1=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0);  %compute values of basis functions on the first face neighbor at Gaussian points
            a2=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0);  %compute values of basis functions on the second face neighbor at Gaussian points
            c1=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,1,0)';
            c2=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,1,0)';
            d1=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,1)';
            d2=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,1)';
            e1=[c1 d1];%compute derivatives of basis functions on the first face neighbor at Gaussian points
            e2=[c2 d2];%compute derivatives of basis functions on the second face neighbor at Gaussian points
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b11(alpha,gamma)=b11(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b12(alpha,gamma)=b12(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)-1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b21(alpha,gamma)=b21(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a1(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b22(alpha,gamma)=b22(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a2(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
        end

    else %一般内部边
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        slope=(end_point(2)-begin_point(2))/(end_point(1)-begin_point(1));
        Jacobi=sqrt(1+slope^2);

        for k=1:length(Gauss_coefficient_local_1D)
            x=Gauss_point_local_1D(k);
            y=slope*(x-begin_point(1))+begin_point(2);
            a1=triangular_local_basis(x,y,basis_type,0,0);  %compute values of basis functions on the first face neighbor at Gaussian points
            a2=triangular_local_basis(x,y,basis_type,0,0);  %compute values of basis functions on the second face neighbor at Gaussian points
            c1=triangular_local_basis(x,y,basis_type,1,0)';
            c2=triangular_local_basis(x,y,basis_type,1,0)';
            d1=triangular_local_basis(x,y,basis_type,0,1)';
            d2=triangular_local_basis(x,y,basis_type,0,1)';
            e1=[c1 d1];%compute derivatives of basis functions on the first face neighbor at Gaussian points
            e2=[c2 d2];%compute derivatives of basis functions on the second face neighbor at Gaussian points
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b11(alpha,gamma)=b11(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b11(alpha,gamma)=b11(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b12(alpha,gamma)=b12(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,gamma)*dot(e1(alpha,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)+eipsilon*0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b12(alpha,gamma)=b12(alpha,gamma)-1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b21(alpha,gamma)=b21(alpha,gamma)+0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a1(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,alpha)*dot(e1(gamma,:),normal_vector);
                    b21(alpha,gamma)=b21(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,alpha)*a1(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b22(alpha,gamma)=b22(alpha,gamma)-0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,gamma)*dot(e2(alpha,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)-eipsilon*0.5*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,alpha)*dot(e2(gamma,:),normal_vector);
                    b22(alpha,gamma)=b22(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*Jacobi*feval(func_term,x,y)*a2(:,alpha)*a2(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
        end
    end

    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for gamma=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)+b11(alpha,gamma);
        end
    end
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for gamma=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+gamma)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+gamma)+b12(alpha,gamma);
        end
    end
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for gamma=1:(basis_type+1)*(basis_type+2)/2
            r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)...
                =r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)+b21(alpha,gamma);
        end
    end
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for gamma=1:(basis_type+1)*(basis_type+2)/2
            r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+gamma)...
                =r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+gamma)+b22(alpha,gamma);
        end
    end
end


%% Boundary edges
for i=1:length(boundary_edges)
    n=boundary_edges(i);
    ele=DGE(5,n);
    begin_point=DGE(1:2,n);
    end_point=DGE(3:4,n);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            a=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0); 
            c=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,1,0)';
            d=triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,1)';
            e=[c d];
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b(alpha,gamma)=b(alpha,gamma)+Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a(:,gamma)*dot(e(alpha,:),normal_vector);
                    b(alpha,gamma)=b(alpha,gamma)+eipsilon*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a(:,alpha)*dot(e(gamma,:),normal_vector);
                    b(alpha,gamma)=b(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,begin_point(1),Gauss_point_local_1D(k))*a(:,alpha)*a(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            a=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0); 
            c=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,1,0)';
            d=triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,1)';
            e=[c d];
            for alpha=1:(basis_type+1)*(basis_type+2)/2  
                for gamma=1:(basis_type+1)*(basis_type+2)/2
                    b(alpha,gamma)=b(alpha,gamma)+Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a(:,gamma)*dot(e(alpha,:),normal_vector);
                    b(alpha,gamma)=b(alpha,gamma)+eipsilon*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a(:,alpha)*dot(e(gamma,:),normal_vector);
                    b(alpha,gamma)=b(alpha,gamma)+1/edge^beta*sigma*Gauss_coefficient_local_1D(k)*feval(func_term,Gauss_point_local_1D(k),begin_point(2))*a(:,alpha)*a(:,gamma)*dot(normal_vector,normal_vector);
                end
            end
        end
    end
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for gamma=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+alpha,(ele-1)*(basis_type+1)*(basis_type+2)/2+gamma)+b(alpha,gamma);
        end
    end 
end