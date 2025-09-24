function r=generate_stiffness_matrix_local_DG_edge(integrand_func_name_1,integrand_func_name_2,DGM,DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,penalty_1,penalty_2,co_beta,interior_edges,boundary_edges)
%integrand_func_name_1,integrand_func_name_2分别用于内部边和边界边

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

%% Interior edges
for i=1:length(interior_edges)
    n=interior_edges(i);%从数组 interior_edges 中获取索引 i 处的元素，并将其赋值给变量 n，n 代表当前处理的内部边的索引(进行这一步操作是因为interior_edges中内部边的排列可能并不是按顺序排列的)
    ele=DGE(5,n);%从 DGE 中获取当前边所属的三角形元素索引记为 ele
    ele_neighbor=DGE(5,DGE(6,n));%从 DGE 中获取公用改边的三角形元素的索引记为 ele_neighbor
    begin_point=DGE(1:2,n);%提取内部边 n 的起点坐标
    end_point=DGE(3:4,n);%提取内部边 n 的终点坐标
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);%计算内部边 n 的长度 即公式中的|e|，用于加罚项系数
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b_2=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    %以上两行初始化了两个单元刚度矩阵 b 和 b_2，它们的大小由形函数的阶数 basis_type 决定，即单个三角形单元上的形函数总数
    if begin_point(1)==end_point(1) % vertical edge竖直边(检查边的起点和终点的 x 坐标是否相等，如果相等，则说明这是一条竖直边)
                                    %对于竖直边，可以将其上的积分视作y维度上的一元积分
        lower_bound=min(begin_point(2),end_point(2));%定义积分下界
        upper_bound=max(begin_point(2),end_point(2));%定义积分上界
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge水平边(检查边的起点和终点的 y 坐标是否相等，如果相等，则说明这是一条水平边)
                                        %对于水平边，可以将其上的积分视作x维度上的一元积分
        lower_bound=min(begin_point(1),end_point(1));%定义积分下界
        upper_bound=max(begin_point(1),end_point(1));%定义积分上界
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);%由标准元上的高斯积分点和积分系数导出一般元上的高斯积分点和积分系数
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
            %尽管b和b_2两个矩阵是同时生成的，但要注意的是每一条边都是要遍历的
            %但在这种编码方式之下，每条边实际上要被计算量词，因此需要做一个拆分
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
            b=b+Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name_1,x,y...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name_1,x,y...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
        end
    end
    
    %dg方法由单元刚度矩阵组合成总刚度矩阵的方式和fem不同
    %dg方法对基函数的连续性没有要求，因此各个单元上的形函数实际上就是整个计算域上的基函数，因此有如下的组合方式
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for beta=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b(alpha,beta);
            r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b_2(alpha,beta);
        end
    end
end


%% Boundary edges
for i=1:length(boundary_edges)
    n=boundary_edges(i);ele=DGE(5,n);
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_2,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty_2,DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_2,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty_2,DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
        end
    end
    
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for beta=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b(alpha,beta);
        end
    end
end
