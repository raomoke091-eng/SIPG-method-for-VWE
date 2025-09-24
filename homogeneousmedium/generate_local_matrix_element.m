%生成积分区域为单元对应的单位刚度矩阵(双线性形式)
function r=generate_local_matrix_element(func_term,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type)


r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);
%创建稀疏方阵 r 
%行数/列数=三角划分中的三角形个数*每个三角形的自由度=总自由度
%解释每个三角形的自由度数=(basis_type+1)*(basis_type+2)/2=T_n：
%basis_type指二元(x,y)多项式的次数/最高项次数，这样的多项式有T_n个自由度/待定系数
%对于Lagrange型形状函数，每个节点上只有一个自由度(每个节点处只定义一个形状函数)
%因此n次形状函数需要T_n个节点与之对应，才能建立T_n个方程来求解T_n个系数
%T_n=二元n次多项式的自由度数=每个单元/元素/element上的节点个数(Lagrange型)=每个多项式对应求解的方程个数=每个单元上的Lagrange型形状函数个数

for n=1:size(DGT,2) %循环迭代处理每个三角形
    vertices=DGM(:,DGT(:,n));
    %获取第n个三角形的三个顶点坐标，存储在2*3大小的矩阵 vertices 中
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    %生成第n个三角形的local高斯积分权重和坐标
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    %创建一个大小为 (basis_type+1)*(basis_type+2)/2 的全零方阵 b，用于存储单元刚度矩阵
    %单元刚度矩阵的大小：单元上的每个形状函数都要两两作用(作用方式为双线性形式)
    %参考李开泰《有限元方法及其应用》P16公式

    %以下循环用于生成单元刚度矩阵 n个三角形单元对应n个单元刚度矩阵
    for k=1:length(Gauss_coefficient_local_triangle)
        a=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,0,0); 
        %compute values of basis function at local Gaussian points  输入函数triangular_local_basis的点坐标(x,y)为当前遍历local高斯点坐标

        c=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,1,0)'; 
        d=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,0,1)';
        e=[c d];
        %compute derivatives of basis function at local Gaussian points
        for i=1:(basis_type+1)*(basis_type+2)/2
            for j=1:(basis_type+1)*(basis_type+2)/2
                b(i,j)=b(i,j)+Gauss_coefficient_local_triangle(k)*feval(func_term,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2))*dot(e(i,:),e(j,:));
            end
        end
    end

    for i=1:(basis_type+1)*(basis_type+2)/2
        for j=1:(basis_type+1)*(basis_type+2)/2
            r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b(i,j);
        end
    end
end