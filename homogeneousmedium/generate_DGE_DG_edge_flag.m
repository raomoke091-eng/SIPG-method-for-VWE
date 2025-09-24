function [DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT)


DGE=zeros(6,3*size(DGT,2));          %创建一个大小为 6x(3*三角形数目) 的全零矩阵 DGE
DG_edge_flag=zeros(4,3*size(DGT,2)); %创建一个大小为 4x(3*三角形数目) 的全零矩阵 DG_edge_flag

%定义矩阵DGE的1-5行
for i=1:size(DGT,2)
    %3*i-2列到3*i列表示第i个三角形的三条边的参数
    %第1、2行为边起点（编码在前）坐标
    %第3、4行为边终点（编码在后）坐标
    %第5行为三角形编码
    DGE(1:5,3*i-2)=[DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) i]; %每个三角形对应的第一列为该三角形第一个顶点到第二个顶点连边
    DGE(1:5,3*i-1)=[DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) i]; %每个三角形对应的第二列为该三角形第二个顶点到第三个顶点连边
    DGE(1:5,3*i)=[DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) i];   %每个三角形对应的第三列为该三角形第三个顶点到第一个顶点连边
end

edge_middle=[(DGE(1,:)+DGE(3,:))/2;(DGE(2,:)+DGE(4,:))/2];  
%size为2*3*size(DGT,2)的矩阵edge_middle存储每个三角形的三个中点坐标
%存储顺序为12连边中点、23连边中点、31连边中点

%接下来的步骤是完成对网格边界的标记，并计算边的方向信息
for i=1:size(DGE,2)  %for 循环遍历 DGE 矩阵的每一列，即遍历网格的每一条边
    distance=sqrt((edge_middle(1,i)*ones(1,size(edge_middle,2))-edge_middle(1,:)).^2+(edge_middle(2,i)*ones(1,size(edge_middle,2))-edge_middle(2,:)).^2);
    %目的：计算当前遍历到的边的中点到所有其他边中点的距离,存储在size为1*size(edge_middle,2)的矩阵distance中  size(edge_middle,2)=size(DGE,2)=3*size(DGT,2)
    %edge_middle(1,i)当前遍历边（第i边）的中点的x坐标，edge_middle(2,i)当前遍历边（第i边）的中点的y坐标
    %edge_middle(1,:)中点坐标矩阵的第一行-全部边中点的x坐标，edge_middle(2,:)中点坐标矩阵的第二行-全部边中点的y坐标
    [un_used,temp]=find(distance<0.000001); %在 distance 中找到距离小于 0.000001 的元素，并返回它们的行索引un_used和列索引temp
    if length(temp)==2
        %这是判断当前遍历边（第i边）是否为内部边的判据
        %原理：两个相邻的三角形共享一条边，把这条边看作两个三角形中的两条边，则这两条边界的中点会非常接近彼此，用小于0.000001来表达这种"接近"
        %若当前边为内部边，则与它的中点非常接近的边中点有两个：一个是它本身，另一个是与它重合的边的中点
        %因此返回的列索引temp包含本边索引和与其重合边的索引
        j=setdiff(temp,i);   
        %将与当前遍历边重合的边的索引保存在变量 j 中   重合边只有一条
        DGE(6,i)=j;
        ele_i=DGE(5,i);ele_j=DGE(5,j);   
        %矩阵 DGE 的第5行表示当前遍历边所在三角形的编码 
        %获取当前边和对应边所在的三角形索引，分别保存在变量 ele_i 和 ele_j 中
        DG_edge_flag(4,i)=sign(ele_i-ele_j); 
        %根据两个三角形索引的差异，确定边的方向（正负号），并将结果保存在 DG_edge_flag 矩阵的第 4 行中
        DG_edge_flag(2:3,i)=sign(ele_i-ele_j)*[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
        %计算各边的单位法向量，并将其归一化。然后根据两个三角形索引的差异，确定边的法方向，并将结果保存在 DG_edge_flag 矩阵的第 2 行和第 3 行中
        %第2行是边的两点坐标之差的y分量，第3行是边的两点坐标之差的x分量
        %即内部边在矩阵 DG_edge_flag 中的第1行元素为0
    else
        DG_edge_flag(1,i)=1;
        % DG_edge_flag 矩阵的第 1 行的相应元素标记为 1，表示当前边是边界边
        %相应的，DG_edge_flag 矩阵的第 1 行的元素标记为 0，则表示当前边是内部边
        DG_edge_flag(2:3,i)=[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
    end
end

