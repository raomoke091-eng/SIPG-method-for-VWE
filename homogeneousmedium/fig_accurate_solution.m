function fig_accurate_solution(t)
    t=0.008;
    Iter_refinement=5;
    DGM=[0,0;1/8,0;0,1/8;1/8,1/8;0,1/8;1/8,0]'; 
    DGT=[1,2,3;4,5,6]';             
    for i=1:Iter_refinement      
        [DGM,DGT]=uniformrefine_triangle(DGM,DGT);                                                   
    end
    
    % 生成示例网格数据
    node = DGM'; % 节点坐标
    element = DGT'; % 单元节点索引

    % 绘制网格图形
    figure(1);
    trimesh(element, node(:,1), node(:,2), 'LineWidth', 1.5); % 使用 trimesh 函数绘制三角形网格
    xlabel('X'); ylabel('Y'); % 添加坐标轴标签
    title('DG Mesh'); % 添加标题
    axis equal; % 设置坐标轴比例相等
 

    if t==0.002
        U=zeros(3*size(DGT,2),1);
        for n=1:size(DGT,2)
            vertices=DGM(:,DGT(:,n));%2*3
            U(3*n-2,1)=accurate_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),t);
            U(3*n-1,1)=accurate_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),t);
            U(3*n,1)=accurate_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),t);
            %分别计算三角单元的三个顶点处的数值解
        end
        figure(2);
        trisurf(DGT',DGM(1,:)',DGM(2,:)',U);  %函数trisurf的用法
        shading interp  %使用插值方式对图形进行着色
        xlabel('x'), ylabel('y')  %设置坐标轴标签
    end

     if t==0.005
        U=zeros(3*size(DGT,2),1);
        for n=1:size(DGT,2)
            vertices=DGM(:,DGT(:,n));%2*3
            U(3*n-2,1)=accurate_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),t);
            U(3*n-1,1)=accurate_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),t);
            U(3*n,1)=accurate_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),t);
            %分别计算三角单元的三个顶点处的数值解
        end
        figure(3);
        trisurf(DGT',DGM(1,:)',DGM(2,:)',U);  %函数trisurf的用法
        shading interp  %使用插值方式对图形进行着色
        xlabel('x'), ylabel('y')  %设置坐标轴标签
     end

     if t==0.008
        U=zeros(3*size(DGT,2),1);
        for n=1:size(DGT,2)
            vertices=DGM(:,DGT(:,n));%2*3
            U(3*n-2,1)=accurate_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),t);
            U(3*n-1,1)=accurate_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),t);
            U(3*n,1)=accurate_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),t);
            %分别计算三角单元的三个顶点处的数值解
        end
        figure(4);
        trisurf(DGT',DGM(1,:)',DGM(2,:)',U);  %函数trisurf的用法
        shading interp  %使用插值方式对图形进行着色
        xlabel('x'), ylabel('y')  %设置坐标轴标签
    end
