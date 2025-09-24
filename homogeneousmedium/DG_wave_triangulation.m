function [L2_error,H1_error,DG_broken_error,r_2,CPU_time]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T,constant,coefficient)
% 04/26/2024
% eipsilon is -1 for SIPG; 
% the coefficients of penalty term are sigma/|e|^beta
% viscoelastic wave equation is dz/dt - div(grad_z) = f, (0,T)*[0,1]^2, let T be 1
%                     with boundary condition g_D and initial condition z_0
% Solve the matrix system: (M + theta*delta_t*A) x_{n+1} = (M + (1-theta)*delta_t*A) x_{n} + delta_t * F(theta*t_{n+1} + (1-theta)*t_{n})
% M is the mass matrix, A is the stiffness matrix

if eipsilon==1 % NIPG
    sigma_1=sigma;sigma_2=2*sigma;
else
    sigma_1=sigma;sigma_2=sigma;
end

%% Using uniform triangular partition, initial triangulation
%% 均匀三角划分创建初始网格（包含两个三角形）（二维情形）
DGM=[0,0;1,0;0,1;1,1;0,1;1,0]'; %DGM矩阵包含网格节点坐标（每一列表示一个节点，每一行代表一个维度）
DGT=[1,2,3;4,5,6]';             %DGT矩阵包含三角形的顶点编码/索引  对应DGM的节点坐标
                                %都是把个数放在列的维度（DGM的列是节点个数，DGT的列是三角形个数）


%% Begin uniform refinement
for i=1:Iter_refinement      %Iter_refinement迭代次数
    [DGM,DGT]=uniformrefine_triangle(DGM,DGT);    %调用自定义函数uniformrefine_triangle(DGM,DGT)对初始DGM、DGT进行细化（分割成更小的三角形，增加网格分辨率）
                                                  %每for循环一次，即每迭代一次，网格细化一次（具体为分割成4个子三角形）
end

%接下来的操作在完成迭代细化之后
[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT); 
[un_used,interior_edges]=find(DGE(6,:));
%通过查找 DGE 矩阵的第 6 行，找到内部边的索引存储在矩阵 interior_edges 中
%find 函数返回DGE(6,:)中非零元素的行索引 un_used 和列索引 interior_edges
boundary_edges=setdiff(1:size(DGE,2),interior_edges);
%1:size(DGE,2) 生成了一个从 1 到 DGE 矩阵列数的向量，表示所有边的索引
%通过 setdiff 函数从中排除内部边的索引，剩下即为边界边索引，这些索引保存在 boundary_edges 变量中
%[rows, cols] = size(DGE);
%for i = 1:rows
%    for j = 1:cols
%        fprintf('%d\t', DGE(i, j));
%    end
%    fprintf('\n');
%end


[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(12);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(12);


%% Assemble stiffness matrix组装刚度矩阵
M=generate_stiffness_matrix_local_DG('mass_constant',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,0);                                        
A1a=generate_local_matrix_element('damping_coefficient',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);

A2a=generate_local_matrix_edge('damping_coefficient',DGE,DGT,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                             ,basis_type,eipsilon,beta,sigma,interior_edges,boundary_edges);                                        
A1b=generate_local_matrix_element('diffusion_coefficient',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);

A2b=generate_local_matrix_edge('diffusion_coefficient',DGE,DGT,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                             ,basis_type,eipsilon,beta,sigma,interior_edges,boundary_edges);                                      
                                             
Aa=A1a+A2a;
Ab=A1b+A2b;
clear A1a A2a A1b A2b 


%% LOOP——initial_steps:

% Because the load vector depends on the time, we must write it into the time cycle
delta_t=T/N_T;  

b_0=generate_load_vector_local_DG('accurate_solution',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_0=M\b_0;%%u_0为u(0)的L2投影

b_1=generate_load_vector_local_DG('accurate_solution_time_derivate',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
v_0=M\b_1;%%v_0为u'(0)的L2投影                                

B_0=generate_load_vector_local_DG('function_f',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_tilde=M\(B_0-Aa*v_0-Ab*u_0);
% b_2=generate_load_vector_local_DG('accurate_solution_time_derivate_2',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
% u_tilde=M\b_2;%%也可以选为u''(0)的L2投影

r_1=u_0+delta_t*v_0+(delta_t)^2/2*u_tilde; %t=1时刻U1
r_0=u_0; %t=0时刻U0

%% LOOP
CPUTime=zeros(1,N_T-1);
for n=0:N_T-2
    tic %计时
    t_0=n*delta_t; 
    t_1=(n+1)*delta_t; 
    t_2=(n+2)*delta_t;
    t=t_1;
    fprintf('当前时间为%f/%f.',t_2,T);
    
    b=generate_load_vector_local_DG('function_f',t,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);% Assemble the load vector
    b=(delta_t)^2*b+2*M*r_1-M*r_0+(delta_t)/2*Aa*r_0-(delta_t)^2/2*Ab*r_0;  
    left_matrix=M+delta_t/2*Aa+(delta_t)^2/2*Ab;
    
    % Solve
    r_2=left_matrix\b;
    clear b
    fprintf('加细第%d次的网格,第%d/%d个时间步',Iter_refinement,n+2,N_T);
    CPUTime(1,n+1) = toc;
    disp(['历时CPUTime: ', num2str(CPUTime(1,n+1)), ' 秒']);
    %% 每一步的Error
    %L2 norm error
    L2_error=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution',t_2,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle); 
    
    %H1 norm error
    temp=L2_error;
    temp1=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution_x_derivative',t_2,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    temp2=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution_y_derivative',t_2,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    H1_error=sqrt(temp1^2+temp2^2+temp^2);  
    
    %DG broken norm error
    temp3=L2_H1_error_3(DGM,DGT,r_2,'mass_constant','accurate_solution_x_derivative','accurate_solution_y_derivative',t_2,basis_type, ...
                            1,0,0,1, ...
                            Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    temp4=L2_H1_error_4(DGM,DGT,r_2,'accurate_solution_x_derivative_2','accurate_solution_y_derivative_2','accurate_solution_x_derivative_y_derivative','accurate_solution_y_derivative_x_derivative',t_2,basis_type, ...
                            2,0,0,2,1,1,1,1, ...
                            Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    temp5=L2_H1_error_5(DGT,DGE,DG_edge_flag,r_2,'mass_constant','accurate_solution',t_2,basis_type,0,0,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
    DG_broken_error=sqrt(temp3^2+temp4^2+temp5^2);
    
    fprintf('当前时间的L^2误差为%f，H^1误差为%f，DG离散误差为%f.\n',L2_error,H1_error,DG_broken_error);
    % Evaluate
    r_0=r_1;
    r_1=r_2;

    
    %% Figure作图
    %绘制三角形网格上的数值解 uh 的三维图形
    uh=zeros(3*size(DGT,2),1);
    %创建一个大小为 3*size(DGT,2) 的列向量 uh，用于存储数值解在每个三角单元的三个顶点处的值（每个三角形分配3个维度）

    for n=1:size(DGT,2)
        vertices=DGM(:,DGT(:,n));%2*3
        uh_local=r_2((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);%当前三角单元上的数值解（基函数系数）
        uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
        uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
        uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
        %分别计算三角单元的三个顶点处的数值解
    end

    figure(2024);
    trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
    shading interp  %使用插值方式对图形进行着色
    xlabel('x'), ylabel('y')  %设置坐标轴标签
    
    %分时刻作图  10个时刻+初始时刻
    for i=1:10
        t_pre_index = i*10;
        delta_t = T/N_T;
        t_a = (i-1)*10*(T/N_T);
        t_b = (i+1)*10*(T/N_T);
        if (t>t_a)&&(t<t_b)
            uh=zeros(3*size(DGT,2),1);
            for n=1:size(DGT,2)
                vertices=DGM(:,DGT(:,n));
                uh_local=r_2((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
                uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
                uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(i)
            %先把画图的数组里面的数据对应起来  三角形顶点的坐标和u_h里面的值要对应正确
            %目前u_h已是一列向量  三角形1顶点1，三角形1顶点2，三角形1顶点3，三角形2顶点1，三角形2定点，...（以此类推）
            %把DGM里面的顶点坐标和u_h对应起来 要结合DGM和DGT两个矩阵的信息
            %注意到这种对应关系是这样的：在DGT中的序号是几，这个点在DGM中的坐标就在第几列
            %因此只需要把DGM转置一下就能得到一个第一列为x坐标数据，第二列为y坐标数据的矩阵
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            %colormap(gray);
            view(30, 50);% 设置视角
            lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            ztickformat('%.2e');   %z轴刻度以科学计数法显示
            title(sprintf('3D images of the viscoelastic wave propagation \n in homofeneous medium at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            saveas(gcf, fileName_3d);  % 保存文件
            

            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            title(sprintf('contours of the viscoelastic wave propagation \n in homofeneous medium at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            saveas(gcf, fileName_cont);  % 保存文件
        end
    end
end
                                             
%% 最终时间T时刻的Error
%L2 norm error
L2_error=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution',T,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);

%H1 norm error
temp=L2_error;
temp1=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution_x_derivative',T,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp2=L2_H1_error_1(DGM,DGT,r_2,'accurate_solution_y_derivative',T,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
H1_error=sqrt(temp1^2+temp2^2+temp^2);
                                    
%DG broken norm error
temp3=L2_H1_error_3(DGM,DGT,r_2,'mass_constant','accurate_solution_x_derivative','accurate_solution_y_derivative',T,basis_type, ...
                        1,0,0,1, ...
                        Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp4=L2_H1_error_4(DGM,DGT,r_2,'accurate_solution_x_derivative_2','accurate_solution_y_derivative_2','accurate_solution_x_derivative_y_derivative','accurate_solution_y_derivative_x_derivative',T,basis_type, ...
                        2,0,0,2,1,1,1,1, ...
                        Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp5=L2_H1_error_5(DGT,DGE,DG_edge_flag,r_2,'mass_constant','accurate_solution',T,basis_type,0,0,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
DG_broken_error=sqrt(temp3^2+temp4^2+temp5^2);
% 输出每次网格划分所用的计算时间cpu time
CPU_time=0;
for i=1:size(CPUTime,2)
    CPU_time=CPU_time+CPUTime(1,i);
end

% 生成示例网格数据
node = DGM'; % 节点坐标
element = DGT'; % 单元节点索引

% 绘制网格图形
figure;
trimesh(element, node(:,1), node(:,2), 'LineWidth', 1.5); % 使用 trimesh 函数绘制三角形网格
xlabel('X'); ylabel('Y'); % 添加坐标轴标签
title('DG Mesh'); % 添加标题
axis equal; % 设置坐标轴比例相等
