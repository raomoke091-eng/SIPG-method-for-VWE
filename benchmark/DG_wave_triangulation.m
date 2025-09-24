function [L2_error,H1_error]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T,constant,coefficient)
% 04/26/2016
% eipsilon is -1 for SIPG; 1 for NIPG; 0 for IIPG
% the coefficients of penalty term are sigma/|e|^beta
% Parabolc problem is dz/dt - div(grad_z) = f, (0,T)*[0,1]^2, let T be 1
%                     with boundary condition g_D and initial condition z_0
% Solve the matrix system: (M + theta*delta_t*A) x_{n+1} = (M + (1-theta)*delta_t*A) x_{n} + delta_t * F(theta*t_{n+1} + (1-theta)*t_{n})
% M is the mass matrix, A is the stiffness matrix

if eipsilon==1 % NIPG
    sigma_1=sigma;sigma_2=2*sigma;
else
    sigma_1=sigma;sigma_2=sigma;
end

%% Using uniform triangular partition, initial triangulation
%DGM=[0,0;1,0;0,1;1,1;0,1;1,0]';
DGM=[0,800;270,800;0,960;270,960;0,960;270,800]';
DGT=[1,2,3;4,5,6]';


%% Begin uniform refinement
for i=1:Iter_refinement
    [DGM,DGT]=uniformrefine_triangle(DGM,DGT);
end
 
[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);


[un_used,interior_edges]=find(DGE(6,:));
boundary_edges=setdiff(1:size(DGE,2),interior_edges);

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(12);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(8);


%% Assemble stiffness matrix
M=generate_stiffness_matrix_local_DG('mass_constant',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,0);                                    
A1e=generate_stiffness_matrix_local_DG('diffusive_coefficient_11e',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,1,0);
A1z=generate_stiffness_matrix_local_DG('diffusive_coefficient_11z',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,1,0);                                        
A2e=generate_stiffness_matrix_local_DG('diffusive_coefficient_22e',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);
A2z=generate_stiffness_matrix_local_DG('diffusive_coefficient_22z',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);                                         
A3e=generate_stiffness_matrix_local_DG_cross('diffusive_coefficient_12e',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);  
A3z=generate_stiffness_matrix_local_DG_cross('diffusive_coefficient_12z',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);                                         
A4e=generate_stiffness_matrix_local_DG_cross('diffusive_coefficient_21e',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);
A4z=generate_stiffness_matrix_local_DG_cross('diffusive_coefficient_21z',DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type...
                                             ,0,1);                                          
A5e=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges_e','integrand_2_on_edges_e',DGM,DGT,DGE,DG_edge_flag...
                                                  ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,sigma_1,sigma_2,beta,interior_edges,boundary_edges);
A5z=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges_z','integrand_2_on_edges_z',DGM,DGT,DGE,DG_edge_flag...
                                                  ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,sigma_1,sigma_2,beta,interior_edges,boundary_edges);                                              
Ae=A1e+A2e+A3e+A4e+A5e;
Az=A1z+A2z+A3z+A4z+A5z;
clear A1e A2e A3e A4e A5e A1z A2z A3z A4z A5z


%% LOOP——initial_steps:

% Because the load vector depends on the time, we must write it into the time cycle
%delta_t=T/N_T;
delta_t=8e-5;

b_0=generate_load_vector_local_DG('accurate_solution',0,DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_0=M\b_0;%%u_0为u(0)的L2投影
b_1=generate_load_vector_local_DG('accurate_solution_time_derivate',0,DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
v_0=M\b_1;%%v_0为u'(0)的L2投影                                
B_0=generate_load_vector_local_DG('function_f',0,DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_tilde=M\(B_0-Ae*v_0-Az*u_0);
%u_tilde=M\(B_0-M1*v_0-Ae*v_0-Az*u_0);
% b_2=generate_load_vector_local_DG('accurate_solution_time_derivate_2',0,DGM,DGT...
%                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
% u_tilde=M\b_2;%%也可以选为u''(0)的L2投影
r_1=u_0+delta_t*v_0+(delta_t)^2/2*u_tilde;
r_0=u_0;

%% LOOP
%创建6个数组存储第12和74个三角元的6个顶点的各个时间步长/时刻的数值解的变化情况
n_1 = 169; %12对应34,35,36
n_2 = 1200; %74对应220,221,222
un_1 = zeros((basis_type+1)*(basis_type+2)/2, N_T-1);  %N_T-1次迭代
un_2 = zeros((basis_type+1)*(basis_type+2)/2, N_T-1);

un_cont = zeros(size(DGT, 2), (basis_type+1)*(basis_type+2)/2, N_T-1); 
for n=0:N_T-2
    tic
    t_0=n*delta_t; 
    t_1=(n+1)*delta_t; 
    t_2=(n+2)*delta_t;
    t=t_1;
    fprintf('当前时间为%f/%f.',t_2,T);
    
    % Assemble the load vector
    b1=generate_load_vector_local_DG('function_f',t,DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
    %{
    b2e=generate_load_vector_local_DG_boundary_edge('integrand_3_on_edges_e',t,DGM,DGT,DGE,DG_edge_flag...
                                                      ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                      ,basis_type,boundary_edges,eipsilon,sigma_2,beta);                                    
    b2z=generate_load_vector_local_DG_boundary_edge('integrand_3_on_edges_z',t,DGM,DGT,DGE,DG_edge_flag...
                                                      ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                      ,basis_type,boundary_edges,eipsilon,sigma_2,beta);
    %}
    %b=b1+b2e+b2z;
    b = b1
    clear b1 
    %b=b*(delta_t)^2+2*M*r_1-((delta_t)^2/2*Az+M-M1*delta_t/2-Ae*delta_t/2)*r_0;
    %left_matrix=M+M1*delta_t/2+Ae*delta_t/2+(delta_t)^2/2*Az;
    b=(delta_t)^2*b+2*M*r_1-M*r_0+(delta_t)/2*Ae*r_0-(delta_t)^2/2*Az*r_0;
    left_matrix=M+delta_t/2*Ae+(delta_t)^2/2*Az;
    % Solve
    r_2=left_matrix\b;
    clear b
    fprintf('加细第%d次的网格,第%d/%d个时间步',Iter_refinement,n+2,N_T);
    toc
    %% 每一步的Error
    L2_error=L2_H1_error(DGM,DGT,r_2,'accurate_solution',t_2,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    temp1=L2_H1_error(DGM,DGT,r_2,'accurate_solution_x_derivative',t_2,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    temp2=L2_H1_error(DGM,DGT,r_2,'accurate_solution_y_derivative',t_2,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    H1_error=sqrt(temp1^2+temp2^2);  
    fprintf('当前时间的L^2误差为%f，H^1误差为%f.\n',L2_error,H1_error);
    % Evaluate
    r_0=r_1;
    r_1=r_2;
    %提取picked三角元的数据并存储
    un_1(:,n+1) = r_2((n_1-1)*(basis_type+1)*(basis_type+2)/2+1:n_1*(basis_type+1)*(basis_type+2)/2,1);
    un_2(:,n+1) = r_2((n_2-1)*(basis_type+1)*(basis_type+2)/2+1:n_2*(basis_type+1)*(basis_type+2)/2,1);
    for j=1:size(DGT, 2)
        un_cont(j,:,n+1) = r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
    end
    %% Figure
    uh=zeros(3*size(DGT,2),1);
    for m=1:size(DGT,2)
           vertices=DGM(:,DGT(:,m));
           uh_local=r_2((m-1)*(basis_type+1)*(basis_type+2)/2+1:m*(basis_type+1)*(basis_type+2)/2,1);
           uh(3*m-2,1)=fe_solution(DGM(1,DGT(1,m)),DGM(2,DGT(1,m)),uh_local,0,0,basis_type);
           uh(3*m-1,1)=fe_solution(DGM(1,DGT(2,m)),DGM(2,DGT(2,m)),uh_local,0,0,basis_type);
           uh(3*m,1)=fe_solution(DGM(1,DGT(3,m)),DGM(2,DGT(3,m)),uh_local,0,0,basis_type);
    end

    figure(2024);
    trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);%%%函数trisurf的用法
    colormap("jet")
    shading interp
    xlabel('x'), ylabel('y')

    figure(2028);
    trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);%%%函数trisurf的用法
    colormap("jet")
    shading interp
    view(2)
    cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
    %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
    % 添加 X 和 Y 轴的名称（对于二维视角）
    xlabel('x');  % 自定义 X 轴名称
    ylabel('y');  % 自定义 Y 轴名称
    
    if n==20 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end
    if n==148 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end
    if n==376
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==415 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end
   

    if n==422 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end
    if n==426 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==435 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==445 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==470 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==550 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==678 
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==808
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==821
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==840
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==870
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end

    if n==920
       uh=zeros(3*size(DGT,2),1);
            for j=1:size(DGT,2)
                vertices=DGM(:,DGT(:,j));
                uh_local=r_2((j-1)*(basis_type+1)*(basis_type+2)/2+1:j*(basis_type+1)*(basis_type+2)/2,1);
                uh(3*j-2,1)=fe_solution(DGM(1,DGT(1,j)),DGM(2,DGT(1,j)),uh_local,0,0,basis_type);
                uh(3*j-1,1)=fe_solution(DGM(1,DGT(2,j)),DGM(2,DGT(2,j)),uh_local,0,0,basis_type);
                uh(3*j,1)=fe_solution(DGM(1,DGT(3,j)),DGM(2,DGT(3,j)),uh_local,0,0,basis_type);
            end
            %%%%%%画三维图像
            figure(20)
            DGT_plot = DGT';
            %利用trisurf函数绘制三维曲面图
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %函数trisurf的用法
                                                   %trisurf用于创建三位表面图(需要4个基本参数)
            shading interp  %使用插值方式对图形进行着色
            colormap(jet);  %hot,summer
            view(30, 30);% 设置视角
            %lighting gouraud;  % 光滑光照
            material shiny;    % 增强光泽
            %camlight;          % 在当前视角添加光源
            xlabel('x'), ylabel('y')  %设置坐标轴标签
            %ztickformat('%.2e');   %z轴刻度以科学计数法显示
            t_pre_index=n*delta_t
            %title(sprintf('3D images of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_3d = fullfile(folder, fileName_3d);  % 合并路径和文件名
            saveas(gcf, fullFileName_3d);  % 保存文件
            %%%%%%画二维等高图
            %%%%%%即把三维图切换到二维视角
            view(2);
            cb = colorbar;% 添加色条，并获取色条句柄，并将其句柄保存到变量 cb
            %set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% 设置色条刻度为科学计数法
            % 添加 X 和 Y 轴的名称（对于二维视角）
            xlabel('x');  % 自定义 X 轴名称
            ylabel('y');  % 自定义 Y 轴名称
            %title(sprintf('contours of the viscoelastic wave propagation \n in wellmodel at t=%d * %.6f',t_pre_index,delta_t));
            %保存图像
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % 将数值插入文件名中
            folder = 'E:\matlab codes\my code\seismicplot\layeredmedium';
            fullFileName_cont = fullfile(folder, fileName_cont);
            saveas(gcf, fullFileName_cont);  % 保存文件
    end
end
%待所有迭代完成后利用提取到的数据作图
%注意保存下来不是最终的数值解  还需要带进线性组合式算一下
uhn_1 = zeros(3,N_T-1);    %每一行是一个顶点每次迭代的数值解
uhn_2 = zeros(3,N_T-1);
for j=1:N_T-1
    uhn_1(1,j)=fe_solution(DGM(1,DGT(1,n_1)),DGM(2,DGT(1,n_1)),un_1(:,j),0,0,basis_type);
    uhn_1(2,j)=fe_solution(DGM(1,DGT(2,n_1)),DGM(2,DGT(2,n_1)),un_1(:,j),0,0,basis_type);
    uhn_1(3,j)=fe_solution(DGM(1,DGT(3,n_1)),DGM(2,DGT(3,n_1)),un_1(:,j),0,0,basis_type);
    uhn_2(1,j)=fe_solution(DGM(1,DGT(1,n_2)),DGM(2,DGT(1,n_2)),un_2(:,j),0,0,basis_type);
    uhn_2(2,j)=fe_solution(DGM(1,DGT(2,n_2)),DGM(2,DGT(2,n_2)),un_2(:,j),0,0,basis_type);
    uhn_2(3,j)=fe_solution(DGM(1,DGT(3,n_2)),DGM(2,DGT(3,n_2)),un_2(:,j),0,0,basis_type);
end
%各个三角形的三个点对应的线采取相同样式
%创建用于画图的时间数据（x轴数据）
t_plot = 2*delta_t : delta_t : T;
y11_plot = uhn_1(1,:);   
y12_plot = uhn_1(2,:); 
y13_plot = uhn_1(3,:); 
y21_plot = uhn_2(1,:); 
y22_plot = uhn_2(2,:); 
y23_plot = uhn_2(3,:); 

% 使用 plot 函数绘制折线图
figure(2025);         % 创建新图窗口
plot(t_plot, y11_plot, '--r', 'DisplayName', 'receiver1');  % 绘制折线图，使用 '--r' 让数据点以红色虚线显示
hold on;
plot(t_plot, y12_plot, '--r', 'DisplayName', 'receiver2');
hold on;
plot(t_plot, y13_plot, '--r', 'DisplayName', 'receiver3');
hold on;
plot(t_plot, y21_plot, '-b', 'DisplayName', 'receiver4');
hold on;
plot(t_plot, y22_plot, '-b', 'DisplayName', 'receiver5');
hold on;
plot(t_plot, y23_plot, '-b', 'DisplayName', 'receiver6');

title('simulation record of viscoelastic wave in submodel', 'FontSize', 14);   % 添加标题
xlabel('time of propagation');          % 添加 x 轴标签
ylabel('amplitude');     % 添加 y 轴标签
legend('show', 'FontSize', 10);
grid on;              % 打开网格

%%%%%%%%画地震炮记录图（合成记录）
%%%%%用meshgrid  把各顶点x坐标或y坐标与时间t来meshgrid
t_plot_cont = 2*delta_t : delta_t : T; %竖轴
y_plot_cont = DGM(2,:); %横轴
x_plot_cont = DGM(1,:);
[Y_plot_cont, T_plot_cont] = meshgrid(y_plot_cont,t_plot_cont); %创建网格
[X_plot_cont, T_plot_cont] = meshgrid(x_plot_cont,t_plot_cont); %创建网格
%for j=1:N_T-1
%    for i=1:size(DGT,2)
%        disp(size(un_cont(i,:,j)))
%    end
%end
%%%%%%最关键的问题是提取温值
%%%%%%即各点处每个迭代时刻的数值解值
uhn_cont = zeros(size(DGT,2), 3 , N_T-1);
for j=1:N_T-1
    for i=1:size(DGT,2)
        uhn_cont(i,1,j)=fe_solution(DGM(1,DGT(1,i)),DGM(2,DGT(1,i)),un_cont(i,:,j)',0,0,basis_type);
        uhn_cont(i,2,j)=fe_solution(DGM(1,DGT(2,i)),DGM(2,DGT(2,i)),un_cont(i,:,j)',0,0,basis_type);
        uhn_cont(i,3,j)=fe_solution(DGM(1,DGT(3,i)),DGM(2,DGT(3,i)),un_cont(i,:,j)',0,0,basis_type);
    end
end
%%%%%uhn_cont  第几个三角形（行）   第几个顶点处的数值解值（列）  第几个迭代时刻（层）
% 展平三维张量 u 为一个二维矩阵，维度为 (顶点数, 时间步数)
num_triangles = size(uhn_cont, 1); % 三角形数量
num_vertices_per_triangle = size(uhn_cont, 2); % 每个三角形的顶点数
num_time_steps = size(uhn_cont, 3); % 时间步数

% 将顶点的数值解展开为二维矩阵：每行是顶点，每列是时间步数
% 例如，每个三角形的三个顶点 (i, j, :) 在不同时刻的数值解都会按顺序存储
uhn_cont_reshaped = reshape(uhn_cont, num_triangles * num_vertices_per_triangle, num_time_steps);
%Y_plot_cont, T_plot_cont = meshgrid(y_plot_cont,t_plot_cont); %创建网格
% 绘制等高线图
figure(2026);
%disp(size(Y_plot_cont))
%disp(size(T_plot_cont))
%disp(size(uhn_cont_reshaped))
contourf(Y_plot_cont', T_plot_cont', uhn_cont_reshaped, 'LineStyle', 'none'); % 20表示等高线的数量
colorbar; % 添加颜色条
xlabel('y');
ylabel('time');
title('contour plot of vicsoelastic wave y_t');

figure(2027);
%disp(size(X_plot_cont))
contourf(X_plot_cont', T_plot_cont', uhn_cont_reshaped, 'LineStyle', 'none'); % 20表示等高线的数量
colorbar; % 添加颜色条
xlabel('x');
ylabel('time');
title('contour plot of vicsoelastic wave x_t');

                                             
%% 最终时间的Error

L2_error=L2_H1_error(DGM,DGT,r_2,'accurate_solution',T,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);

temp1=L2_H1_error(DGM,DGT,r_2,'accurate_solution_x_derivative',T,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp2=L2_H1_error(DGM,DGT,r_2,'accurate_solution_y_derivative',T,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
H1_error=sqrt(temp1^2+temp2^2);
                                    

