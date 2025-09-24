function [L2_error,H1_error, DG_broken_error,r_2,CPU_time]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T,constant,coefficient)
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
%DGM=[0,0;pi/10,0;0,pi/10;pi/10,pi/10;0,pi/10;pi/10,0]';
DGM=[0,0;1,0;0,1;1,1;0,1;1,0]';
%DGM=[0,0;1/5,0;0,1/5;1/5,1/5;0,1/5;1/5,0]';
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
delta_t=T/N_T;

b_0=generate_load_vector_local_DG('accurate_solution',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_0=M\b_0;%%u_0为u(0)的L2投影
%u_0=Ae\b_0;%%u_0为u(0)的pi_a投影


b_1=generate_load_vector_local_DG('accurate_solution_time_derivate',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
v_0=M\b_1;%%v_0为u'(0)的L2投影   
%v_0=Ae\b_1;%%v_0为u'(0)的pi_a投影  


%B_0=generate_load_vector_local_DG('function_f',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
%u_tilde=M\(B_0-Ae*v_0-Az*u_0);
%u_tilde=M\(B_0-M1*v_0-Ae*v_0-Az*u_0);

b_2=generate_load_vector_local_DG('accurate_solution_time_derivate_2',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_tilde=M\b_2;%%也可以选为u''(0)的L2投影
%u_tilde=Ae\b_2;%%也可以选为u''(0)的pi_a投影


r_1=u_0+delta_t*v_0+(delta_t)^2/2*u_tilde;
r_0=u_0;

%% LOOP
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
    b = b1;
    clear b1 
    %b=b*(delta_t)^2+2*M*r_1-((delta_t)^2/2*Az+M-M1*delta_t/2-Ae*delta_t/2)*r_0;
    %left_matrix=M+M1*delta_t/2+Ae*delta_t/2+(delta_t)^2/2*Az;
    b=(delta_t)^2*b+2*M*r_1-M*r_0+(delta_t)/2*Ae*r_0-(delta_t)^2/2*Az*r_0;
    left_matrix=M+delta_t/2*Ae+(delta_t)^2/2*Az;
    % Solve
    r_2=left_matrix\b;
    clear b
    CPUTime(n+1,1) = toc;
    fprintf('加细第%d次的网格,第%d/%d个时间步',Iter_refinement,n+2,N_T);
    disp(['历时CPUTime: ', num2str(CPUTime(n+1,1)), ' 秒']);
    %% 每一步的Error
    %L2_error=L2_H1_error(DGM,DGT,r_2,'accurate_solution',t_2,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    %temp1=L2_H1_error(DGM,DGT,r_2,'accurate_solution_x_derivative',t_2,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    %temp2=L2_H1_error(DGM,DGT,r_2,'accurate_solution_y_derivative',t_2,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    %H1_error=sqrt(temp1^2+temp2^2); 
    
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
    fprintf('当前时间的L^2误差为%f，H^1误差为%f.\n，DG离散误差为%f.\n',L2_error,H1_error,DG_broken_error);
    % Evaluate
    r_0=r_1;
    r_1=r_2;
    

    %% Figure
    uh=zeros(3*size(DGT,2),1);
    for n=1:size(DGT,2)
           vertices=DGM(:,DGT(:,n));
           uh_local=r_2((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
           uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
           uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
           uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
    end

    figure(2024);
    trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);%%%函数trisurf的用法
    %colormap("jet")
    shading interp
    xlabel('x'), ylabel('y')
end                                             
%% 最终时间的Error
%{
L2_error=L2_H1_error(DGM,DGT,r_2,'accurate_solution',T,basis_type,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);

temp1=L2_H1_error(DGM,DGT,r_2,'accurate_solution_x_derivative',T,basis_type,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp2=L2_H1_error(DGM,DGT,r_2,'accurate_solution_y_derivative',T,basis_type,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
H1_error=sqrt(temp1^2+temp2^2);
%}

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
for i=1:size(CPUTime)
    CPU_time=CPU_time+CPUTime(i,1);
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

