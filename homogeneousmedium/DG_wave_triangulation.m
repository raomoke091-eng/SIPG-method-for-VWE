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
%% �������ǻ��ִ�����ʼ���񣨰������������Σ�����ά���Σ�
DGM=[0,0;1,0;0,1;1,1;0,1;1,0]'; %DGM�����������ڵ����꣨ÿһ�б�ʾһ���ڵ㣬ÿһ�д���һ��ά�ȣ�
DGT=[1,2,3;4,5,6]';             %DGT������������εĶ������/����  ��ӦDGM�Ľڵ�����
                                %���ǰѸ��������е�ά�ȣ�DGM�����ǽڵ������DGT�����������θ�����


%% Begin uniform refinement
for i=1:Iter_refinement      %Iter_refinement��������
    [DGM,DGT]=uniformrefine_triangle(DGM,DGT);    %�����Զ��庯��uniformrefine_triangle(DGM,DGT)�Գ�ʼDGM��DGT����ϸ�����ָ�ɸ�С�������Σ���������ֱ��ʣ�
                                                  %ÿforѭ��һ�Σ���ÿ����һ�Σ�����ϸ��һ�Σ�����Ϊ�ָ��4���������Σ�
end

%�������Ĳ�������ɵ���ϸ��֮��
[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT); 
[un_used,interior_edges]=find(DGE(6,:));
%ͨ������ DGE ����ĵ� 6 �У��ҵ��ڲ��ߵ������洢�ھ��� interior_edges ��
%find ��������DGE(6,:)�з���Ԫ�ص������� un_used �������� interior_edges
boundary_edges=setdiff(1:size(DGE,2),interior_edges);
%1:size(DGE,2) ������һ���� 1 �� DGE ������������������ʾ���бߵ�����
%ͨ�� setdiff ���������ų��ڲ��ߵ�������ʣ�¼�Ϊ�߽����������Щ���������� boundary_edges ������
%[rows, cols] = size(DGE);
%for i = 1:rows
%    for j = 1:cols
%        fprintf('%d\t', DGE(i, j));
%    end
%    fprintf('\n');
%end


[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(12);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(12);


%% Assemble stiffness matrix��װ�նȾ���
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


%% LOOP����initial_steps:

% Because the load vector depends on the time, we must write it into the time cycle
delta_t=T/N_T;  

b_0=generate_load_vector_local_DG('accurate_solution',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_0=M\b_0;%%u_0Ϊu(0)��L2ͶӰ

b_1=generate_load_vector_local_DG('accurate_solution_time_derivate',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
v_0=M\b_1;%%v_0Ϊu'(0)��L2ͶӰ                                

B_0=generate_load_vector_local_DG('function_f',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
u_tilde=M\(B_0-Aa*v_0-Ab*u_0);
% b_2=generate_load_vector_local_DG('accurate_solution_time_derivate_2',0,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);
% u_tilde=M\b_2;%%Ҳ����ѡΪu''(0)��L2ͶӰ

r_1=u_0+delta_t*v_0+(delta_t)^2/2*u_tilde; %t=1ʱ��U1
r_0=u_0; %t=0ʱ��U0

%% LOOP
CPUTime=zeros(1,N_T-1);
for n=0:N_T-2
    tic %��ʱ
    t_0=n*delta_t; 
    t_1=(n+1)*delta_t; 
    t_2=(n+2)*delta_t;
    t=t_1;
    fprintf('��ǰʱ��Ϊ%f/%f.',t_2,T);
    
    b=generate_load_vector_local_DG('function_f',t,DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type);% Assemble the load vector
    b=(delta_t)^2*b+2*M*r_1-M*r_0+(delta_t)/2*Aa*r_0-(delta_t)^2/2*Ab*r_0;  
    left_matrix=M+delta_t/2*Aa+(delta_t)^2/2*Ab;
    
    % Solve
    r_2=left_matrix\b;
    clear b
    fprintf('��ϸ��%d�ε�����,��%d/%d��ʱ�䲽',Iter_refinement,n+2,N_T);
    CPUTime(1,n+1) = toc;
    disp(['��ʱCPUTime: ', num2str(CPUTime(1,n+1)), ' ��']);
    %% ÿһ����Error
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
    
    fprintf('��ǰʱ���L^2���Ϊ%f��H^1���Ϊ%f��DG��ɢ���Ϊ%f.\n',L2_error,H1_error,DG_broken_error);
    % Evaluate
    r_0=r_1;
    r_1=r_2;

    
    %% Figure��ͼ
    %���������������ϵ���ֵ�� uh ����άͼ��
    uh=zeros(3*size(DGT,2),1);
    %����һ����СΪ 3*size(DGT,2) �������� uh�����ڴ洢��ֵ����ÿ�����ǵ�Ԫ���������㴦��ֵ��ÿ�������η���3��ά�ȣ�

    for n=1:size(DGT,2)
        vertices=DGM(:,DGT(:,n));%2*3
        uh_local=r_2((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);%��ǰ���ǵ�Ԫ�ϵ���ֵ�⣨������ϵ����
        uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
        uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
        uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
        %�ֱ�������ǵ�Ԫ���������㴦����ֵ��
    end

    figure(2024);
    trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);  %����trisurf���÷�
    shading interp  %ʹ�ò�ֵ��ʽ��ͼ�ν�����ɫ
    xlabel('x'), ylabel('y')  %�����������ǩ
    
    %��ʱ����ͼ  10��ʱ��+��ʼʱ��
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
            %%%%%%����άͼ��
            figure(i)
            %�Ȱѻ�ͼ��������������ݶ�Ӧ����  �����ζ���������u_h�����ֵҪ��Ӧ��ȷ
            %Ŀǰu_h����һ������  ������1����1��������1����2��������1����3��������2����1��������2���㣬...���Դ����ƣ�
            %��DGM����Ķ��������u_h��Ӧ���� Ҫ���DGM��DGT�����������Ϣ
            %ע�⵽���ֶ�Ӧ��ϵ�������ģ���DGT�е�����Ǽ����������DGM�е�������ڵڼ���
            %���ֻ��Ҫ��DGMת��һ�¾��ܵõ�һ����һ��Ϊx�������ݣ��ڶ���Ϊy�������ݵľ���
            DGT_plot = DGT';
            %����trisurf����������ά����ͼ
            trisurf(DGT_plot,DGM(1,:)',DGM(2,:)',uh);  %����trisurf���÷�
                                                   %trisurf���ڴ�����λ����ͼ(��Ҫ4����������)
            shading interp  %ʹ�ò�ֵ��ʽ��ͼ�ν�����ɫ
            %colormap(gray);
            view(30, 50);% �����ӽ�
            lighting gouraud;  % �⻬����
            material shiny;    % ��ǿ����
            camlight;          % �ڵ�ǰ�ӽ���ӹ�Դ
            xlabel('x'), ylabel('y')  %�����������ǩ
            ztickformat('%.2e');   %z��̶��Կ�ѧ��������ʾ
            title(sprintf('3D images of the viscoelastic wave propagation \n in homofeneous medium at t=%d * %.6f',t_pre_index,delta_t));
            fileName_3d = sprintf('3D_surface_%d.png', t_pre_index);  % ����ֵ�����ļ�����
            saveas(gcf, fileName_3d);  % �����ļ�
            

            %%%%%%����ά�ȸ�ͼ
            %%%%%%������άͼ�л�����ά�ӽ�
            view(2);
            cb = colorbar;% ���ɫ��������ȡɫ������������������浽���� cb
            set(cb, 'YTickLabel', sprintfc('%.2e', cb.Ticks));% ����ɫ���̶�Ϊ��ѧ������
            % ��� X �� Y ������ƣ����ڶ�ά�ӽǣ�
            xlabel('x');  % �Զ��� X ������
            ylabel('y');  % �Զ��� Y ������
            title(sprintf('contours of the viscoelastic wave propagation \n in homofeneous medium at t=%d * %.6f',t_pre_index,delta_t));
            %����ͼ��
            fileName_cont  = sprintf('2d_cont_%d.png', t_pre_index);  % ����ֵ�����ļ�����
            saveas(gcf, fileName_cont);  % �����ļ�
        end
    end
end
                                             
%% ����ʱ��Tʱ�̵�Error
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
% ���ÿ�����񻮷����õļ���ʱ��cpu time
CPU_time=0;
for i=1:size(CPUTime,2)
    CPU_time=CPU_time+CPUTime(1,i);
end

% ����ʾ����������
node = DGM'; % �ڵ�����
element = DGT'; % ��Ԫ�ڵ�����

% ��������ͼ��
figure;
trimesh(element, node(:,1), node(:,2), 'LineWidth', 1.5); % ʹ�� trimesh ������������������
xlabel('X'); ylabel('Y'); % ����������ǩ
title('DG Mesh'); % ��ӱ���
axis equal; % ����������������
