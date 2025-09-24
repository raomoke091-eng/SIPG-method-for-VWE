function check_h_2
% An example: Backward Euler, NIPG, P1(E) basis function, ordinary penalty type.
% One can modify following parameters for Backward Euler(theta=1), Forward Euler(theta=0), Crank-Nicolson discretizations(theta=1/2), 
%                                         different schemes such as NIPG(eipsilon=1), IIPG(eipsilon=0), SIPG(eipsilon=-1), OBB(eipsilon=1,sigma=0,k>=2),
%                                         basis functions of different basis types(basis_type=1,2,3...), 
%                                         penalty(beta=1) and super-penalty(beta=3) methods
% with appropriate proportional relation between time step delta_t and mesh size h (i.e., 1/2^Iter_refinement)!

basis_type=2;
eipsilon=-1;
sigma=1000;
beta=1;
theta=1;
T=0.1;
% Relation between delta_t and h
Iter_refinement=1;N_T=2^(Iter_refinement*2);%delta_t=h^2
[L2_error_1,H1_error_1,DG_broken_error_1,data_1,CPU_time_1]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
save('data_1.mat','data_1');

Iter_refinement=2;N_T=2^(Iter_refinement*2);
[L2_error_2,H1_error_2,DG_broken_error_2,data_2,CPU_time_2]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
save('data_2.mat','data_2');

Iter_refinement=3;N_T=2^(Iter_refinement*2);
[L2_error_3,H1_error_3,DG_broken_error_3,data_3,CPU_time_3]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
save('data_3.mat','data_3');

Iter_refinement=4;N_T=2^(Iter_refinement*2);
[L2_error_4,H1_error_4,DG_broken_error_4,data_4,CPU_time_4]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
save('data_4.mat','data_4');

Iter_refinement=5;N_T=2^(Iter_refinement*2);
[L2_error_5,H1_error_5,DG_broken_error_5,data_5,CPU_time_5]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
save('data_5.mat','data_5');

%Iter_refinement=6;N_T=2^(Iter_refinement*2);
%[L2_error_6,H1_error_6,DG_broken_error_6,data_6]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
%save('data_6.mat','data_6');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������primal��ʽ������
%%
h=[1/2;1/4;1/8;1/16;1/32];%1/64];
eL2=[L2_error_1;L2_error_2;L2_error_3;L2_error_4;L2_error_5];%L2_error_6];%������
eH1=[H1_error_1;H1_error_2;H1_error_3;H1_error_4;H1_error_5];%H1_error_6];
eDGbroken=[DG_broken_error_1;DG_broken_error_2;DG_broken_error_3;DG_broken_error_4;DG_broken_error_5];%DG_broken_error_6];
CPU_time=[CPU_time_1;CPU_time_2;CPU_time_3;CPU_time_4;CPU_time_5];

%eInf=[uh_infinity_error_8;uh_infinity_error_16;uh_infinity_error_32;uh_infinity_error_64];%%%%%��ڵ�������(�����)

L=length(h);
eL2_Order=zeros(L-1,1);%%%%%%%%%%%%n-1ά��������Ϊn������ֻ�ܼ���n-1��������  ������һ�����������
eH1_Order=zeros(L-1,1);
eDGbroken_Order=zeros(L-1,1);
%eInfOrder=zeros(L-1,1);

for i=1:L-1
    eL2_Order(i)=log(eL2(i)/eL2(i+1))/log(2);
    eH1_Order(i)=log(eH1(i)/eH1(i+1))/log(2);
    eDGbroken_Order(i)=log(eDGbroken(i)/eDGbroken(i+1))/log(2);
%  eInfOrder(i)=log(eInf(i)/eInf(i+1))/log(2);
end

T1=table(h,eL2,eH1,eDGbroken,CPU_time);   %����һ�����T1������������ h��eL2 �� eDGbroken ���ݣ��������������ֱ��Ӧ��������
disp('the errors for example 2:');   %��ʾ�ַ���
format short e;   %���������ʽ��������������ʾ��ʽ����Ϊ�̵Ŀ�ѧ������
disp(T1);  %��ʾ��� T1 �е�����

h1=h(2:L,1);   %����һ���µ������� h1���� h �ĵڶ���Ԫ�ؿ�ʼ�����һ��Ԫ�ص��Ӽ���ȥ���˵�һ��Ԫ�أ�
T2=table(h1,eL2_Order,eH1_Order,eDGbroken_Order);  %%����һ�����T2������������ h1��eL2_Order �� eDGbroken_Order ���ݣ��������������ֱ��Ӧ��������
disp('the spatial convergence rate for example 2:');
format short;  %�������ʽ�ָ�ΪĬ�����ã�����ʾ�����͸�����
disp(T2);
disp('p=2,N_T=2^(Iter_refinement*2),at T=0.01,sigma=1000');


%% Figure
figure(2024);

L=length(h); %ע��h������ı��˱���  h���ȱ�ʾ�������������ݵ�����
X=linspace(1,L,L); %������һ����1�� L �ĵȼ������� X��������ʾ��������ݵ㣬���������������ͼ���ϻ�������

slopeL2=[1/2*eL2(1);(1/2)^4*eL2(1);(1/2)^7*eL2(1);(1/2)^10*eL2(1);(1/2)^13*eL2(1)];%(1/2)^16*eL2(1)];
slopeDGbroken=[1/2*eDGbroken(1);1/8*eDGbroken(1);1/32*eDGbroken(1);1/128*eDGbroken(1);1/512*eDGbroken(1)];%1/2048*eDGbroken(1)];
%�������ڻ���б�ʲο��ߵ����ݣ�slopeL2 �� slopeDGbroken �ֱ��Ӧ�� eL2 �� eDGbroken ��б�����У�����Щб�ʲο������Ƚ�ʵ����������б��
semilogy(X,eL2,'b-s',X,eDGbroken,'g-O',X,slopeL2,'k--',X,slopeDGbroken,'r--'); %��ͼ�δ����л����������ߣ��ַ��������ֱ�ָ��������ʽ����ɫ����ɫʵ�ߴ��� eL2����ɫԲȦ���� eDGbroken����ɫ���ߴ��� slopeL2����ɫ���ߴ��� slopeDGbroken)
% ,N,N.^(-5/2),'k--');
h2=legend('$||u-u_h||_0$','$||u-u_h||_{\star}$','$ch^3$','$ch^2$','Location','southwest'); %��ͼ�ε����½Ǵ���ͼ�� h2��������ʶͼ�е�����
set(h2,'Interpreter','latex')  %����ͼ�����ַ����Ľ�����Ϊ LaTeX���Ա�֧�� LaTeX ��ʽ���ı�
set(h2,'fontsize',15);  %����ͼ�����ַ����������СΪ15
set(gca,'XTick',1:L);  %����ͼ�εĺ���̶ȣ�������1�� L ������
xlabel('refinement iteration');
ylabel('Error');
%grid on;


