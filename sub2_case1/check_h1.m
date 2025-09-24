function check
% An example: Backward Euler, NIPG, P1(E) basis function, ordinary penalty type.
% One can modify following parameters for Backward Euler(theta=1), Forward Euler(theta=0), Crank-Nicolson discretizations(theta=1/2), 
%                                         different schemes such as NIPG(eipsilon=1), IIPG(eipsilon=0), SIPG(eipsilon=-1), OBB(eipsilon=1,sigma=0,k>=2),
%                                         basis functions of different basis types(basis_type=1,2,3...), 
%                                         penalty(beta=1) and super-penalty(beta=3) methods
% with appropriate proportional relation between time step delta_t and mesh size h (i.e., 1/2^Iter_refinement)!
clc
clear
cd('E:/matlab codes/my code/sub2_case1_ex1_explicit_min_CFL')
basis_type=1;
eipsilon=-1;
sigma=200*(basis_type+1)^2;
beta=1;
theta=1;
T=1;
% Relation between delta_t and h
Iter_refinement=1;%N_T=2^(Iter_refinement*2);
[L2_error_1,H1_error_1,DG_broken_error_1,data_1,CPU_time_1]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,T);

Iter_refinement=2;%N_T=2^(Iter_refinement*2);
[L2_error_2,H1_error_2,DG_broken_error_2,data_2,CPU_time_2]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,T);

Iter_refinement=3;%N_T=2^(Iter_refinement*2);
[L2_error_3,H1_error_3,DG_broken_error_3,data_3,CPU_time_3]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,T);

Iter_refinement=4;%N_T=2^(Iter_refinement*2);
[L2_error_4,H1_error_4,DG_broken_error_4,data_4,CPU_time_4]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,T);

%Iter_refinement=5;N_T=2^(Iter_refinement*2);
%[L2_error_5,H1_error_5,DG_broken_error_5,data_5,CPU_time_5]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%  Iter_refinement=6;N_T=2^(Iter_refinement*2);
%  [L2_error_6,H1_error_6]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������primal��ʽ������
%%
h=[1/2;1/4;1/8;1/16];%1/32]%;1/64];
eL2=[L2_error_1;L2_error_2;L2_error_3;L2_error_4];%L2_error_5]%;L2_error_6];%%%%%%%%%%%%%%nά����
eH1=[H1_error_1;H1_error_2;H1_error_3;H1_error_4];%H1_error_5]%;H1_error_6];
%eInf=[uh_infinity_error_8;uh_infinity_error_16;uh_infinity_error_32;uh_infinity_error_64];%%%%%��ڵ�������(�����)
eDGbroken=[DG_broken_error_1;DG_broken_error_2;DG_broken_error_3;DG_broken_error_4];%DG_broken_error_5];%DG_broken_error_6];
CPU_time=[CPU_time_1;CPU_time_2;CPU_time_3;CPU_time_4];%CPU_time_5];

L=length(h);
eL2Order=zeros(L-1,1);%%%%%%%%%%%%n-1ά��������Ϊn������ֻ�ܼ���n-1��������
eH1Order=zeros(L-1,1);
%eInfOrder=zeros(L-1,1);
eDGbroken_Order=zeros(L-1,1);

for i=1:L-1
    eL2Order(i)=log(eL2(i)/eL2(i+1))/log(2);
    eH1Order(i)=log(eH1(i)/eH1(i+1))/log(2);
    eDGbroken_Order(i)=log(eDGbroken(i)/eDGbroken(i+1))/log(2);
%  eInfOrder(i)=log(eInf(i)/eInf(i+1))/log(2);
end

T1=table(h,eL2,eH1,eDGbroken,CPU_time);
disp('the errors for example 1_implicit_min:');
format short e;
disp(T1);
h=h(2:L,1);
disp('the convergence rate for example 1_implicit_min:');
format short;
T2=table(h,eL2Order,eH1Order,eDGbroken_Order);
disp(T2);
disp('p=1,N_T=2^(Iter_refinement*2)');

figure;
h=[1/2;1/4;1/8;1/16;1/32]%;1/64];
L=length(h); %ע��h������ı��˱���  h���ȱ�ʾ�������������ݵ�����
%X=linspace(1,L,L); %������һ����1�� L �ĵȼ������� X��������ʾ��������ݵ㣬���������������ͼ���ϻ�������
X = 1:L;
slopeL2=[1/2*eL2(1);1/8*eL2(1);1/32*eL2(1);1/128*eL2(1);1/512*eL2(1)];%1/2048*eL2(1)];
slopeH1=[1/2*eH1(1);1/4*eH1(1);1/8*eH1(1);1/16*eH1(1);1/32*eH1(1)];
slopeDGbroken=[1/2*eDGbroken(1);1/4*eDGbroken(1);1/8*eDGbroken(1);1/16*eDGbroken(1);1/32*eDGbroken(1)];%1/64*eDGbroken(1)];
%�������ڻ���б�ʲο��ߵ����ݣ�slopeL2 �� slopeDGbroken �ֱ��Ӧ�� eL2 �� eDGbroken ��б�����У�����Щб�ʲο������Ƚ�ʵ����������б��
semilogy(X,eL2,'b-s',X,eH1,'c-x',X,eDGbroken,'g-O',X,slopeL2,'k--',X,slopeH1,'y--',X,slopeDGbroken,'r--','LineWidth',1.9); %��ͼ�δ����л����������ߣ��ַ��������ֱ�ָ��������ʽ����ɫ����ɫʵ�ߴ��� eL2����ɫԲȦ���� eDGbroken����ɫ���ߴ��� slopeL2����ɫ���ߴ��� slopeDGbroken)
% ,N,N.^(-5/2),'k--');
h2=legend('$||u-u_h||_0$','$||u-u_h||_1$','$||u-u_h||_{\star}$','$ch^2$','$ch$','$ch$','Location','southwest'); %��ͼ�ε����½Ǵ���ͼ�� h2��������ʶͼ�е�����
set(h2,'Interpreter','latex')  %����ͼ�����ַ����Ľ�����Ϊ LaTeX���Ա�֧�� LaTeX ��ʽ���ı�
set(h2,'fontsize',15);  %����ͼ�����ַ����������СΪ15
set(gca,'XTick',1:L,'XLim',[0,L+1],'Ylim',[1e-5 1e2]);  %����ͼ�εĺ���̶ȣ�������1�� L ������
xlabel('refined iteration');
ylabel('log(error)');
grid on;
set(gca,'GridLineStyle','--'); %������������ʽ
set(gca,'MarkerSize',6);  %���ӱ�Ǵ�С
%title('');
%grid on;

