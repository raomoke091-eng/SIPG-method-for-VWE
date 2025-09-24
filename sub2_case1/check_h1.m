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
%这里给出primal格式的误差表：
%%
h=[1/2;1/4;1/8;1/16];%1/32]%;1/64];
eL2=[L2_error_1;L2_error_2;L2_error_3;L2_error_4];%L2_error_5]%;L2_error_6];%%%%%%%%%%%%%%n维向量
eH1=[H1_error_1;H1_error_2;H1_error_3;H1_error_4];%H1_error_5]%;H1_error_6];
%eInf=[uh_infinity_error_8;uh_infinity_error_16;uh_infinity_error_32;uh_infinity_error_64];%%%%%求节点最大误差(无穷范数)
eDGbroken=[DG_broken_error_1;DG_broken_error_2;DG_broken_error_3;DG_broken_error_4];%DG_broken_error_5];%DG_broken_error_6];
CPU_time=[CPU_time_1;CPU_time_2;CPU_time_3;CPU_time_4];%CPU_time_5];

L=length(h);
eL2Order=zeros(L-1,1);%%%%%%%%%%%%n-1维向量，因为n步计算只能计算n-1步收敛率
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
L=length(h); %注意h在上面改变了变量  h长度表示迭代次数或数据点数量
%X=linspace(1,L,L); %创建了一个从1到 L 的等间距的向量 X，用来表示横轴的数据点，这个向量将用于在图形上绘制数据
X = 1:L;
slopeL2=[1/2*eL2(1);1/8*eL2(1);1/32*eL2(1);1/128*eL2(1);1/512*eL2(1)];%1/2048*eL2(1)];
slopeH1=[1/2*eH1(1);1/4*eH1(1);1/8*eH1(1);1/16*eH1(1);1/32*eH1(1)];
slopeDGbroken=[1/2*eDGbroken(1);1/4*eDGbroken(1);1/8*eDGbroken(1);1/16*eDGbroken(1);1/32*eDGbroken(1)];%1/64*eDGbroken(1)];
%计算用于绘制斜率参考线的数据，slopeL2 和 slopeDGbroken 分别对应于 eL2 和 eDGbroken 的斜率序列，用这些斜率参考线来比较实际误差和理论斜率
semilogy(X,eL2,'b-s',X,eH1,'c-x',X,eDGbroken,'g-O',X,slopeL2,'k--',X,slopeH1,'y--',X,slopeDGbroken,'r--','LineWidth',1.9); %在图形窗口中绘制四条曲线，字符串参数分别指定数据样式和颜色（蓝色实线代表 eL2，绿色圆圈代表 eDGbroken，黑色虚线代表 slopeL2，红色虚线代表 slopeDGbroken)
% ,N,N.^(-5/2),'k--');
h2=legend('$||u-u_h||_0$','$||u-u_h||_1$','$||u-u_h||_{\star}$','$ch^2$','$ch$','$ch$','Location','southwest'); %在图形的左下角创建图例 h2，用来标识图中的曲线
set(h2,'Interpreter','latex')  %设置图例中字符串的解释器为 LaTeX，以便支持 LaTeX 格式的文本
set(h2,'fontsize',15);  %设置图例中字符串的字体大小为15
set(gca,'XTick',1:L,'XLim',[0,L+1],'Ylim',[1e-5 1e2]);  %设置图形的横轴刻度，包含从1到 L 的整数
xlabel('refined iteration');
ylabel('log(error)');
grid on;
set(gca,'GridLineStyle','--'); %设置网格线样式
set(gca,'MarkerSize',6);  %增加标记大小
%title('');
%grid on;

