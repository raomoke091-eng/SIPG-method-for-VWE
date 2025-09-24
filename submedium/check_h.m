function check
% An example: Backward Euler, NIPG, P1(E) basis function, ordinary penalty type.
% One can modify following parameters for Backward Euler(theta=1), Forward Euler(theta=0), Crank-Nicolson discretizations(theta=1/2), 
%                                         different schemes such as NIPG(eipsilon=1), IIPG(eipsilon=0), SIPG(eipsilon=-1), OBB(eipsilon=1,sigma=0,k>=2),
%                                         basis functions of different basis types(basis_type=1,2,3...), 
%                                         penalty(beta=1) and super-penalty(beta=3) methods
% with appropriate proportional relation between time step delta_t and mesh size h (i.e., 1/2^Iter_refinement)!

basis_type=2;
eipsilon=-1;
sigma=800;
beta=1;
theta=1;
T=0.5;
Iter_refinement=4;N_T=2^(Iter_refinement*2);
[L2_error_4,H1_error_4]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);
%{
%  Iter_refinement=6;N_T=2^(Iter_refinement*2);
%  [L2_error_6,H1_error_6]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%这里给出primal格式的误差表：
%%
h=[1/2;1/4;1/8;1/16]%;1/32]%;1/64];
eL2=[L2_error_1;L2_error_2;L2_error_3;L2_error_4]%;L2_error_5]%;L2_error_6];%%%%%%%%%%%%%%n维向量
eH1=[H1_error_1;H1_error_2;H1_error_3;H1_error_4]%;H1_error_5]%;H1_error_6];
%eInf=[uh_infinity_error_8;uh_infinity_error_16;uh_infinity_error_32;uh_infinity_error_64];%%%%%求节点最大误差(无穷范数)



L=length(h);
eL2Order=zeros(L-1,1);%%%%%%%%%%%%n-1维向量，因为n步计算只能计算n-1步收敛率
eH1Order=zeros(L-1,1);
eInfOrder=zeros(L-1,1);

for i=1:L-1
    eL2Order(i)=log(eL2(i)/eL2(i+1))/log(2);
    eH1Order(i)=log(eH1(i)/eH1(i+1))/log(2);
%  eInfOrder(i)=log(eInf(i)/eInf(i+1))/log(2);
end

T1=table(h,eL2,eH1);
disp('the errors for example 1:');
format short e;
disp(T1);
h=h(2:L,1);
disp('the convergence rate for example 1:');
format short;
T2=table(h,eL2Order,eH1Order);
disp(T2);
disp('p=1,N_T=2^(Iter_refinement*2),sigma=200');
%}


