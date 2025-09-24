function check_t
% An example: Backward Euler, NIPG, P1(E) basis function, ordinary penalty type.
% One can modify following parameters for Backward Euler(theta=1), Forward Euler(theta=0), Crank-Nicolson discretizations(theta=1/2), 
%                                         different schemes such as NIPG(eipsilon=1), IIPG(eipsilon=0), SIPG(eipsilon=-1), OBB(eipsilon=1,sigma=0,k>=2),
%                                         basis functions of different basis types(basis_type=1,2,3...), 
%                                         penalty(beta=1) and super-penalty(beta=3) methods
% with appropriate proportional relation between time step delta_t and mesh size h (i.e., 1/2^Iter_refinement)!

basis_type=2;
eipsilon=-1;
sigma=20;
beta=1;
theta=1;
T=0.1;
% Relation between delta_t and h
Iter_refinement=6;N_T=2;
[L2_error_1,H1_error_1,DG_broken_error_1]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

Iter_refinement=6;N_T=4;
[L2_error_2,H1_error_2,DG_broken_error_2]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

Iter_refinement=6;N_T=8;
[L2_error_3,H1_error_3,DG_broken_error_3]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

Iter_refinement=6;N_T=16;
[L2_error_4,H1_error_4,DG_broken_error_4]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%Iter_refinement=5;N_T=32;
%[L2_error_5,H1_error_5,DG_broken_error_5]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%Iter_refinement=5;N_T=64;
%[L2_error_6,H1_error_6,DG_broken_error_6]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%  Iter_refinement=6;N_T=2^(Iter_refinement*2);
%  [L2_error_6,H1_error_6]=DG_wave_triangulation(Iter_refinement,basis_type,eipsilon,sigma,beta,theta,N_T,T);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������primal��ʽ������
%%
delta_t = [1/2;1/4;1/8;1/16]%;1/32]%;1/64];
eL2=[L2_error_1;L2_error_2;L2_error_3;L2_error_4]%;L2_error_5]%;L2_error_6];%������
eH1=[H1_error_1;H1_error_2;H1_error_3;H1_error_4]%;H1_error_5]%;H1_error_6];
eDGbroken=[DG_broken_error_1;DG_broken_error_2;DG_broken_error_3;DG_broken_error_4]%;DG_broken_error_5]%;DG_broken_error_6];
% eInf=[uh_infinity_error_8;uh_infinity_error_16;uh_infinity_error_32;uh_infinity_error_64];%%%%%��ڵ�������(�����)


L=length(delta_t);
eL2_Order=zeros(L-1,1);%%%%%%%%%%%%n-1ά��������Ϊn������ֻ�ܼ���n-1��������
eH1_Order=zeros(L-1,1);
eDGbroken_Order=zeros(L-1,1);
% eInfOrder=zeros(L-1,1);

for i=1:L-1
    eL2_Order(i)=log(eL2(i)/eL2(i+1))/log(2);
    eH1_Order(i)=log(eH1(i)/eH1(i+1))/log(2);
    eDGbroken_Order(i)=log(eDGbroken(i)/eDGbroken(i+1))/log(2);
% eInfOrder(i)=log(eInf(i)/eInf(i+1))/log(2);
end

T1=table(delta_t,eL2,eDGbroken);
disp('the errors for example 1:');
format short e;
disp(T1);

delta_t_2=delta_t(2:L,1);
T2=table(delta_t_2,eL2_Order,eDGbroken_Order);
disp('the convergence rate for example 1:');
format short;
disp(T2);
disp('p=2,h=2^(-6),N_T=2^i,sigma=20');







