function [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference,left_lower_point,h_partition)
%Generate the Gauss coefficients and Gauss points on the local retangular element by using affine tranformation
%Gauss_coefficient_local,Gauss_point_local:the Gauss coefficients and Gauss points on the local retangular element
%Some details are in my standard 2-dimension FE Tool boxes notes part 6
%Gauss_coefficient_reference,Gauss_point_reference: the Gauss coefficients and Gauss points on the reference square 
%left_lower_point: The coordinates of the left-lower vertice of the current element whose local basis we are evaluating
%h_partition: the step size of the partition

Jacobi=h_partition(1)*h_partition(2)/4;
Gauss_coefficient_local=Gauss_coefficient_reference*Jacobi;
Gauss_point_local(:,1)=h_partition(1)*Gauss_point_reference(:,1)/2+left_lower_point(1)+h_partition(1)/2;
Gauss_point_local(:,2)=h_partition(2)*Gauss_point_reference(:,2)/2+left_lower_point(2)+h_partition(2)/2;