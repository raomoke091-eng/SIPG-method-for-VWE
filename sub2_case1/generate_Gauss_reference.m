function [Gauss_coefficient_reference,Gauss_point_reference]=generate_Gauss_reference(Gauss_point_number)
%Generate the Gauss coefficients and Gauss points on the reference square whose vertices are (-1,-1),(1,-1),(1,1),(-1,1)
%Gauss_point_number:the number of Gauss points in the formula. The Gauss formula depends on it.
%Gauss_coefficient_reference,Gauss_point_reference: the Gauss coefficients and Gauss points on the reference sqrare
if Gauss_point_number==4
    Gauss_coefficient_reference=[1,1,1,1];
    Gauss_point_reference=[1/sqrt(3),1/sqrt(3);1/sqrt(3),-1/sqrt(3);-1/sqrt(3),1/sqrt(3);-1/sqrt(3),-1/sqrt(3)];
elseif Gauss_point_number==9
    Gauss_coefficient_reference=[64/81,100/324,100/324,100/324,100/324,40/81,40/81,40/81,40/81];
    Gauss_point_reference=[0,0;sqrt(3/5),sqrt(3/5);sqrt(3/5),-sqrt(3/5);-sqrt(3/5),sqrt(3/5);-sqrt(3/5),-sqrt(3/5);0,sqrt(3/5);0,-sqrt(3/5);sqrt(3/5),0;-sqrt(3/5),0];
end