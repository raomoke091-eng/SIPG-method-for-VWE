function d=element_diameter(vertices)
x12=vertices(1,2)-vertices(1,1);
x23=vertices(1,3)-vertices(1,2);
y12=vertices(2,2)-vertices(2,1);
y23=vertices(2,3)-vertices(2,2);

cen1=(vertices(:,1)+vertices(:,2))/2;       %12点连边中点
cen2=(vertices(:,2)+vertices(:,3))/2;       %23点连边中点

%分情况(7)计算两条边的垂直平分线的交点坐标
if x12~=0 && y12~=0 && x23~=0 && y23~=0 
    k1=-1/((vertices(2,1)-vertices(2,2))/(vertices(1,1)-vertices(1,2)));    %12点连边垂直平分线
    b1=cen1(2)-k1*cen1(1);
    k2=-1/((vertices(2,2)-vertices(2,3))/(vertices(1,2)-vertices(1,3)));    %23点连边垂直平分线
    b2=cen2(2)-k2*cen2(1);
    x0=-(b1-b2)/(k1-k2);        %求两直线交点x坐标
    y0=-(-b2*k1+b1*k2)/(k1-k2); %求两直线交点y坐标
elseif x12~=0 && y12==0 && x23~=0 && y23~=0
    k2=-1/((vertices(2,2)-vertices(2,3))/(vertices(1,2)-vertices(1,3)));    %23点连边垂直平分线
    b2=cen2(2)-k2*cen2(1);
    x0=cen1(1);
    y0=k2*x0+b2;
elseif x12==0 && y12~=0 && x23~=0 && y23~=0
    k2=-1/((vertices(2,2)-vertices(2,3))/(vertices(1,2)-vertices(1,3)));    %23点连边垂直平分线
    b2=cen2(2)-k2*cen2(1);
    y0=cen1(2);
    x0=(y0-b2)/k2;
elseif x12~=0 && y12~=0 && x23==0 && y23~=0 
    k1=-1/((vertices(2,1)-vertices(2,2))/(vertices(1,1)-vertices(1,2)));    %12点连边垂直平分线
    b1=cen1(2)-k1*cen1(1);
    y0=cen2(2);
    x0=(y0-b1)/k1;        
elseif x12~=0 && y12==0 && x23==0 && y23~=0  %交点即为第二个顶点
    x0=vertices(1,2);
    y0=vertices(2,2);
elseif x12~=0 && y12~=0 && x23~=0 && y23==0 
    k1=-1/((vertices(2,1)-vertices(2,2))/(vertices(1,1)-vertices(1,2)));    %12点连边垂直平分线
    b1=cen1(2)-k1*cen1(1);
    x0=cen2(1);
    y0=k1*x0+b1; 
elseif x12==0 && y12~=0 && x23~=0 && y23==0
    x0=vertices(1,2);
    y0=vertices(2,2);
end
                                  
d=sqrt((y0-vertices(2,1))^2+(x0-vertices(1,1))^2); 


