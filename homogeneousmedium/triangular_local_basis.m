function r=triangular_local_basis(x,y,basis_type,derivative_degree_x,derivative_degree_y)
%本函数用途：计算三角形单元上的形状函数及其导数
% (x,y)表示要计算形函数的相关值的点的坐标
% basis_type表示形函数的类型/形函数最高项次数
% derivative_degree_x表示形函数关于x求导的阶数 关于x求几阶导
% derivative_degree_y表示形函数关于y求导的阶数 关于y求几阶导
% 返回值为向量r，包含全部形函数(注意这里的形函数全部取为单项式 参考Beatrice P35)

if basis_type==1 %形函数最高项次数为1，即形函数为线性二元多项式
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1 x y];  %返回值为在输入的点处三个形函数1,x,y的具体值
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0 1 0];
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0 0 1];   
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0 0 0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0 0 0]; 
    elseif derivative_degree_x==1&&derivative_degree_y==1
        r=[0 0 0]; 
    end
elseif basis_type==2
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1 x y x.^2 x.*y y.^2];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0 1 0 2.*x y 0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0 0 1 0 x 2.*y];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0 0 0 2 0 0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0 0 0 0 0 2];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0 0 0 0 1 0];
    end
elseif basis_type==3
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6.*x    2.*y     0        0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2.*x     6.*y];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0];
    end
elseif basis_type==4
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2];
    end
elseif basis_type==5
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3];
    end
elseif basis_type==6
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4];
    end
elseif basis_type==7
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6      x^7     x^6*y     x^5*y^2     x^4*y^3       x^3*y^4     x^2*y^5    x*y^6     y^7];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0        7*x^6   6*x^5*y   5*x^4*y^2   4*x^3*y^3     3*x^2*y^4   2*x*y^5    y^6       0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5    0       x^6       2*x^5*y     3*x^4*y^2     4*x^3*y^3   5*x^2*y^4  6*x*y^5   7*y^6];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0        0       6*x^5     10*x^4*y    12*x^3*y^2    12*x^2*y^3  10*x*y^4   6*y^5     0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0        42*x^5  30*x^4*y  20*x^3*y^2  12*x^2*y^3    6*x*y^4     2*y^5      0         0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4   0       0         2*x^5       6*x^4*y       12*x^3*y^2  20*x^2*y^3 30*x*y^4  42*y^5];
    end
elseif basis_type==8
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6      x^7     x^6*y     x^5*y^2     x^4*y^3       x^3*y^4     x^2*y^5    x*y^6     y^7      x^8     x^7*y      x^6*y^2      x^5*y^3       x^4*y^4      x^3*y^5      x^2*y^6      x*y^7     y^8];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0        7*x^6   6*x^5*y   5*x^4*y^2   4*x^3*y^3     3*x^2*y^4   2*x*y^5    y^6       0        8*x^7   7*x^6*y    6*x^5*y^2    5*x^4*y^3     4*x^3*y^4    3*x^2*y^5    2*x*y^6      y^7       0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5    0       x^6       2*x^5*y     3*x^4*y^2     4*x^3*y^3   5*x^2*y^4  6*x*y^5   7*y^6    0       x^7        2*x^6*y      3*x^5*y^2     4*x^4*y^3    5*x^3*y^4    6*x^2*y^5    7*x*y^6   8*y^7];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0        0       6*x^5     10*x^4*y    12*x^3*y^2    12*x^2*y^3  10*x*y^4   6*y^5     0        0       7*x^6      12*x^5*y     15*x^4*y^2    16*x^3*y^3   15*x^2*y^4   12*x*y^5     7*y^6     0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0        42*x^5  30*x^4*y  20*x^3*y^2  12*x^2*y^3    6*x*y^4     2*y^5      0         0        56*x^6  42*x^5*y   30*x^4*y^2   20*x^3*y^3    12*x^2*y^4   6*x*y^5      2*y^6        0         0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4   0       0         2*x^5       6*x^4*y       12*x^3*y^2  20*x^2*y^3 30*x*y^4  42*y^5   0       0          2*x^6        6*x^5*y       12*x^4*y^2   20*x^3*y^3   30*x^2*y^4   42*x*y^5  56*y^6];
    end
end
