function r=diffusive_coefficient(x,y,t)

%{
if 0 <= y < 0.4
    r = 1;
elseif 0.4 <= y <0.6
    r = 2.5;
elseif 0.6 <= y <0.8
    r = 2;
elseif 0.8 <= y <=1
    r = 10;
end
%}
% 定义常数
    a_1 = 1;
    a_2 = 2.5;
    a_3 = 2;
    a_4 = 10;
    transition_width = 100.0;  % 过渡区域的宽度

    % 计算 sigmoid 函数
    sigmoid = @(t) 1 / (1 + exp(-t));
    r = a_1 + (a_2 - a_1) * sigmoid(transition_width * (y - 0.4)) ...
            + (a_3 - a_2) * sigmoid(transition_width * (y - 0.6)) ...
            + (a_4 - a_3) * sigmoid(transition_width * (y - 0.8));