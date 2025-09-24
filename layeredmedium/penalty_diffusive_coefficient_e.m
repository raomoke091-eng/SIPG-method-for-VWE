function r=penalty_diffusive_coefficient_e(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。
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