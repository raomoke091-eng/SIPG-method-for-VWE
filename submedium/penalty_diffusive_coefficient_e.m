function r=penalty_diffusive_coefficient_e(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。
sigmoid = @(t) 1 / (1 + exp(-t));
    % 定义常量
    a_1 = 1.0;
    a_2 = 5.0;
    a_3 = 10.0;
    
    % 定义过渡宽度
    transition_width1 = 100.0;
    transition_width2 = 100.0;
    transition_width3 = 100.0;
    transition_width4 = 100.0;
    transition_width5 = 100.0;
    
    % 为每个条件分配值
    if 0 <= x < 0.2
        r = a_1 + (a_2 - a_1) * sigmoid(transition_width1 * (y - 0.35)) + ...
                    (a_3 - a_2) * sigmoid(transition_width1 * (y - 0.85));
    elseif 0.2 <= x < 0.4
        r = a_1 + (a_2 - a_1) * sigmoid(transition_width2 * (y - 0.35)) + ...
                    (a_3 - a_2) * sigmoid(transition_width2 * (y + x - 1.05));
    elseif 0.4 <= x < 0.6
        r= a_1 + (a_2 - a_1) * sigmoid(transition_width3 * (y - 0.35)) + ...
                    (a_3 - a_2) * sigmoid(transition_width3 * (y - 0.65));
    
    elseif 0.6 <= x < 0.7
        r = a_1 + (a_2 - a_1) * sigmoid(transition_width4 * (y- 0.35)) + ...
                    (a_3 - a_2) * sigmoid(transition_width4 * (y - 2 * x + 0.55));
    
    elseif 0.7 <= x <= 1
        r = a_1 + (a_2 - a_1) * sigmoid(transition_width5 * (y - 0.35)) + ...
                    (a_3 - a_2) * sigmoid(transition_width5 * (y - 0.85));
    end
