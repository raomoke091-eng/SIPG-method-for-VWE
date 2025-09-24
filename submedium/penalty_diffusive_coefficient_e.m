function r=penalty_diffusive_coefficient_e(x,y,t)
%% ע����������ǲ�ͬ��Beatrice���ϵģ��ڷ���Ҳ�������ɢϵ�������Σ�������Grote�Ĳ����̽ⷨ��
sigmoid = @(t) 1 / (1 + exp(-t));
    % ���峣��
    a_1 = 1.0;
    a_2 = 5.0;
    a_3 = 10.0;
    
    % ������ɿ��
    transition_width1 = 100.0;
    transition_width2 = 100.0;
    transition_width3 = 100.0;
    transition_width4 = 100.0;
    transition_width5 = 100.0;
    
    % Ϊÿ����������ֵ
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
