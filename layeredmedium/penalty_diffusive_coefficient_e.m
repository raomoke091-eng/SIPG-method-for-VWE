function r=penalty_diffusive_coefficient_e(x,y,t)
%% ע����������ǲ�ͬ��Beatrice���ϵģ��ڷ���Ҳ�������ɢϵ�������Σ�������Grote�Ĳ����̽ⷨ��
    a_1 = 1;
    a_2 = 2.5;
    a_3 = 2;
    a_4 = 10;
    transition_width = 100.0;  % ��������Ŀ��

    % ���� sigmoid ����
    sigmoid = @(t) 1 / (1 + exp(-t));
    r = a_1 + (a_2 - a_1) * sigmoid(transition_width * (y - 0.4)) ...
            + (a_3 - a_2) * sigmoid(transition_width * (y - 0.6)) ...
            + (a_4 - a_3) * sigmoid(transition_width * (y - 0.8));