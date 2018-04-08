% Mann-Kendallͻ�����
% ��ة
% 2018/4/6

% tsΪʱ������
% tsDataΪ��ts��Ӧ�Ĵ�����Mann-Kendall���Ƽ��������
% figureTitleΪ��ͼ���⣨ʡ�Լ�����ͼ��

% recordΪ2�о���
% ��һ��Ϊ�����򼴽�����ͻ�䣨���㣩��ʱ��㣬�ڶ���Ϊ����������
% UF��UBΪ�������ߵ�ֵ

function [record, UF, UB] = MannKendall_mutation(ts, tsData, figureTitle)
    if nargin < 3
        figureTitle = '';
    end

    %% Mann-Kendallͻ��������岿��
    len = length(tsData);
    s = zeros(len, 1);
    UF = s;
    tsData_reverse = s;
    s_reverse = s;
    UB = s;
    tsData_reverse(:) = tsData(end : -1 : 1);  % ��������
    for k = 2 : len
        s(k) = sum(tsData(k) > tsData(1 : k)) + s(k-1);
        s_reverse(k) = ...
        sum(tsData_reverse(k) > tsData_reverse(1 : k)) + s_reverse(k-1);
        E = k * (k - 1) / 4;  % s(k)�ľ�ֵ
        Var = k * (k - 1) * (2 * k + 5) / 72;  % s(k)�ķ���
        UF(k) = (s(k) - E) / sqrt(Var);
        UB(len + 1 -k) = -(s_reverse(k) - E) / sqrt(Var);
    end

    %% �жϽ���λ��
    d = UF - UB;
    record = [];
    crossPoint = [];
    for k = 1 : len
        if d(k) == 0
            record = [record; ts(k), UF(k)];
            crossPoint = [crossPoint; k, UF(k)];
        elseif k == len
            break
        elseif d(k) * d(k + 1) < 0
            A = [UF(k+1)-UF(k), -1
                UB(k+1)-UB(k), -1];
            b = [k*(UF(k+1)-UF(k)) - UF(k)
                k*(UB(k+1)-UB(k)) - UB(k)];
            X = A\b;  % ���Ԫһ�η�����
            record = [record; ts(k), X(2)];
            crossPoint = [crossPoint; X(1), X(2)];
        end
    end
    if ~isempty(record)
        record(abs(record(:, 2))>1.96, :) = [];  % ɾ����ͻ���Ľ���
        crossPoint(abs(crossPoint(:, 2))>1.96, :) = [];% ɾ����ͻ���Ľ���
    end    
    %% ��ͻ����ͼʱ��ʹ��UF��UB_reverse
    if isempty(figureTitle)
        return
    end
    figure()
    hold on
    plot(UF, 'r-', 'linewidth', 2)
    plot(UB, 'b-.', 'linewidth', 2)
    plot([1, len, nan, len, 1], [1.96, 1.96, nan, -1.96, -1.96], ...
        'k:', 'linewidth', 1)
    legend('UF', 'UB', '\alpha=0.05', 'AutoUpdate', 'off')
    legend('boxoff')
    plot([1, len], [0, 0], 'k', 'LineWidth', 0.5)
    axis([1, len, -3, 6])
    xticks(1 : len)
    xticklabels(ts)
    xlabel('���', 'FontName', 'TimesNewRoman', 'FontSize', 12)
    ylabel('MKͳ����', 'FontName', 'TimesNewRoman', 'Fontsize', 12)
    title(figureTitle)
    plot(crossPoint(:, 1), crossPoint(:, 2), 'go', 'MarkerFaceColor', 'g')
    for k = 1 : size(crossPoint, 1)
        plot([1 1]*crossPoint(k, 1), [-3 crossPoint(k, 2)], ...
            'k--', 'linewidth', 1)
    end
end
        