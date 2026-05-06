function color = rainbow
% 生成绿色在中间的明亮彩虹色 colormap (红->黄->绿->青->蓝)
n = 256; % 颜色数量（确保中间点是绿色）

% 定义关键颜色点：红 -> 橙 -> 黄 -> 绿(正中间) -> 青 -> 蓝
% 关键：绿色必须在 50% 位置
key_colors = [
    1.0, 0.25, 0.25;  % 1. 明亮红 (位置 0%)
    1.0, 0.65, 0.15;  % 2. 橙 (位置 25%)
    0.9, 0.9, 0.2;    % 3. 亮黄 (位置 37.5%)
    0.15, 1.0, 0.25;  % 4. 纯绿 (位置 50% - 正中间!)
    0.15, 1.0, 0.9;   % 5. 青 (位置 62.5%)
    0.25, 0.55, 1.0;  % 6. 明亮蓝 (位置 75%)
    0.2, 0.35, 0.95   % 7. 深蓝端 (位置 100%，但不偏暗)
];

% 关键位置定义（归一化 0-1）
key_positions = [0, 0.25, 0.375, 0.5, 0.625, 0.75, 1];

% 插值生成完整 colormap
rainbow_green_center = zeros(n, 3);
for i = 1:n
    t = (i-1) / (n-1); % 当前位置 0-1
    
    % 找到当前位置所在的区间
    idx = find(t >= key_positions, 1, 'last');
    if idx >= length(key_positions)
        rainbow_green_center(i, :) = key_colors(end, :);
    else
        % 计算区间内的比例
        t_local = (t - key_positions(idx)) / (key_positions(idx+1) - key_positions(idx));
        % 线性插值
        rainbow_green_center(i, :) = (1-t_local) * key_colors(idx, :) + ...
                                      t_local * key_colors(idx+1, :);
    end
end

color = flipud(rainbow_green_center) ;
end