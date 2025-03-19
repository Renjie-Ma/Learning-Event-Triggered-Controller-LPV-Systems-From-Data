function [r] = getRef_RtEv_8(tmax,R,pi)
for k = 1:(tmax-1)
    theta = 2 * pi * k / (tmax-1); % 角度参数
    x_ref = R * sin(theta);        % x坐标
    y_ref = R * sin(2 * theta);    % y坐标
    r{k} = [x_ref; y_ref];         % 存储当前点的坐标
end
end