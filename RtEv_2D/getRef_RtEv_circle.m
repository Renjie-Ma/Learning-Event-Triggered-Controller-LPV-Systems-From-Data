function [r] = getRef_RtEv_circle(tmax,R,pi)
for k=1:1:(tmax-1)
    theta=2*pi*k/(tmax-1);
    x_ref= R*sin(theta);
    y_ref= R*cos(theta);
    r{k}=[x_ref;y_ref];
end
end