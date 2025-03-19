function [r] = getRef_RtEv_sin(tmax,R,pi)
for k=1:1:(tmax-1)
    theta=2*pi*k/(tmax-1);
    x_ref= R*sin(4*theta);
    r{k}=x_ref;
end
end