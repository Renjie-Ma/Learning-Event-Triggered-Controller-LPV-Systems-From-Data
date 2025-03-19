function [r] = getRef_RtEv1(tmax)
for k=1:1:(tmax-1)
    if k<(tmax-1)/4
        r{k}=-1;
    elseif k>=(tmax-1)/4 && k<(tmax-1)/2
        r{k}=1;
    elseif k>=(tmax-1)/2 && k<3*(tmax-1)/4
        r{k}=-1;
    elseif k>=3*(tmax-1)/4 && k<=(tmax-1)
        r{k}=1;
    end
end
end

