function [y] = h(xk) % single-fidelity function
%function [y] = h(xk_exp,xk_imrv2) % multi-fidelity function

%returns radii (observed variable) from ensemble

y = [xk(1,:)]; % single-fidelity y
%y = [xk_exp(1,:),] % multi-fidelity y

%y = xk(1:4,:);
%y = [xk(1,:);xk(4,:)];
end

