function xdot = NREL5MW_RHSFunc(x,u,M,C,K)

pos = x(1:length(x)/2);
vel = x(length(x)/2+1:end);

acc = M\(u - C*vel - K*pos);

xdot = [vel;acc];


end

