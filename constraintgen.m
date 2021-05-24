function [A,b]=constraintgen(dim, predmod, x0)

%Input Constraints
A = [ eye(dim.N) ; -1*eye(dim.N) ];
b = [20*ones(dim.N,1) ; zeros(dim.N,1)];

%State Constraints
A = [A ; predmod.S ; -predmod.S];
b = [b ; 7.2*ones(size(predmod.S,1),1)-predmod.T*x0 ; predmod.T*x0];


end