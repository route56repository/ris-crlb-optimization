function [x,comp_x] = generate_random_prs_matrix(n, m, p)
x=sqrt(p/(n*m))*hadamard(m);
x=x(randperm(m),:);
x(n+1:end,:)=[];
comp_x=x*x';


end
