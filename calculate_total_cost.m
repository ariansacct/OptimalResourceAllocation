function [ total_cost,cost ] = calculate_total_cost( p,m,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cost=-1*ones(size(p,1),1);
total_cost=-1*ones(size(p,1),1);

% lambda_hw = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
% lambda_vmm = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
% processor_s = 100 + (200 - 100) .* rand(1,n);
% processor_v = 1 + (50 - 1) .* rand(1,m);
% mem_s = 100 + (200 - 100) .* rand(1,n);
% mem_v = 1 + (50 - 1) .* rand(1,m);
% st_s = 100 + (200 - 100) .* rand(1,n);    %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE
% st_v = 1 + (50 - 1) .* rand(1,m);         %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE

global lambda_hw;
global lambda_vmm;
global processor_s;
global processor_v;
global mem_s;
global mem_v;
global st_s;
global st_v;

for row=1:size(p,1)
    argument=p(row,:);
    x = zeros(m,n);
    for i=1:m
        x(i,argument(i))=1;
    end
%     disp(x);
    omega=any(x);
    
    
    c=0;
    for j=1:n
        c=c + omega(j)*(lambda_hw(j)+lambda_vmm(j));
    end
%     disp(c);
    
    p_p=0;
    for j=1:n
        p_p = p_p + max( 0 , dot(x(:,j),processor_v)-processor_s(j) );
    end
    
    p_mem=0;
    for j=1:n
        p_mem = p_mem + max( 0 , dot(x(:,j),mem_v)-mem_s(j) );
    end
    
    p_st=0;
    for j=1:n
        p_st = p_st + max( 0 , dot(x(:,j),st_v)-st_s(j) );
    end
    
    eta = 0.01;
    tc=c+eta*(p_p+p_mem+p_st);
    
    cost(row,:)=c;
    total_cost(row,:)=tc;
    
end

return

