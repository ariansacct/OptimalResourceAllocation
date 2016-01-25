
% function [bestnest,fmin]=cuckoo_search_new2(n)
function [bestnest,fmin]=cuckoo_search_new2()
% if nargin<1,
%     % Number of nests (or different solutions)
%     n=100;
% end

noofnests=100;

% Discovery rate of alien eggs/solutions
pa=0.25;

global m;
global n;
m = 100;
n = 30;

global lambda_hw;
global lambda_vmm;
global processor_s;
global processor_v;
global mem_s;
global mem_v;
global st_s;
global st_v;

lambda_hw = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
lambda_vmm = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
processor_s = 100 + (200 - 100) .* rand(1,n);
processor_v = 1 + (50 - 1) .* rand(1,m);
mem_s = 100 + (200 - 100) .* rand(1,n);
mem_v = 1 + (50 - 1) .* rand(1,m);
st_s = 100 + (200 - 100) .* rand(1,n);    %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE
st_v = 1 + (50 - 1) .* rand(1,m);         %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE


% global cost;
% cost = 10000;   % initialize to a big number
% global best_cost;
% best_cost=50000;    % initialize to a big number
% global X;

low = 1;
up = n;

%% Change this if you want to get better results
% N_IterTotal=1000;
N_IterTotal=10;
Iter=0;
%% Simple bounds of the search domain
% Lower bounds
nd=m; % same as dim
Lb=low*ones(1,nd);
% Upper bounds
Ub=up*ones(1,nd);

% Random initial solutions
nest=ones(noofnests,nd);
for i=1:noofnests,
    %nest(i,:)=floor( Lb+(Ub-Lb).*rand(size(Lb)) )
    nest(i,:)=round( Lb+(Ub-Lb).*rand(size(Lb)) );
end
determine_x(nest);

% Get the current best
fitness=10^10*ones(noofnests,1);
cost=10^10*ones(noofnests,1);
[cmin,fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,cost);

N_iter=0;
%% Starting iterations
for iter=1:N_IterTotal,
    % Generate new solutions (but keep the current best)
    new_nest=get_cuckoos(nest,bestnest,Lb,Ub);
    
    determine_x(new_nest);
    
    
    [cnew,fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,cost);
    % Update the counter
    N_iter=N_iter+noofnests;
    % Discovery and randomization
    new_nest=empty_nests(nest,Lb,Ub,pa);
    
    determine_x(new_nest);
    
    % Evaluate this set of solutions
    [cnew,fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,cost);
    % Update the counter again
    N_iter=N_iter+noofnests;
    % Find the best objective so far
    if fnew<fmin,
        fmin=fnew;
        cmin=cnew;
        bestnest=best;
%         best_cost=cost;
    end
    
    Iter=Iter+1;
    Iter
    
end %% End of iterations

%% Post-optimization processing
%% Display all the nests
disp(strcat('Total number of iterations=',num2str(N_iter)));
fmin
cmin
bestnest
end





%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    
    % In the next equation, the difference factor (s-best) means that
    % when the solution is the best solution, it remains unchanged.
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale;
    % otherwise, Levy flights may become too aggresive/efficient,
    % which makes new solutions (even) jump out side of the design domain
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
    % Apply simple bounds/limits
    nest(j,:)= floor( simplebounds(s,Lb,Ub) );
end
end

%% Find the current best nest
function [cmin,fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,cost)
% Evaluating all new solutions
for j=1:size(nest,1),
    [fnew,cnew]=fobj(newnest(j,:));
    if fnew<=fitness(j),
        fitness(j)=fnew;
        cost(j)=cnew;
        nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=nest(K,:);
cmin=cost(K,:);
end

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then
% this cuckoo's egg is less likely to be discovered, thus the fitness should
% be related to the difference in solutions.  Therefore, it is a good idea
% to do a random walk in a biased way with some random step sizes.
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
    new_nest(j,:)= floor( simplebounds(s,Lb,Ub) );
end
end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;
end

%% You can replace the following by your own functions
% A d-dimensional objective function
function [ total_cost,cost ]=fobj(p)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2.
%  with a minimum at (1,1, ...., 1);
% z=sum((u-1).^2);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% cost=-1*ones(size(p,1),1);
% total_cost=-1*ones(size(p,1),1);


% lambda_hw = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
% lambda_vmm = 0.00005 + (0.00010 - 0.00005) .* rand(1,n);
% processor_s = 100 + (200 - 100) .* rand(1,n);
% processor_v = 1 + (50 - 1) .* rand(1,m);
% mem_s = 100 + (200 - 100) .* rand(1,n);
% mem_v = 1 + (50 - 1) .* rand(1,m);
% st_s = 100 + (200 - 100) .* rand(1,n);    %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE
% st_v = 1 + (50 - 1) .* rand(1,m);         %   FELAN GOZASHTIMESH MESLE MEMORY; IN BAYAD AVAZ SHE

global m;
global n;

global lambda_hw;
global lambda_vmm;
global processor_s;
global processor_v;
global mem_s;
global mem_v;
global st_s;
global st_v;

% for row=1:size(p,1)
%     argument=p(row,:);
argument=p;
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
    
%     cost(row,:)=c;
%     total_cost(row,:)=tc;

cost=c;
total_cost=tc;
    
% end


end

function omega=findomega(j)
global X;
omega=max( max( 0,X(:,j) ) );
end

