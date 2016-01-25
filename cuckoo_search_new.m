
function [bestnest,fmin]=cuckoo_search_new(n)
if nargin<1,
    % Number of nests (or different solutions)
    n=100;
end

% Discovery rate of alien eggs/solutions
pa=0.25;

global M;
global N;
global cost;
cost = 10000;   % initialize to a big number
global best_cost;
best_cost=50000;    % initialize to a big number
global x;
M = 10;
N = 3;
low = 1;
up = N;

%% Change this if you want to get better results
N_IterTotal=1000;
% N_IterTotal=3;
Iter=0;
%% Simple bounds of the search domain
% Lower bounds
nd=M; % same as dim
Lb=low*ones(1,nd);
% Upper bounds
Ub=up*ones(1,nd);

% Random initial solutions
nest=ones(n,nd);
for i=1:n,
    %nest(i,:)=floor( Lb+(Ub-Lb).*rand(size(Lb)) )
    nest(i,:)=round( Lb+(Ub-Lb).*rand(size(Lb)) );
end
determine_x(nest);

% Get the current best
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness);

N_iter=0;
%% Starting iterations
for iter=1:N_IterTotal,
    % Generate new solutions (but keep the current best)
    new_nest=get_cuckoos(nest,bestnest,Lb,Ub);
    
    %determine_x(new_nest);
    
    
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter
    N_iter=N_iter+n;
    % Discovery and randomization
    new_nest=empty_nests(nest,Lb,Ub,pa);
    
    %determine_x(new_nest);
    
    % Evaluate this set of solutions
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter again
    N_iter=N_iter+n;
    % Find the best objective so far
    if fnew<fmin,
        fmin=fnew;
        bestnest=best;
        best_cost=cost;
    end
    
    Iter=Iter+1;
    Iter
    
end %% End of iterations

%% Post-optimization processing
%% Display all the nests
disp(strcat('Total number of iterations=',num2str(N_iter)));
fmin
best_cost
bestnest

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

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness)
% Evaluating all new solutions
for j=1:size(nest,1),
    fnew=fobj(newnest(j,:));
    if fnew<=fitness(j),
        fitness(j)=fnew;
        nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=nest(K,:);

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

%% You can replace the following by your own functions
% A d-dimensional objective function
function z=fobj(sol)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2.
%  with a minimum at (1,1, ...., 1);
% z=sum((u-1).^2);
x = determine_x(sol);
global M;
global N;
global cost;
M = 10;
N = 3;
PHR = 0.00005 + (0.00010 - 0.00005) .* rand(1,N);
E = 15 + (25 - 15) .* rand(M,N);
CBW = 1 + (4 - 1) .* rand(N,N);
CR = 0 + (25 - 0) .* rand(M,M);
PL = 100 + (200 - 100) .* rand(N,N);
CHR = 0.00015 + (0.00030 - 0.00015) .* rand(N,N);
Mem = 100 + (200 - 100) .* rand(1,N);
mem = 1 + (50 - 1) .* rand(1,M);
L = 100 + (200 - 100) .* rand(1,N);
l = 1 + (50 - 1) .* rand(1,M);

eta = 0.01;    %   this should be set appropriately

C1 = 0;
C2 = 0;
for k=1:N
    for i=1:M
        C1 = C1 + PHR(1,k) .* x(i,k) .* CR(i,k);
    end
end

for p = 1:N-1
    for q = p+1:N
        for i = 1:M-1
            for j = i+1:M
                C2 = C2 + CHR(p,q) .* x(i,p) .* x(j,q) .* CR(i,j) / CBW(p,q);
            end
        end
    end
end

C = C1 + C2;

Pm = 0;
for k = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + mem(1,i) .* x(i,k) - Mem(1,k);
    end
    Pm = Pm + max(0, temp);
end

Pl = 0;
for K = 1:N
    temp = 0;
    for i = 1:M
        temp = temp + l(1,i) .* x(i,k) - L(1,k);
    end
    Pl = Pl + max(0, temp);
end

Pc = 0;
for p = 1:N
    for q = p+1:N
        
        temp = 0;
        for i = 1:M
            for j = 1:M
                temp = temp + x(i,p) .* x(j,q) .* CR(i,j);
            end
        end
        maxsecondarg = temp - PL(p,q);
        
        
        
        
        Pc = Pc + max(0, maxsecondarg);
        
    end
end

cost=C;
z = C + eta .* (Pm + Pl + Pc);

