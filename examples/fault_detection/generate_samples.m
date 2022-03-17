%% generate set of random samples for the fault detection examples
%
% use to compare conzono and ellisotope versions. use conzono sampling
% method

rng(1);

N_sims = 10; % number of simulations to run
N = 100; % number of iterations 
n = 2; % system dimension

% noise sets
W = conZonotope([0;0],eye(2));
V = conZonotope([0;0],[0.06 0; 0 0.6]);

% initial set of states
X0 = conZonotope([0.6;70],[0.06 0; 0 0.6]);

x0 = randPoint(X0,N_sims);
w = zeros(N_sims,n,N);
v = zeros(N_sims,n,N+1);

for i = 1:N_sims
    w(i,:,:) = randPoint(W,N);
    v(i,:,:) = randPoint(V,N+1);
end

samples.x0 = x0;
samples.w = w;
samples.v = v;

save('samples.mat','samples');
