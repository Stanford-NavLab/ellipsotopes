%% user parameters
% center
c = [0; 0];
% generator matrix 
G = [1 0 0 0;
     0 1 0 0] ;
% %  % constraints
% A = [1 0 -1  0;
%      0 1  0 -1] ;
% b = [1; 0];
% tol = 1e-2 ;
% % index set
% J = {[1,2],[3,4]};

% norm to consider
norm_p = 2 ;

E1 = ellipsotope(norm_p,c,G);
plot(E);
axis equal

 
%% automated from here
% get number of generators


%% plotting
