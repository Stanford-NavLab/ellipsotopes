classdef ellipsotope_box_world < handle
    %ellipsotope_world 
    %   World with box obstacles represented as ellipsotopes, and start and
    %   goal position for agent
    
    properties
        bounds
        n_obs
        n_dim
        start
        goal
        obs 
        buffer = 1;
        w_obs_min = 0.5 ; % minimum obstacle width [m] 
        w_obs_max = 1.0 ; % maximum obstacle width [m]
    end
    
    methods
        %% constructor
        function W = ellipsotope_box_world()
            % default properties
            W.bounds = 5.*[-1,1,-1,1] ; % 2-D world
            W.n_obs = 10; % number of obstacles
            W.n_dim = 2;
        end
        
        %% setup
        function setup(W)
            
            % generate uniformly-distributed random obstacles
            O_ctr = rand_in_bounds(W.bounds,[],[],W.n_obs) ; % obstacle positions
            O_wid = rand_range(W.w_obs_min,W.w_obs_max,[],[],2,W.n_obs) ; % obstacle widths

            % obstacle zonotope array
            W.obs = cell(1,W.n_obs) ;
            for i = 1:W.n_obs
                % create generators from widths and randomly rotate 
                O_gen = O_wid(:,i) .* eye(n_dim) ;
                O_gen = rotation_matrix_2D(pi*rand()) * O_gen;
                W.obs{i} = ellipsotope(2,O_ctr(:,i),O_gen,[],[],{1,2}) ;
            end
            
            % generate start position on left side of room with initial
            % heading of 0, and make sure it's not too close to the walls
            B = W.bounds ; b = W.buffer ;
            xlo = B(1) ; xhi = B(2) ; ylo = B(3) ; yhi = B(4) ;
            
            xlo = xlo + 2*b ;
            xhi = xhi - 2*b ;
            ylo = ylo + 2*b ;
            yhi = yhi - 2*b ;
            
            if isempty(W.start)
                s = [xlo ;
                     rand_range(ylo, yhi) ;
                     0 ] ;
                W.start = s ;   
            end
            
            % generate goal position on right side of room
            if isempty(W.goal)
                g = [xhi ;
                     rand_range(ylo, yhi)] ;
                W.goal = g ;
            end
        end
        
        function world_info = get_world_info(W)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

