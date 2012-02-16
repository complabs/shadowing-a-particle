%% Lab1 - Shadowing a Particle, Reference Trajectory Solver
% Implements lab1_odeSolver function that simulates trajectory for the single particle 
% (without any visualizations) and tracks total energy of the system. Simulation can be
% performed using different numerical integrator algorithms with different time-step
% sizes.
%
% Filename: lab1_odeSolver.m
% Date:     2012-02-11
% Author:   Mikica B Kocic 

function result = lab1_odeSolver( solver, N_pts, relError, h )

%% Arguments:
%
%  solver    -- ode45, rk4, forwardEuler or semiEuler
%  N_pts     -- number of points (solution) to return
%  relError  -- relative error tolerance (used by ode45)
%  h         -- integration time-step (used by RK4, ForwardEuler and SemiEuler)
%
% Result: returns matrix Nx6 in the following format:
%
%  [  t_0   x_0  y_0   vx_0  vy_0   E_tot_0  ;
%     ...
%     t_f   x_f  y_f   vx_f  vy_f   E_tot_f  ]

%=========================================================================================
%% Default arguments

    if exist( 'solver', 'var' ) == 0
        solver = 'RK4';               % Default ode solver
    end

    if exist( 'N_pts', 'var' ) == 0
        N_pts = 80;                   % Default number of points returned in result
    end

    if exist( 'relError', 'var' ) == 0
        relError = 1e-7;              % Default relative error (for ode45)
    end

    if exist( 'h', 'var' ) == 0
        h = 1e-3;                     % Default time-step
    end

%=========================================================================================
%% Definitions / Physical model parameters

    %%% Physical constants

    g = [ 0, -9.81 ];  % Acceleration of Gravity, m/s^2

    %%% Particle system parameters

    m    = 1.0;        % Particle mass, kg
    L    = 5.0;        % Characteristic length of the system, m
    x_t0 = [ -L, -L ]; % Particle initial position, m
    v_t0 = [  5, 10 ]; % Particle intital velocity, m/s
    
    %%% Parameters of the spherically symmetric attractive force fields

    fk = [  32,   40,  28,  16,  20 ]'     ;  % fk(k)   = strength of the force 'k'
    rk = [  0.3, 0.2, 0.4, 0.5, 0.3 ]' * L ;  % rk(k)   = radius of the force 'k'
    pk = [ -0.2,  0.8 ;  ... % pk1            % pk(k,:) = center position of the force 'k'
           -0.3, -0.8 ;  ... % pk2
           -0.6,  0.1 ;  ... % pk3
            0.4,  0.7 ;  ... % pk4
            0.8, -0.3    ... % pk5
         ] * L;

    %%% Integration parameters

    t_0  = 0;  % Initial time, s
    t_f  = 3;  % Final (simulation) time, s

%=========================================================================================
%% Solve the differential equation
% See nested functions: MyRK4(), MyForwardEuler() and MySemiEuler()

    if strcmpi( solver, 'RK4' ) % Use our 4-th order Runge-Kutta algorithm

        [ T, Y ] = MyRK4( ...
            @odeFunc,                    ... % handle that evaluates the equations
            t_0, t_f,                    ... % the interval of integration
            [ x_t0, v_t0 ],              ... % initial conditions
            h,                           ... % integration time-step
            N_pts                        ... % number of points in the result
        );
    
    elseif strcmpi( solver, 'ForwardEuler' ) % Use our forward Euler integration

        [ T, Y ] = MyForwardEuler( ...
            @odeFunc,                    ... % handle that evaluates the equations
            t_0, t_f,                    ... % the interval of integration
            [ x_t0, v_t0 ],              ... % initial conditions
            h,                           ... % integration time-step
            N_pts                        ... % number of points in the result
        );
    
    elseif strcmpi( solver, 'SemiEuler' ) % Use our semi-implicit Euler integration

        [ T, Y ] = MySemiEuler( ...
            @odeFunc,                    ... % handle that evaluates the equations
            t_0, t_f,                    ... % the interval of integration
            [ x_t0, v_t0 ],              ... % initial conditions
            h,                           ... % integration time-step
            N_pts                        ... % number of points in the result
        );
    
    elseif strcmpi( solver, 'ode45' ) % Use MATLAB Runge-Kutta (4,5) formula
        
        [ T, Y ] = ode45( ...
            @odeFunc,                    ... % handle that evaluates the equations
            linspace( t_0, t_f, N_pts ), ... % the interval of integration
            [ x_t0, v_t0 ],              ... % initial conditions
            odeset( 'RelTol', relError ) ... % relative error tolerance
        );
    
    else
        error([
            'lab1_odeSolver: ''solver'' must be one of: ', ...
            'ode45, ForwardEuler, SemiEuler or RK4' ...
        ]);
    end
    
%=========================================================================================
%% Finaly, calculate the total energy at each point and return the result
    
    E_tot = zeros( size(Y,1), 1 );
    for i = 1:size(Y,1)
        E_tot(i,:) = totalEnergy( Y(i,:) );
    end
    
    result = [ T, Y, E_tot ];
    
%=========================================================================================
%% MyRK4: 4-th Order Runge-Kutta ODE Solver
% For algorithm see: Mathews, John H., Numerical Methods for Mathematics, 
% Science and Engineering, 2nd Ed, Prentice Hall, 1992; 
% Section 9.5, Runge-Kutta Methods (p. 460)
    
    function [ T, Y ] = MyRK4( foo, t_0, t_f, y_0, h, N_pts )

        %%% Arguments
        %   foo:    function handle
        %   t_0:    initial time
        %   t_f:    final time
        %   y_0:    initial value
        %   h:      step-interval
        %   N_pts:  number of points in result

        if h <= 0
            error( 'RK4: ''h'' must be positive real value' );
        end
        
        NT = ( t_f - t_0 ) / h;
        
        if N_pts > NT + 1
            N_pts = NT + 1;
        end

        %%% Initialize state variables
        
        tj = t_0;
        yj = y_0';
        
        T = zeros( NT+1, 1 );
        Y = zeros( NT+1, size( yj, 1 ) );

        T(1) = tj;
        Y(1,:) = yj;

        %%% Loop, solving equation and stepping-forward
        
        for j = 1:NT

            k1 = h * foo( tj,       yj        );
            k2 = h * foo( tj + h/2, yj + k1/2 );
            k3 = h * foo( tj + h/2, yj + k2/2 );
            k4 = h * foo( tj + h,   yj + k3   );

            yj = yj + ( k1 + 2 * k2 + 2 * k3 + k4 )./6;
            tj = t_0 + h * j;
            
            Y(j+1,:) = yj;
            T(j+1) = tj;
        end

        % Return requested number of points
        ix = floor( linspace( 1, NT+1, N_pts ) );
        T = T(ix);
        Y = Y(ix,:);
    end

%=========================================================================================
%% MyEuler: Forward Euler ODE Solver
    
    function [ T, Y ] = MyForwardEuler( foo, t_0, t_f, y_0, h, N_pts )

        %%% Arguments
        %   foo:    function handle
        %   t_0:    initial time
        %   t_f:    final time
        %   y_0:    initial value
        %   h:      step-interval
        %   N_pts:  number of points in result

        if h <= 0
            error( 'Euler: ''h'' must be positive real value' );
        end
        
        NT = ( t_f - t_0 ) / h;
        
        if N_pts > NT + 1
            N_pts = NT + 1;
        end

        %%% Initialize state variables
        
        tj = t_0;
        yj = y_0';
        
        T = zeros( NT+1, 1 );
        Y = zeros( NT+1, size( yj, 1 ) );

        T(1) = tj;
        Y(1,:) = yj;

        %%% Loop, solving equation and stepping-forward
        
        for j = 1:NT

            yj = yj + h * foo( tj, yj );
            tj = t_0 + h * j;
            
            Y(j+1,:) = yj;
            T(j+1) = tj;
        end

        % Return requested number of points
        ix = floor( linspace( 1, NT+1, N_pts ) );
        T = T(ix);
        Y = Y(ix,:);
    end

%=========================================================================================
%% MySemiEuler: Semi-Implicit Euler ODE Solver
% WARNING: Solver requires that:
%   1) the dimension of the vector Y is even 
%   2) Y contains positions in the first half and velocities in the second half

    function [ T, Y ] = MySemiEuler( foo, t_0, t_f, y_0, h, N_pts )

        %%% Arguments
        %   foo:    function handle
        %   t_0:    initial time
        %   t_f:    final time
        %   y_0:    initial value
        %   h:      step-interval
        %   N_pts:  number of points in result

        if h <= 0
            error( 'SemiEuler: ''h'' must be positive real value' );
        end

        NT = ( t_f - t_0 ) / h;
        
        if N_pts > NT + 1
            N_pts = NT + 1;
        end

        %%% Initialize state variables
        
        tj = t_0;
        yj = y_0';
        
        T = zeros( NT+1, 1 );
        Y = zeros( NT+1, size( yj, 1 ) );

        T(1) = tj;
        Y(1,:) = yj;

        % Split indices into X (lower) and V (top) parts
        
        yDim = size( yj, 1 );
        if mod( yDim, 2 ) ~= 0 || yDim < 2
            error( 'SemiEuler: dimension of ''Y'' vector must be even and >= 2' );
        end
        
        xIndices = 1 : yDim/2;        % Bottom part holding positions
        vIndices = yDim/2 + 1 : yDim; % Top part holding velocities

        %%% Loop, solving equation and stepping-forward
        
        for j = 1:NT

            k1 = foo( tj, yj );
            
            % Integrate velocities first and than positions
            yj(vIndices,:) = yj(vIndices,:) + h * k1(vIndices,:);
            yj(xIndices,:) = yj(xIndices,:) + h * yj(vIndices,:);
            tj = t_0 + h * j;
            
            Y(j+1,:) = yj;
            T(j+1) = tj;
        end

        % Return requested number of points
        ix = floor( linspace( 1, NT+1, N_pts ) );
        T = T(ix);
        Y = Y(ix,:);
    end

%=========================================================================================
%% Differential Equation of Motion

    function dY = odeFunc( ~, Y )

        %%% Extract position and velocity from 'Y'
        % Note that the first argument of the function (the time variable) is ignored.

        x = Y(1:2)';
        v = Y(3:4)';

        %%% Calculate the total force

        f_tot = m * g;   % Start with gravitational force

        % Add-up forces of all sorces in the force field

        NS = size( fk, 1 );  % Number of sources

        for k = 1:NS

            % Note: According to profiler, it is slower to use norm( X_p ) when finding 
            % distance from the source.
            X_p = x - pk(k,:);              % displacement from the source
            sq_X_p = X_p(1)^2 + X_p(2)^2;   % squared distance from the source

            f_tot = f_tot            ...
                  - fk(k)            ...
                    / rk(k)^2        ...
                    * exp( -0.5 * sq_X_p / rk(k)^2 ) ...
                    * X_p;
        end

        %%% Differential Equation

        dY = [ v, f_tot / m ]';

    end % function odeFunc

%=========================================================================================
%% totalEnergy - Computes Total Energy of the Particle

    function result = totalEnergy( Y )
        
        %%% Extract position and velocity from 'Y'
        % Note that the first argument of the function (the time variable) is ignored.

        x = Y(1:2);
        v = Y(3:4);

        %%% Calculate total energy as the sum of kinetic and potential energies
        
        E_k = 0.5 * m * ( v * v' );
        
        E_p = - m * ( g * x' );
        for k = 1:size( fk, 1 )
            % Note: According to profiler, it is slower to use norm( X_p ) when finding 
            % distance from the source.
            X_p = x - pk(k,:);              % displacement from the source
            sq_X_p = X_p(1)^2 + X_p(2)^2;   % squared distance from the source
            E_p = E_p - fk(k) * exp( -0.5 * sq_X_p / rk(k)^2 );
        end

        result = E_k + E_p;

    end % function totalEnergy

end % function lab1_odeSolver