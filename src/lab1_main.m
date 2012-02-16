%% Lab1 - Shadowing a Particle, Main Simulation and Animation Engine
% Implements lab1_main function that simulates and visualizes trajectories of larger 
% number of particles and calculates/tracks/plots total energy of the whole system. 
% Compares also simulated trajectories to the reference trajectory.
%
% Filename: lab1_main.m
% Date:     2012-02-14
% Author:   Mikica B Kocic 

function result = lab1_main( testSuiteNo )

%% Function Arguments:
%
%  testSuiteNo -- ID of the test profile corresponding to specific parameter setup.
%                 If omitted, parameters are configured as required for the 'basic level' 
%                 simulation in assignment specification: N_p = 500, forward Euler 
%                 integration with time-step h = 0.01 s and with trajectory animation. 
%                 See 'Test Suite Setup' section bellow for overview of ID values.
%
% Result: Returns matrix 2x6 in the following format:
%
%  [  t_0   x_0  y_0   vx_0  vy_0   E_tot_0  ;
%     t_n   x_n  y_n   vx_n  vy_n   E_tot_n  ]
%
% where the last row contains final position and velocities of the first particle and
% average total energy of the system per particle i.e. E_tot_sys / N_p.

%=========================================================================================
%% Definitions / Simulation Parameters
% Properties that defines the mathematical model and simulation constants,
% e.g., gravity, friction, mass and geometry of the objects etc.

    %%% Simulation and visualization control flags

    flags.Disturbed    = true;  % Disturb position and velocities
    flags.UndisturbN1  = true;  % Do not disturb particle #1
    flags.Animate      = true;  % Render particles in real-time
    flags.SaveFigures  = false; % Save figures (plots) as EPS files
    flags.PlotEnergy   = true;  % Plot kinetic, potential and total energy over time

    %%% Physical constants

    g = [ 0, -9.81 ];           % Acceleration of Gravity, m/s^2

    %%% Integration parameters

    t_0  = 0;                   % Initial time, s
    t_f  = 3;                   % Final (simulation) time, s
    h    = 0.01;                % Time-step, s

    stepper = 'ForwardEuler';   % Forward Euler time-integration
    % stepper = 'SemiEuler';    % Semi-implicit Euler time-integration

    %%% General parameters of the particle system

    N_p  = 500;                 % Number of particles
    m    = 1.0;                 % Particle mass, kg
    L    = 5.0;                 % Characteristic length of the system, m
    x_t0 = [ -L, -L ];          % Particle initial position, m
    v_t0 = [  5, 10 ];          % Particle intital velocity, m/s

    %%% Disturbances of the initial position and velocity

    % Position disturbance, uniformly distributed between bounds
    delX.lower = -0.02 * [L,L]; % lower bound for [x,y]
    delX.upper =  0.02 * [L,L]; % upper bound for [x,y]

    % Velocity disturbance, uniformly distributed between bounds
    delV.lower = -0.02 * v_t0;  % lower bound for [Vx,Vy]
    delV.upper =  0.02 * v_t0;  % upper bound for [Vx,Vy]

    %%% Parameters of the spherically symmetric attractive force fields

    fk = [  32,   40,  28,  16,  20 ]'     ;  % fk(k)   = strength of the force 'k'
    rk = [  0.3, 0.2, 0.4, 0.5, 0.3 ]' * L ;  % rk(k)   = radius of the force 'k'
    pk = [ -0.2,  0.8 ;  ... % pk1            % pk(k,:) = center position of the force 'k'
           -0.3, -0.8 ;  ... % pk2
           -0.6,  0.1 ;  ... % pk3
            0.4,  0.7 ;  ... % pk4
            0.8, -0.3    ... % pk5
         ] * L;

%-----------------------------------------------------------------------------------------
%% Test Suite Setup
% Configure different parameter setups deppending on the test suite id.

    if exist( 'testSuiteNo', 'var' ) == 0  % Not specified test suite id
        testSuiteNo = 0;  % Use default values for parameters and constants
    end
       
    % Overview of test suites:
    %                  default  11   12   21   22   31   32   41   42
    % Save Figures     -        yes  -    -    -    -    -    -    -
    % Animate          yes      yes  yes  yes  yes  -    -    -    -
    % Plot Energy      yes      yes  yes  yes  yes  -    -    -    -
    % Disturbed        yes      yes  yes  -    -    -    -    -    -
    % N_p              500      500  1    1    1    1    500  1    500
    % Stepper          FE       FE   FE   FE   SE   FE   SE   FE   SE
    % 
    switch testSuiteNo
        
        case 11
            flags.SaveFigures = true;
            N_p               = 500;
            
        case 12
            flags.SaveFigures = true;
            N_p               = 1;
            
        case 21
            flags.Animate     = true;
            flags.PlotEnergy  = true;
            flags.Disturbed   = false;
            N_p               = 1;
            stepper           = 'ForwardEuler';
            
        case 22
            flags.Animate     = true;
            flags.PlotEnergy  = true;
            flags.Disturbed   = false;
            N_p               = 1;
            stepper           = 'SemiEuler';
       
        case 31
            flags.Animate     = false;
            flags.PlotEnergy  = false;
            flags.Disturbed   = false;
            N_p               = 1;
            stepper           = 'ForwardEuler';

        case 32
            flags.Animate     = false;
            flags.PlotEnergy  = false;
            flags.Disturbed   = false;
            N_p               = 1;
            stepper           = 'SemiEuler';

        case 41
            flags.Animate     = false;
            flags.PlotEnergy  = false;
            flags.Disturbed   = false;
            N_p               = 500;
            stepper           = 'ForwardEuler';

        case 42
            flags.Animate     = false;
            flags.PlotEnergy  = false;
            flags.Disturbed   = false;
            N_p               = 500;
            stepper           = 'SemiEuler';
        
        case 100 % Have a fun without gravity :)
            flags.PlotEnergy = false;
            N_p  = 50;
            t_f  = 10; % longer simulation
            g    = [ 0, 0 ]; % no gravity, just spherical force fields
            x_t0 = [ 0, 0 ]; % initial position in center
            v_t0 = [ 0, 0 ]; % no initial velocity
            delV.lower = [0,0]; % no velocity disturbance
            delV.upper = [0,0];
            
        case 101 % Generate title image for lab rapport
            flags.PlotEnergy = false;
            N_p  = 20;
            t_f  = 10; % 10 sec simulation
            g    = [ 0, 0 ]; % no gravity, just spherical force fields
            x_t0 = [ 0, 0 ]; % initial position in center
            v_t0 = [ 0, 0 ]; % no initial velocity
    end

%=========================================================================================
%% Initialization
% Create the data structures for the simulation and assigning them initial data.

    %%% Number of time steps

    NT = t_f / h;

    %%% Number of the sources of the force field

    NS = size( fk, 1 );

    %%% Initialize the state of the particles
    % State variables X, V and M are kept as N_p*2 matrices that can be easily 
    % 'vectorized' as vectors with N_p*2 degrees of freedom.

    X = repmat( x_t0, N_p, 1 );
    V = repmat( v_t0, N_p, 1 );
    M = repmat( m,    N_p, 2 );

    if flags.Disturbed % Induce position and veolicty displacement disturbances

        % Calculate position and velocity uniform distribution intervals
        delX.range = delX.upper - delX.lower; % position disturbance range
        delV.range = delV.upper - delV.lower; % velocity disturbance range
        
        % Induce position disturbances in delX inteval
        X = X + bsxfun( @plus, bsxfun( @times, rand( N_p, 2 ), delX.range ), delX.lower );

        % Induce velocity disturbances in delV interval
        V = V + bsxfun( @plus, bsxfun( @times, rand( N_p, 2 ), delV.range ), delV.lower );

        % Nullify disturbances for particle #1 if it's flagged as an 'undisturbed'
        if flags.UndisturbN1
            X(1,:) = x_t0;
            V(1,:) = v_t0;
        end
    end

    X_0 = X; % Take a snapshot of initial X at t0
    V_0 = V; % Take a snapshot of initial V at t0
    
    %%% Initialize computed and derived quantities

    if flags.PlotEnergy
        E_k   = zeros( 1, NT );  % Allocate vector holding kinetic energy snapshots
        E_p   = zeros( 1, NT );  % Allocate vector holding potential energy snapshots
        E_tot = zeros( 1, NT );  % Allocate vector holding total energy snapshots
    end

    %%% Find reference solution using ODE numerical solver

    refSolution = lab1_odeSolver ();

    %%% Extract final position and initial and final total energies.
    %
    % ODE solver returns the result as matrix N_pts * 6 in the following format:
    %  [  t_0   x_0  y_0   vx_0  vy_0   E_tot_0  ;
    %     ...
    %     t_f   x_f  y_f   vx_f  vy_f   E_tot_f  ]

    X_f     = refSolution( end, 2:3 ); % final position
    V_f     = refSolution( end, 4:5 ); % final velocity
    E_tot_0 = refSolution(   1, 6   ); % initial total energy
    E_tot_f = refSolution( end, 6   ); % final total energy

    %%% Render movement of the particles in real-time

    if flags.Animate

        % Get size of the screen

        screenRect    = get( 0, 'ScreenSize' ); % = [ left, bottom, width, height ]
        screen.Width  = screenRect(3);
        screen.Height = screenRect(4);

        % Open new figure centered on screen displaying particles in real-time

        mainFig = figure( ...
            'Name', 'Shadowing a Particle', ... % window title
            'DoubleBuffer', 'on', ... % flicker-free rendering 
            'Position', ... % as [ left, bottom, width, height ]:
            [ screen.Width/4, screen.Height/4, screen.Height/2, screen.Height/2 ], ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperPosition', [ 0, 0, 17, 17 ] ... % left, bottom, width, height
        );

        hold on; % Retain subsequent plots in figure

        % Title, labels and axes:
        title( 'Trajectories',    'FontSize', 11 );
        xlabel( '\it{x} \rm{/m}', 'FontSize', 12, 'FontName', 'Times' );
        ylabel( '\it{y} \rm{/m}', 'FontSize', 12, 'FontName', 'Times' );
        axis( [ -L L -L L ] * 1.1 )  % sets limits for the x- and y-axis
        daspect( [ 1 1 1 ] )   % sets equal scaling on all axes

        % Initialize structure holding trajectories for all particles
        xPath = zeros( [ NT + 1, N_p, 2 ] );
        xPath(1,:,:) = X; % initial positions
        trajPlot = zeros( N_p, 1 ); % trajectory plot handles (fixed bug: N_p instead NT)
        
        % Create trajectory plots for particles 2:N_p (initially invisible)
        for i = 2:N_p
            trajPlot(i) = plot( X(i,1), X(i,2), ':', ...
                'Color', [ 0.87 0.87 0.87 ], ...
                'Visible', 'off' );
        end

        % Draw force field sources as circles with annotated intensities in centra
        for k = 1:size(fk)
            xc = pk(k,1); 
            yc = pk(k,2); 
            rc = rk(k);
            text( xc, yc, [ 'f_', num2str(k), ' = ', num2str( fk(k) ), ' N' ], ...
                'HorizontalAlignment', 'center', ...
                'Color', [ 0.8 0.2 0.2 ] );
            rectangle( ...
                'Position',  [ xc - rc, yc - rc, 2*rc, 2*rc ], ...
                'Curvature', [ 1 1 ], ...
                'LineStyle', ':', ...
                'EdgeColor', [ 0.9 0.4 0.4 ], ...
                'FaceColor', 'none' );
        end

        % Plot all particles as circles but the trajectory only of the first particle
        if flags.UndisturbN1
            % Plot shadow particles as gray circles
            restPlot = plot( X(2:end,1), X(2:end,2), 'o', ...
                'Color', [ 0.65 0.65 0.65 ], 'MarkerSize', 5 );
            % Plot the trajectory of the undisturbed particle in red
            trajPlot(1) = plot( X(1,1), X(1,2), 'r-', 'LineWidth', 2 );
            mainPlot = plot( X(1,1), X(1,2), 'ko', 'MarkerFaceColor', 'r' );
        else
            % Plot all particles as gray circles
            trajPlot(1) = plot( X(1,1), X(1,2), '-', ...
                'Color', [ 0.87 0.87 0.87 ] );
            mainPlot = plot( X(:,1), X(:,2), 'o', ...
                'Color', [ 0.65 0.65 0.65 ], 'MarkerSize', 5 );
        end

        % Plot reference solution in blue smaller square
        refPlot = plot( refSolution(:,2), refSolution(:,3), 'bs', 'MarkerSize', 5 );
        
        % Legend:
        if flags.UndisturbN1 && N_p > 1
            legend( [ refPlot mainPlot restPlot  ], 'Reference', 'Undisturbed', 'Shadow' );
        elseif flags.UndisturbN1
            legend( [ refPlot mainPlot  ], 'Reference', 'Undisturbed'  );
        else
            legend( [ refPlot mainPlot ], 'Reference', 'Shadow'  );
        end
        
    end % flags.Animate

%=========================================================================================
%% Main Simulation Loop

    t = t_0; % Set simulated time
    tic; % Start measuring elapsed (real) time
    tNextRedraw = 0; % Schedule next animation redraw
    
for n = 1:NT
    
    %-------------------------------------------------------------------------------------
    %% Compute Forces and Constraints
    % Compute external forces on the objects and internal forces for interacting objects. 
    % Forces are accumulated to a total force for each object.

    %%% Start with gravitational force
    
    F_tot = bsxfun( @times, M, g );  % (valid also for particles with different masses)

    %%% Add the influence of the force field
    
    for k = 1:NS          % For each source in the force field
        for i = 1:N_p     % For each particle

            % PERFORMANCE NOTES:   (See lab report for profiling)
            % 1) It is ca 10% faster to have inner loop for 'particles' inside 'sources'
            % 2) It is ca 3 times (!) faster *NOT* to use norm( X_p ) when finding
            %    distance from the source
            
            X_p = X(i,:) - pk(k,:);   % displacement from the source
            sq_X_p = X_p(1)^2 + X_p(2)^2;  % squared distance from the source

            F_tot(i,:) = F_tot(i,:)                     ...
                       - fk(k)                          ...
                         / rk(k)^2                      ...
                         * exp( -0.5 * sq_X_p / rk(k)^2 ) ...
                         * X_p;
        end
    end

    %-------------------------------------------------------------------------------------
    %% Step Forward and Update State Variables
    % Advances the simulation data from the current point of time to the next and
    % also computes the new state variables (solving a system of equations).

    if strcmpi( stepper, 'ForwardEuler' ) % Forward Euler time-integration

        X = X + h * V;            % solve position
        V = V + h * F_tot ./ M;   % solve velocity
        t = t_0 + h * n;          % time step (avod t += h that accumulates error)

    elseif strcmpi( stepper, 'SemiEuler' ) % Semi-implicit Euler integration

        V = V + h * F_tot ./ M;   % solve velocity
        X = X + h * V;            % solve position
        t = t_0 + h * n;          % time step (avod t += h that accumulates error)

    else
        error( 'Stepper must be one of: ForwardEuler or SemiEuler!' )
    end
    
    %-------------------------------------------------------------------------------------
    %% Update Derived Quantities
    % After the state variables are updated, we can update all derived quantities.

    %%% Compute kinetic, potential and total energy
    
    E_k_n = 0.5 * ( ( M(:) .* V(:) )' * V(:) );   % Kinetic energy

    % Potential energy is the sum of gravitational potential energy
    % and the potential energy of the force field (calculated in the loop)

    E_p_n = - trace( X' * bsxfun( @times, M, g ) );

    for k = 1:NS          % For each source in the force field
        for i = 1:N_p     % For each particle
            % PERFORMANCE NOTES:
            % 1) It is ca 10% faster to have inner loop for 'particles' inside 'sources'
            % 2) It is ca 3 times (!) faster *NOT* to use norm( X_p ) when finding
            %    distance from the source
            X_p = X(i,:) - pk(k,:);         % displacement from the source
            sq_X_p = X_p(1)^2 + X_p(2)^2;   % squared distance from the source
            E_p_n = E_p_n - fk(k) * exp( -0.5 * sq_X_p / rk(k)^2 );
        end
    end

    E_tot_n = E_k_n + E_p_n;  % total energy

    % Store energy levels for later display
    
    if flags.PlotEnergy
        E_k(n)   = E_k_n;
        E_p(n)   = E_p_n;
        E_tot(n) = E_tot_n;
    end
    
    %-------------------------------------------------------------------------------------
    %% Simulation I/O (handle input and output of the simulation)
    % Input can be user interaction or data streaming from another simulation or hardware 
    % in the loop. Output can be data storage for post-processing, realtime graphics 
    % rendering, signals to haptic force-feedback unit.

    %%% Real-Time Graphics Rendering
    
    if flags.Animate
        
        % Remember current positions of all particles at each integration step
        xPath(n+1,:,:) = X;
    
        % Update the X/Y-data of the real-time plot
        if flags.UndisturbN1
            set( mainPlot,    'XData', X(1,1),         'YData', X(1,2)         );
            set( restPlot,    'XData', X(2:end,1),     'YData', X(2:end,2)     );
            set( trajPlot(1), 'XData', xPath(1:n,1,1), 'YData', xPath(1:n,1,2) );
        else
            set( mainPlot, 'XData', X(:,1), 'YData', X(:,2) );
        end

        % Pause if rendering is too fast for the real-time animation.
        % The following code skips frames if loop steps are lengthy.
        
        toSleep = ( t - t_0 ) - toc; % sleep time = simulation time - real time
        
        if toSleep > 0 
            % Pause, if required to slow-down the simulation (to be in real-time)
            pause( toSleep );
            tNextRedraw = toc + 0.040;
        elseif toc >= tNextRedraw
            % Update every 40 ms (25 Hz) in lengthy calculations
            drawnow;
            tNextRedraw = toc + 0.040;
        end
    end

end % of the Main Simulation Loop --------------------------------------------------------

%=========================================================================================
%% Post-processing
% Processing of data stored during the simulation into quantities as required
% e.g., energy, temperature, velocity fields, time-averaged force. Data is stored in
% file or made into graphs, tables or animation. The system state variables may be saved 
% and used for initialization data for another simulation.

    fprintf( '\n---------------------- lab1_main\n\n' );
    fprintf( 'Number of particles .: %d\n',    N_p             );
    fprintf( 'Disturbed x0/v0? ....: %d\n',    flags.Disturbed );
    fprintf( 'Integration method ..: %s\n',    stepper         );
    fprintf( 'Time step ...........: %g s\n',  h               );
    fprintf( 'Number of steps .....: %d\n',    NT              );
    fprintf( 'Simulation time .....: %g s\n',  t               );
    fprintf( 'Elapsed time ........: %g s\n',  toc             );
    fprintf( '\nResults:\n\n' );
    
    %%% Make visible trajectories of all particles and print the figure

    if flags.Animate
        
        for i = 1 : N_p
            set( trajPlot(i), ...
                'XData', xPath(:,i,1), 'YData', xPath(:,i,2), ...
                'Visible', 'on' );
        end
        
        drawnow; % Flush event queue and update figure window
        
        % Save the figure as EPS file
        if flags.SaveFigures
            print( mainFig, '-r300', '-depsc2', ...
                [ 'lab1_fig1a_', datestr( now, 'YYYYmmddHHMMSS' ), '.eps' ] );
        end
    end

    %%% Plot kinetic, potential and total energyover time

    if flags.PlotEnergy

        fig = figure( 'name', 'Energies', ... %Create a new figure window
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperPosition', [ 0, 0, 17, 10 ] ... % left, bottom, width, height
        );
    
        hold on;           % Retain subsequent plots in figure
        T = h * (1:NT);    % Create a time-vector for plotting of energies

        % Title, axes and grid:
        title( 'Kinetic, Potential and Total Energy', 'FontSize', 11 );
        xlabel( '\it{t} \rm{/s}', 'FontSize', 12, 'FontName', 'Times' );
        ylabel( '\it{E} \rm{/J}', 'FontSize', 12, 'FontName', 'Times' );
        grid on;

        % Plots:
        plot( T, E_k,   'b-',  'LineWidth', 0.5 ); % kinetic energy: blue, thin
        plot( T, E_p,   'r-.', 'LineWidth', 0.5 ); % potential energy: red, dot-sash, thin
        plot( T, E_tot, 'k-',  'LineWidth', 1.5 ); % total energy: black, thick

        % Legend:
        hleg = legend( '\it{E}\rm{_k}', '\it{E}\rm{_p}', '\it{E}\rm{_{tot}}' );
        set( hleg, 'FontSize', 11, 'FontName', 'Times' );
        
        drawnow; % Flush event queue and update figure window

        % Save the figure as EPS file
        if flags.SaveFigures
            print( fig, '-painters', '-depsc2', ...
                [ 'lab1_fig1b_', datestr( now, 'YYYYmmddHHMMSS' ), '.eps' ] );
        end
    end

    %%% Show the displacement of the first particle at time t_f

    fprintf( 'X_f:     %14.9f, %0.9f\n', X_f(1), X_f(2) );
    fprintf( 'X_n:     %14.9f, %0.9f\n', X(1,1), X(1,2) );
    fprintf( 'Error:   %14.9f\n\n', norm( X(1,:) - X_f ) );

    %%% Compute the minimum, average, standard deviation and maximum position and
    %%% velocity displacements of particles at times t_0 and t_f

    if N_p > 1
        %%% INITIAL POSITION STATISTICS -------------------------------------------
        % Calculate position displacement array X_p(i) = || X(i)_0 - X_0 ||
        X_p = bsxfun( @minus, X_0, repmat( x_t0, N_p, 1 ) ).^2;
        X_p = ( X_p(:,1) + X_p(:,2) ).^0.5; % calculate norm of individual X_p(i,:)
        % Calculate velocity displacement array V_p(i) = || V(i)_0 - V_0 ||
        V_p = bsxfun( @minus, V_0, repmat( v_t0, N_p, 1 ) ).^2;
        V_p = ( V_p(:,1) + V_p(:,2) ).^0.5; % calculate norm of individual V_p(i,:)
        
        % Ignore undisturbed (first) particle, if any:
        if flags.UndisturbN1; X_p(1,:) = []; V_p(1,:) = []; end 

        % Display minimum, maximum, average and standard deviation of X_p and V_p
        fprintf( 'Statistics for initial disturbances:\n' );
        fprintf( '--------\t%10s\t%10s\n', 'del_X', 'del_V' );
        fprintf( 'Minimum:\t%0.9f\t%0.9f\n',   min  ( X_p ), min  ( V_p ) );
        fprintf( 'Maximum:\t%0.9f\t%0.9f\n',   max  ( X_p ), max  ( V_p ) );
        fprintf( 'Average:\t%0.9f\t%0.9f\n',   mean ( X_p ), mean ( V_p ) );
        fprintf( 'Std Dev:\t%0.9f\t%0.9f\n\n', std  ( X_p ), std  ( V_p ) );
        
        %%% FINAL POSITION STATISTICS ---------------------------------------------
        % Calculate position displacement array X_p(i) = || X(i)_n - X_f ||
        X_p = bsxfun( @minus, X, X_f ).^2; % from final positions
        X_p = ( X_p(:,1) + X_p(:,2) ).^0.5; % calculate norm of individual X_p(i,:)
        % Calculate velocity displacement array V_p(i) = || V(i)_n - V_f ||
        V_p = bsxfun( @minus, V, V_f ).^2; % from final positions
        V_p = ( V_p(:,1) + V_p(:,2) ).^0.5; % calculate norm of individual V_p(i,:)

        % Ignore undisturbed (first) particle, if any:
        if flags.UndisturbN1; X_p(1,:) = []; V_p(1,:) = []; end 
        
        % Display minimum, maximum, average and standard deviation of X_p and V_p
        fprintf( 'Statistics for final disturbances:\n' );
        fprintf( '--------\t%10s\t%10s\n', 'del_X', 'del_V' );
        fprintf( 'Minimum:\t%0.9f\t%0.9f\n',   min  ( X_p ), min  ( V_p ) );
        fprintf( 'Maximum:\t%0.9f\t%0.9f\n',   max  ( X_p ), max  ( V_p ) );
        fprintf( 'Average:\t%0.9f\t%0.9f\n',   mean ( X_p ), mean ( V_p ) );
        fprintf( 'Std Dev:\t%0.9f\t%0.9f\n\n', std  ( X_p ), std  ( V_p ) );
    end

    %%% Calculate and display total energy absolute error

    fprintf( 'Reference trajectory energy:\n' );
    fprintf( 'E_tot_0: %14.9f\n', E_tot_0 );
    fprintf( 'E_tot_f: %14.9f\n', E_tot_f );
    fprintf( 'Abs Err: %14.9f\n', abs( E_tot_f - E_tot_0) );
    fprintf( '\n' );
    fprintf( 'Simulated trajectory energy (average per particle):\n' );
    fprintf( 'E_tot_0: %14.9f\n', E_tot_0  );
    fprintf( 'E_tot_n: %14.9f\n', abs( E_tot_n / N_p ) );
    fprintf( 'Abs Err: %14.9f\n', abs( E_tot_n / N_p - E_tot_0 ) );
    fprintf( '\n' );

    %%% Return the result
    
    result = [ t_0, x_t0,   v_t0,   E_tot_0       ; ...
               t,   X(1,:), V(1,:), E_tot_n / N_p ];
            
end % function lab1_main
