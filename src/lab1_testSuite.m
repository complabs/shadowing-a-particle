%% Lab1 - Shadowing a Particle, Test Suite
% Runs all simulation required for basic- and advanced level taks.
%
% Filename: lab1_testSuite.m
% Date:     2012-02-14
% Author:   Mikica B Kocic 

%=========================================================================================
%% Restart Simulation 

    clear all;   % Remove all functions, variables and global variables from workspace
    close all;   % Delete all figures whose handles are not hidden
    clc;         % Clear command window

%=========================================================================================
%% Definitions / Simulation Parameters

    flags.SaveFigures  = false;   % Save figures (plots) as EPS files

%=========================================================================================
%% Main Simulation
% - Calculate reference trajectory (using ODE solver, for single particle)
% - Disturb initial conditions
% - Integrate equations (for large N_p, in vector form)
% - Visualize particle positions at each time-step and the reference trajectory
% - Compute the minimum, average and maximum displacement of the particles' positions
% - Calculate and plot energy level (kinetic, potential and total)

    fprintf( '====================== Main Simulation: Mandatory & Extra Task 1 ====\n' );

    if flags.SaveFigures
        lab1_main( 11 );  % Configuration suite #11 saves figures
    else
        lab1_main;
    end

%=========================================================================================
%% Advanced

    fprintf( '====================== Advanced Levels: Extra Tasks 2-4 =============\n' );
    
    % Start progress bar
    wb = waitbar( 0, 'Running advanced tasks...', 'Name', 'Lab1: Advanced tasks' );
    
    % Get the reference solution using Runge-Kutta 4 with time-step 1e-4
    
    fprintf( '\nIntegrating reference solution; this make take a while... ' );
    drawnow; % Flush event queue and update figure window
    
    waitbar( 0.1, wb, 'Finding reference solution...' );
    
    refSol  = lab1_odeSolver( 'RK4', 2, 0, 1e-4 );
    X_f     = refSol( 2, 2:3 );
    E_tot_f = refSol( 2, 6 );
    
    % Array of integration time-steps
    H = [ 0.5, 0.1, 0.05, 0.01, 0.005, 0.001 ];
    
    % Integration methods
    algos = { 'ForwardEuler', 'SemiEuler', 'RK4' };
    
    % Collect results in matrix 'R' where the first colum contains 'h' value and
    % the rest columns simulation results error X and error E_tot (2 columns per 
    % integration algorithm)
    
    R = zeros( length(H), 1 + 2 * length(algos) );  
    R(:,1) = H;   % put h values into the first column
    
    fprintf( 'Done.\n' );
    
    % Calculate trajectory using different integration algorithms
    
    fprintf( 'Running simulations using different integration methods... ' );
    
    for j = 1 : length(algos)
        for i = 1 : length(H)

            waitbar( 0.1 + 0.9 * ( (j-1) * length(H) + i )/length(algos)/length(H), ...
                wb, [ 'Running: ', algos{j}, ', h = ', num2str(H(i)) ] );
    
            % Solve lab1 ODE using given algorithm with given time-step
            
            sol = lab1_odeSolver( algos{j}, 2, 0, H(i) );

            % Collect results
            
            X_n = sol(2,2:3);   % get final X(n)
            E_tot_0 = sol(1,6); % get initial E_tot(0)
            E_tot_n = sol(2,6); % get final E_tot(n)

            R( i, j*2     ) = norm( X_n - X_f );        % err X
            R( i, j*2 + 1 ) = abs( E_tot_n - E_tot_0 ); % err E_tot
        end
    end

    fprintf( 'Done.\n' );
    
    % Show & sync final status in progress bar and then close it
    waitbar( 1, wb, 'Done!' ); drawnow; pause( 0.1 ); close( wb );
    
    %%% Display results ------------------------------------------------------------------
    
    % Display displacement position errors (tabular)
    
    fprintf( '\nPosition Absolute Errors\n\n' );
    
    fprintf( 'h\t%s\t%s\t%s\n', algos{1}, algos{2}, algos{3} );
    for i = 1:size(R,1)
        fprintf( '%g', R(i,1) );
        for j = 2:2:size(R,2)
            fprintf( '\t%0.12f', R(i,j) )
        end
        fprintf( '\n' );
    end
    
    % Display energy errors (tabular)
    
    fprintf( '\nEnergy Absolute Errors\n\n' );
    
    fprintf( 'h\t%s\t%s\t%s\n', algos{1}, algos{2}, algos{3} );
    for i = 1:size(R,1)
        fprintf( '%g', R(i,1) );
        for j = 3:2:size(R,2)
            fprintf( '\t%0.12f', R(i,j) )
        end
        fprintf( '\n' );
    end
    
    %%% Plot results ---------------------------------------------------------------------

    % Plot position displacement errors ---------
    
    fig = figure( 'Name', 'Displacement Abs Errors', ...
        'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
        'PaperPosition', [ 0, 0, 17, 10 ] ... % left, bottom, width, height
    );
    
    linedefs = { 'ko-', 'r^-', 'bs-' };
    for j = 1 : length(algos)
        loglog( R(:,1), R(:,j*2), linedefs{j} );
        hold on;
    end
    
    % Title, legend, labels and grid:
    grid;
    title( 'Displacement Absolute Errors' );
    legend( algos, 'Location', 'Best' );
    xlabel( '\it{h} \rm{/s}', 'FontSize', 11, 'FontName', 'Times' );
    ylabel( '\Delta\it{X} \rm{/m}', 'FontSize', 11, 'FontName', 'Times' );

    % Save figure as EPS file
    if flags.SaveFigures
        print( fig, '-depsc2', ...
            [ 'lab1_fig6a_', datestr( now, 'YYYYmmddHHMMSS' ), '.eps' ] );
    end
    
    %%% Plot energy errors ---------
    
    fig = figure( 'Name', 'Energy Abs Errors', ...
        'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
        'PaperPosition', [ 0, 0, 17, 10 ] ... % left, bottom, width, height
    );
    
    for j = 1 : length(algos)
        loglog( R(:,1), R(:,j*2+1), linedefs{j} );
        hold on;
    end
    
    % Title, legend, labels and grid:
    grid;
    title( 'Energy Absolute Errors' );
    legend( algos, 'Location', 'Best' );
    xlabel( '\it{h} \rm{/s}', 'FontSize', 11, 'FontName', 'Times' );
    ylabel( '\Delta\it{E}_{\rm{tot}} \rm{/J}', 'FontSize', 11, 'FontName', 'Times' );

    % Save figure as EPS file
    if flags.SaveFigures
        print( fig, '-depsc2', ...
            [ 'lab1_fig6b_', datestr( now, 'YYYYmmddHHMMSS' ), '.eps' ] );
    end

    fprintf( '\n\n' );

    drawnow; % Flush event queue and update figure window
    
    % Cascade existing figures so that they don't directly overlap
    
    figs = flipud( findobj( 0, 'Type', 'figure' ) );
    for n = 2:length(figs)
        pPrev = get( figs(n-1), 'Position' );
        pos = get( figs(n), 'Position' );
        pos(1:2) = pPrev(1:2) + [ 30, pPrev(4) - 30 ] - [ 0, pos(4) ];
        set( figs(n), 'Position',  pos );
    end
    
    fprintf( '====================== COMPLETED ====================================\n' );

%=========================================================================================