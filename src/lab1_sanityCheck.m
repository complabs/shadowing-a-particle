%% Lab1 - Shadowing a Particle, Sanity-Check Tests
% Tests sanity of numerical integration algorithms in lab1_odeSolver and lab1_main.
%
% Filename: lab1_sanityCheck.m
% Date:     2012-02-14
% Author:   Mikica B Kocic 

%=========================================================================================
%% Restart Simulation 

    clear all;   % Remove all functions, variables and global variables from workspace
    close all;   % Delete all figures whose handles are not hidden
    clc;         % Clear command window

%=========================================================================================
%% sanity check of the numerical integration algorithms lab1_odeSolver vs lab1_main 

    hRef   = 1e-4;
    h      = 1e-2;
    relTol = 1e-7;
    
    fprintf( '====================== SANITY CHECKS ================================\n' );
    
    format long;
    
    fprintf( '\n---------------------- Reference Solution: RK4, h = %g\n\n', hRef );

    refSol  = lab1_odeSolver( 'RK4', 2, 0, hRef );
    disp( refSol' );

    fprintf( '\n---------------------- odeSolver: ode45, RelTol = %g\n\n', relTol );
    
    sol = lab1_odeSolver( 'ode45', 3, relTol, 0 );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- odeSolver: RK4, h = %g \n\n', h );
    
    sol = lab1_odeSolver( 'RK4', 2, 0, h );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- odeSolver: ForwardEuler, h = %g\n\n', h );
    
    sol = lab1_odeSolver( 'ForwardEuler', 2, 0, h );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- main: ForwardEuler, h = %g, N_p = 1\n\n', h );
    
    [ ~, sol ] = evalc( 'lab1_main( 31 )' );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- main: ForwardEuler, h = %g, N_p = 500\n\n', h );
    
    [ ~, sol ] = evalc( 'lab1_main( 41 )' );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- odeSolver: SemiEuler, h = %g\n\n', h );
    
    sol = lab1_odeSolver( 'SemiEuler', 2, 0, h );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- main: SemiEuler, h = %g, N_p = 1\n\n', h );
    
    [ ~, sol ] = evalc( 'lab1_main( 32 )' );
    disp( sol(end,:)' - refSol(end,:)' );
    
    fprintf( '\n---------------------- main: SemiEuler, h = %g, N_p = 500\n\n', h );
    
    [ ~, sol ] = evalc( 'lab1_main( 42 )' );
    disp( sol(end,:)' - refSol(end,:)' );
   
    format;
    
    fprintf( '====================== COMPLETED ====================================\n' );
    
%=========================================================================================