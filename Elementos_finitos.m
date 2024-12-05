% Main Program for Finite Element Method (FEM) for Axial Displacement in Bars

clear; clc; close all;

% Step 1: Define Section Properties
sections = defineSectionProperties();

% Step 2: User Input for Interpolation and Elements
disp('Finite Element Method for Bars');
fprintf('\n');
n = getValidInput('Select interpolation type (1 for linear, 2 for quadratic): ', [1, 2]);
k = getValidInput('Enter the number of elements per section: ', 1:100);

% Step 3: Compute Element and Global Matrices
[K_elements, M_elements] = computeElementMatrices(sections, n, k);
[K_global, M_global] = assembleGlobalMatrices(K_elements, M_elements, n, k);

% Step 4: Apply Boundary Conditions
boundary_nodes = 1; % Example: First node fixed
K_reduced = applyBoundaryConditions(K_global, boundary_nodes);
M_reduced = applyBoundaryConditions(M_global, boundary_nodes);

% Step 5: Solve for Modes and Frequencies
disp('Modes and Natural Frequencies');
[eigenVectors, eigenValues] = solveEigenProblem(K_reduced, M_reduced);
disp('Eigenvectors (Mode Shapes):');
disp(eigenVectors);
disp('Eigenvalues (Natural Frequencies):');
disp(eigenValues);

%% Supporting Functions (Embedded Below)

% Function: Define Section Properties
function sections = defineSectionProperties()
    sections = struct();
    sections.lengths = [4/12, 6/12, 2/12, 3/12]; % Lengths of each section
    sections.areas = {@(x) pi*(-2*x+14)^2, @(x) pi*(20/2)^2, ...
                      @(x) pi*(-x+20)^2, @(x) pi*8^2}; % Area functions
    sections.elasticity = 10152642; % Elastic modulus (psi)
    sections.density = 0.0975437; % Density (lbm/in^3)
end

% Function: Get Valid Input
function value = getValidInput(prompt, valid_range)
    while true
        value = input(prompt);
        if ismember(value, valid_range)
            break;
        else
            disp('Invalid input. Try again.');
        end
    end
end

% Function: Compute Element Matrices
function [K_elements, M_elements] = computeElementMatrices(sections, n, k)
    syms x
    lengths = sections.lengths;
    areas = sections.areas;
    E = sections.elasticity;
    rho = sections.density;

    % Define interpolation functions
    if n == 1
        N = [((x - 1)/(2 - 1)), ((x - 2)/(1 - 2))];
    elseif n == 2
        x_mid = 1.5; % Midpoint for quadratic interpolation
        N = [((x - 2)*(x - x_mid))/((1 - 2)*(1 - x_mid)), ...
             ((x - 1)*(x - x_mid))/((2 - 1)*(2 - x_mid)), ...
             ((x - 1)*(x - 2))/((x_mid - 1)*(x_mid - 2))];
    end

    K_elements = cell(length(lengths), 1);
    M_elements = cell(length(lengths), 1);

    for s = 1:length(lengths)
        % Divide section into elements
        elem_length = lengths(s) / k;
        lim_elem = (0:k) * elem_length;

        % Initialize stiffness and inertia matrices
        K_section = cell(k, 1);
        M_section = cell(k, 1);

        % Compute matrices for each element
        for e = 1:k
            x1 = lim_elem(e);
            x2 = lim_elem(e + 1);

            % Stiffness Matrix
            A = areas{s};
            K_elem = zeros(n + 1);
            for i = 1:(n + 1)
                for j = 1:(n + 1)
                    integrand = A(x) * diff(N(i), x) * diff(N(j), x);
                    K_elem(i, j) = E * double(int(integrand, x, x1, x2));
                end
            end
            K_section{e} = K_elem;

            % Inertia Matrix
            M_elem = zeros(n + 1);
            for i = 1:(n + 1)
                for j = 1:(n + 1)
                    integrand = A(x) * N(i) * N(j);
                    M_elem(i, j) = rho * double(int(integrand, x, x1, x2));
                end
            end
            M_section{e} = M_elem;
        end
        K_elements{s} = K_section;
        M_elements{s} = M_section;
    end
end

% Function: Assemble Global Matrices
function [K_global, M_global] = assembleGlobalMatrices(K_elements, M_elements, n, k)
    num_nodes = n * k + 1; % Total number of nodes
    K_global = zeros(num_nodes);
    M_global = zeros(num_nodes);

    % Loop through each section and assemble matrices
    for s = 1:length(K_elements)
        K_section = K_elements{s};
        M_section = M_elements{s};

        for e = 1:length(K_section)
            for i = 1:(n + 1)
                for j = 1:(n + 1)
                    global_i = i + (e - 1);
                    global_j = j + (e - 1);

                    % Add to global stiffness and inertia matrices
                    K_global(global_i, global_j) = K_global(global_i, global_j) + K_section{e}(i, j);
                    M_global(global_i, global_j) = M_global(global_i, global_j) + M_section{e}(i, j);
                end
            end
        end
    end
end

% Function: Apply Boundary Conditions
function K_reduced = applyBoundaryConditions(K_global, boundary_nodes)
    K_reduced = K_global;
    K_reduced(boundary_nodes, :) = [];
    K_reduced(:, boundary_nodes) = [];
end

% Function: Solve Eigenvalue Problem
function [eigenVectors, eigenValues] = solveEigenProblem(K_reduced, M_reduced)
    [eigenVectors, eigenValues] = eig(K_reduced, M_reduced);
end

