# Finite Element Method (FEM) for Axial Displacement in Bars

This MATLAB program implements the **Finite Element Method (FEM)** to analyze the **natural mode frequencies** and **mode shapes** of a segmented bar with varying dimensions. The program computes the stiffness and inertia matrices of each section of the bar and assembles them into global matrices to determine the system's eigenvalues and eigenvectors.  

## Features
- User input for interpolation type (linear or quadratic) and the number of elements per section.
- Analysis of bars with segmented properties such as length, cross-sectional area, elasticity, and density.
- Computation of:
  - **Stiffness Matrices**
  - **Mass/Inertia Matrices**
  - **Modes (Eigenvectors)**
  - **Natural Frequencies (Eigenvalues)**
- Application of boundary conditions to solve for reduced matrices.

---

## Table of Contents
1. [Requirements](#requirements)
2. [How to Run](#how-to-run)
3. [Input Details](#input-details)
4. [Program Workflow](#program-workflow)
5. [Key MATLAB Functions](#key-matlab-functions)
6. [Example Use Case](#example-use-case)
7. [Output](#output)
8. [License](#license)

---

## Requirements
- MATLAB R2021a or later.
- Basic understanding of Finite Element Analysis (FEA) and Vibrations.

---

## How to Run
1. Download the script `fem_axial_displacement.m`.
2. Open MATLAB and navigate to the folder containing the script.
3. Run the script:
   ```matlab
   fem_axial_displacement
   ```
4. Follow the prompts to input the required data.

---

## Input Details
The program requires the following user inputs:
1. **Interpolation Type**:
   - `1` for Linear Interpolation.
   - `2` for Quadratic Interpolation.
2. **Number of Elements per Section**:
   - A number between `1` and `100`.

The bar's dimensions, including lengths, cross-sectional areas, elasticity, and density, are pre-defined within the program but can be modified in the `defineSectionProperties` function.

---

## Program Workflow
1. **Define Section Properties**:
   - Define the segmented bar's length, cross-sectional area, elasticity, and density.

2. **Input Parameters**:
   - User selects interpolation type and number of elements.

3. **Compute Element Matrices**:
   - Local stiffness and inertia matrices are calculated for each element.

4. **Assemble Global Matrices**:
   - Global stiffness and inertia matrices are assembled using local matrices.

5. **Apply Boundary Conditions**:
   - Boundary conditions are applied to reduce matrices.

6. **Solve Eigenvalue Problem**:
   - The program solves for eigenvalues (natural frequencies) and eigenvectors (mode shapes).

---

## Key MATLAB Functions
The program uses the following key functions:
- `defineSectionProperties`: Defines the bar's segmented properties.
- `getValidInput`: Ensures valid user input.
- `computeElementMatrices`: Computes stiffness and inertia matrices for each element.
- `assembleGlobalMatrices`: Assembles global matrices from element matrices.
- `applyBoundaryConditions`: Applies boundary conditions to global matrices.
- `solveEigenProblem`: Solves for eigenvalues and eigenvectors.

---

## Example Use Case
1. A segmented bar with varying cross-sectional areas and lengths:
   - **Section Lengths**: [4 in, 6 in, 2 in, 3 in]
   - **Cross-Sectional Area**:
     - Section 1: Parabolic variation.
     - Section 2: Constant.
     - Section 3: Linear variation.
     - Section 4: Constant.
   - **Elasticity**: `10,152,642 psi`
   - **Density**: `0.0975 lbm/in³`

2. User selects:
   - Linear interpolation (`1`).
   - Number of elements: `4`.

3. Output includes:
   - **Mode Shapes (Eigenvectors)**.
   - **Natural Frequencies (Eigenvalues)**.

---

## Output
The program displays the following:
1. **Eigenvectors**: Mode shapes of the bar.
2. **Eigenvalues**: Natural frequencies of the bar.

These results can be used to analyze the bar's vibration characteristics under the given dimensions and properties.

---

## License
This project is for educational purposes. Feel free to use or modify it for non-commercial applications.

---

