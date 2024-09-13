# Mars Orbit: Assignment 2
 
This project focuses on the analysis of Mars' orbit using historical opposition data. The aim is to determine Mars' orbit by fitting a model based on opposition observations. This repository includes the Python code used for the analysis, as well as a detailed report of the findings.
 
## Project Description
 
The Mars Orbit assignment involves determining Mars' orbital parameters, such as the position of the equant, angular velocity, and orbit radius, to minimize the angular error between observed and predicted positions. By leveraging Python, we conducted an exhaustive search to optimize these parameters and created visualizations to represent the orbit.
 
### Key Features:
- **Mars Equant Model**: Models Mars' orbit around the Sun, taking into account the position of the equant.
- **Optimization Functions**: Iterative search algorithms to find the best parameters for Mars' orbit.
- **Data Visualization**: Graphical representation of Mars' orbit, with the positions of oppositions relative to the Sun and Aries.
- **Detailed Report**: A comprehensive analysis, including results, parameters, and visual outputs.
 
## Repository Contents
 
- `Assignment2.py`: Contains Python script for calculating Mars orbit parameters, modeling and optimization.
- `mars_data.csv`: The dataset of Mars opposition observations.
- `Report2.pdf`: The final report that summarizes the results of the analysis, including best-fit parameters, visualizations, and implementation details.
- `Plots.png`: Plot1 and Plot2 for visualisation. 
- `README.md`: This file, providing an overview of the project.
- `Mars_Assignment.pdf`: Contains the questions and directions related to the assignment.
 
## Key Results
 
After running the `bestMarsOrbitParams` function, the following best-fit parameters were determined:
 
- **Best Error**: 0.0806 degrees
- **Optimum Center Angle \(c\)**: 149.0 degrees
- **Optimum Equant Distance \(e1\)**: 1.54 units
- **Optimum Equant Angle \(e2\)**: 93.0 degrees
- **Optimum Equant Zero Angle \(z\)**: 55.8 degrees
- **Optimum Orbit Radius \(r\)**: 8.29 units
- **Optimum Angular Speed \(s\)**: 0.524 degrees per day
