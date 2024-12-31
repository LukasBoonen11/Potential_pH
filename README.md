# Potential_pH
Pourbaix Diagram Construction - Python Scripts

# Overview

The Python scripts contained in this set have been developed for the construction and analysis of Pourbaix diagrams in electrochemical systems. The scripts' primary objective is to visualize the stability regions of different species in the system, a task that is achieved through thermodynamic calculations that involve Gibbs free energy changes. The utilization of these diagrams facilitates the assessment of the potential-pH relationships of metal species as well as different solid and dissolved metal species, and their respective stability, particularly in the presence of various environmental anions, such as chloride, sulfate, and carbonate.

# Requirements

To run the scripts, you need the following Python libraries:

NumPy: for numerical calculations and array handling.

Matplotlib: for plotting the diagrams.

# Functionality

The Python scripts function as follows:
  
  Thermodynamic calculations:
  
  The Gibbs free energy change for each reaction is calculated using the standard Gibbs free energy of formation (ΔfG) of the reactants and products.
  The overall Gibbs free energy change takes into account proton (H⁺) and electron (e⁻) exchange reactions.
  The thermodynamic stability of each species is determined by calculating the Gibbs free energy change for all possible reactions and identifying the most thermodynamically favorable path.

  Diagram construction:
  
  The diagrams are constructed using Matplotlib, with pH values on the x-axis and electrochemical potential (in volts versus SHE) on the y-axis.
  Stability regions are depicted for various species as regions enclosed by boundary lines.
  The stability of water is highlighted with sloping dashed lines representing the boundaries for hydrogen evolution and oxygen reduction.

  Customization:
  
  The user can modify the input concentrations, and potential/pH ranges as required.
  Additional species and reactions can be added to extend the system of interest.

# How to Use

  To use the scripts for generating Pourbaix diagrams:

  Run the Script: Execute the Python script via the command line.
  View the Output: The script will output a Pourbaix diagram as a visual plot, which can be saved to a file (e.g., PNG, PDF) or displayed interactively.
  Analyze the Results: The diagram will show the stability regions of different species in the system.

# Conclusion
  
  These scripts offer a powerful tool for visualizing the thermodynamic stability of metal species and their corrosion behavior across various pH and potential conditions.
  By incorporating environmental factors like chloride, sulfate, and bicarbonate ions, the diagrams provide a more realistic representation of the system under investigation, making them highly applicable for studies on corrosion and   electrochemical reactions in diverse environments.
