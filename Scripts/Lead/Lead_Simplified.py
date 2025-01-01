import numpy as np
import matplotlib.pyplot as plt

## System parameters

T = 298.15          # Temperature (K)

k = (1 / np.log10(np.exp(1)))           # Constant k

C = np.array((96485.3321, 8.31446261815324, 1, 1))          # (0) Faraday constant (C/mol), (1) Universal gas constant (J/K mol), (2) Hydrogen partial pressure (atm), (3) Oxygen partial pressure (atm)

## Gibbs free energy change

def Gibbs_Free_Energy_Proton(pH):
    dG_Proton = -((k * C[1] * T * pH) *1E-03)        # Proton Gibbs free energy (kJ/mol)

    return dG_Proton

def Gibbs_Free_Energy_Electron(E):
    dG_Electron = -((C[0] * E) * 1E-03)          # Electron Gibbs free energy (kJ/mol)

    return dG_Electron

## Potential - pH grid

pH = np.linspace(-2, 16, 5000)         # pH

E = np.linspace(2, -2, 5000)           # redox potential (V)

dG_Proton_Grid, dG_Electron_Grid = np.meshgrid(Gibbs_Free_Energy_Proton(pH), Gibbs_Free_Energy_Electron(E))

## Hydrogen evolution reaction

E_HER = 0.000           # Standard reduction potential HER (V)
n_HER = 2           # Number of electrons exchanged in the HER reaction

E_pH_HER = (E_HER + (((k *C[1] * T) / (n_HER * C[0])) * (- np.log10(C[2]))) + (((-2 * k * C[1] * T) / (n_HER * C[0])) * (pH)))

## Oxygen reduction reaction

E_ORR = 1.229           # Standard reduction potential ORR (V)
n_ORR = 4           # Number of electrons exchanged in the ORR reaction

E_pH_ORR = (E_ORR + (((k *C[1] * T) / (n_ORR * C[0])) * (- np.log10(C[3]))) + (((-4 * k * C[1] * T) / (n_ORR * C[0])) * (pH)))

## Thermodynamic data

M = np.array((1E-06))

dG_H2O = np.array((-237.129))         # Thermodynamic data - (0) H2O

dG_Metal_Species = np.array((0, -24.4, -338.4, -277.6, -226.3, -575.6, -452.2, -187.9, -217.3, -411.8, -601.2))         # Thermodynamic data - (0) Pb, (1) Pb2+, (2) HPbO2-, (3) PbO32-, (4) PbOH+, (5) Pb(OH)3-, (6) Pb(OH)2, (7) PbO, (8) PbO2, (9) Pb2O3, (10) Pb3O4

dG_Metal_Species[1] = (dG_Metal_Species[1] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Pb2+

dG_Metal_Species[2] = (dG_Metal_Species[2] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment HPbO2-

dG_Metal_Species[3] = (dG_Metal_Species[3] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment PbO32-

dG_Metal_Species[4] = (dG_Metal_Species[4] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment PbOH+

dG_Metal_Species[5] = (dG_Metal_Species[5] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Pb(OH)3-

## Gibbs free energy mimimization

dG_H2O_Grid = np.full_like(dG_Proton_Grid, dG_H2O)

dG_Pb_Grid_0 = np.full_like(dG_Proton_Grid, dG_Metal_Species[0])

dG_Pb_Grid_1 = np.full_like(dG_Proton_Grid, dG_Metal_Species[1])

dG_Pb_Grid_2 = np.full_like(dG_Proton_Grid, dG_Metal_Species[2])

dG_Pb_Grid_3 = np.full_like(dG_Proton_Grid, dG_Metal_Species[3])

dG_Pb_Grid_4 = np.full_like(dG_Proton_Grid, dG_Metal_Species[4])

dG_Pb_Grid_5 = np.full_like(dG_Proton_Grid, dG_Metal_Species[5])

dG_Pb_Grid_6 = np.full_like(dG_Proton_Grid, dG_Metal_Species[6])

dG_Pb_Grid_7 = np.full_like(dG_Proton_Grid, dG_Metal_Species[7])

dG_Pb_Grid_8 = np.full_like(dG_Proton_Grid, dG_Metal_Species[8])

dG_Pb_Grid_9 = np.full_like(dG_Proton_Grid, dG_Metal_Species[9])

dG_Pb_Grid_10 = np.full_like(dG_Proton_Grid, dG_Metal_Species[10])

### Reaction 0 - Pb

dG_Reaction_0_Grid = (dG_Pb_Grid_0 - dG_Pb_Grid_0)

### Reaction 1 - Pb2+

dG_Reaction_1_Grid = ((dG_Pb_Grid_1 + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0))

### Reaction 2 - HPbO2-

dG_Reaction_2_Grid = ((dG_Pb_Grid_2 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 3 - PbO32-

dG_Reaction_3_Grid = ((dG_Pb_Grid_3 + (6 * dG_Proton_Grid) + (4 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (3 * dG_H2O_Grid)))

### Reaction 4 - PbOH+

dG_Reaction_4_Grid = ((dG_Pb_Grid_4 + dG_Proton_Grid + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + dG_H2O_Grid))

### Reaction 5 - Pb(OH)3-

dG_Reaction_5_Grid = ((dG_Pb_Grid_5 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (3 * dG_H2O_Grid)))

### Reaction 6 - Pb(OH)2

dG_Reaction_6_Grid = ((dG_Pb_Grid_6 + (2 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 7 - PbO

dG_Reaction_7_Grid = ((dG_Pb_Grid_7 + (2 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + dG_H2O_Grid))

### Reaction 8 - PbO2

dG_Reaction_8_Grid = ((dG_Pb_Grid_8 + (4 * dG_Proton_Grid) + (4 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 9 - Pb2O3

dG_Reaction_9_Grid = (((1/2 * dG_Pb_Grid_9) + (3 * dG_Proton_Grid) + (3 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (3/2 * dG_H2O_Grid)))

### Reaction 10 - Pb3O4

dG_Reaction_10_Grid = (((1/3 * dG_Pb_Grid_10) + (8/3 * dG_Proton_Grid) + (8/3 * dG_Electron_Grid)) -(dG_Pb_Grid_0 + (4/3 * dG_H2O_Grid)))

### Minimum Gibbs free energy

dG_Energy_Minimization = np.stack([dG_Reaction_0_Grid, dG_Reaction_1_Grid, dG_Reaction_2_Grid, dG_Reaction_3_Grid, dG_Reaction_4_Grid, dG_Reaction_5_Grid, dG_Reaction_6_Grid, dG_Reaction_7_Grid, dG_Reaction_8_Grid, dG_Reaction_9_Grid, dG_Reaction_10_Grid], axis = 0)

dG_Pourbaix_Diagram = np.argmin(dG_Energy_Minimization, axis = 0)

dG_pH_Diff = np.diff(dG_Pourbaix_Diagram, axis = 1)

dG_E_Diff = np.diff(dG_Pourbaix_Diagram, axis = 0)

transition_pH = dG_pH_Diff != 0
transition_E = dG_E_Diff != 0

## Pourbaix diagram

plt.rcParams['font.family'] = 'times new roman'         # Set font type
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

plt.rcParams['axes.labelsize'] = 7          # Set font sizes
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['font.size'] = 7

plt.rcParams['xtick.direction'] = 'in'          # Set x-axis properties
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.size'] = 1.5
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['xtick.top'] = True

plt.rcParams['ytick.direction'] = 'in'          # Set y-axis properties
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.size'] = 1.5
plt.rcParams['ytick.minor.width'] = 0.5
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['ytick.right'] = True

plt.rcParams['axes.linewidth'] = 0.5            # Set line/marker properties
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['lines.markersize'] = 0.5

plt.rcParams['legend.frameon'] = False          # Set legend frame

plt.rcParams['savefig.bbox'] = 'tight'          # Set save options
plt.rcParams['savefig.pad_inches'] = 0.01

plt.rcParams['figure.figsize'] = ((100 / 25.4), (1 * 3.5))          # Set the figure size

fig_Pourbaix_Diagram, ax1 = plt.subplots()

def Plot_Pourbaix_Diagram():
    ax1.plot(pH, E_pH_HER, '--', color = '#7E7E7E', alpha = 0.25)
    ax1.plot(pH, E_pH_ORR, '--', color = '#7E7E7E', alpha = 0.25)

    ax1.plot(np.repeat(7, 100), np.linspace(-2, 2, 100), '--', color = '#7E7E7E', alpha = 0.25)

    ax1.plot(np.repeat(6.75, 100), np.linspace(0.75, 1.40, 100), '-', color = '#141414')
    ax1.plot(np.repeat(15, 100), np.linspace(-1.25, -0.50, 100), '-', color = '#141414')

    rows, cols = dG_Pourbaix_Diagram.shape

    for i in range(0, (rows - 1)):
            for ii in range(0, (cols - 1)):
                if (transition_pH[i,ii] == 1) or (transition_E[i,ii] == 1):
                    ax1.plot(pH[ii], E[i], '.', color = '#141414')
                    
    ax1.set_xlabel(r'pH')
    ax1.set_xlim([-2, 16])
    ax1.set_xticks(np.arange(-2, 18, 2))
    ax1.set_ylabel(r'Potential E / V vs. SHE')
    ax1.set_ylim([-2, 2])
    ax1.set_yticks(np.arange(-2, 2.5, 0.5))
    ax1.set_title(r'(b)', loc = 'left')

    plt.text(7.15, 1.85, r'N', fontsize = 6)

    plt.text(0.00, 1.30, r'O$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00, 1.15, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')

    plt.text(0.00, 0.08, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00,-0.08, r'H$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    
    plt.text(4.00, -1.50, r'Pb$_{(s)}$', fontsize = 10)
    plt.text(2.00, 0.50, r'Pb$_{(aq)}^{2+}$', fontsize = 10)
    plt.text(5.00, 1.50, r'PbOH$_{(aq)}^{+}$', fontsize = 10)
    plt.text(8.50, 1.25, r'PbO$_{2 (s)}$', fontsize = 10)
    plt.text(8.75, 0.05, r'Pb(OH)$_{2 (s)}$', fontsize = 10)
    plt.text(11.75, -1.50, r'Pb(OH)$_{3 (aq)}^{-}$', fontsize = 10)
    plt.text(12.75, 1.00, r'PbO$_{3 (aq)}^{2-}$', fontsize = 10)

    fig_Pourbaix_Diagram.tight_layout()
Plot_Pourbaix_Diagram()

def Export_Plot(fig, filename, dpi):
    fig.savefig(f"{filename}.png", dpi = dpi)
    fig.savefig(f"{filename}.jpg", dpi = dpi)
    fig.savefig(f"{filename}.pdf")
    fig.savefig(f"{filename}.svg")

plt.show()
