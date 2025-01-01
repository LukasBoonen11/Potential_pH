import numpy as np
import matplotlib.pyplot as plt

## System parameters

T = 298.15          # Temperature (K)

k = (1 / np.log10(np.exp(1)))           # Constant k

C = np.array((96485.3321, 8.31446261815324, 1, 1))          # (0) Faraday constant (C/mol), (1) Universal gas constant (J/K mol), (2) Hydrogen partial pressure (bar), (3) Oxygen partial pressure (bar)

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


N = np.array((1E-09))

Cl = np.array((2E-03))


SO42 = np.array((0.75E-03))


CO32 = np.array((1.15E-03))


dG_H2O = np.array((-237.1))         # Thermodynamic data - (0) H2O

dG_Metal_Species = np.array((0, -147.1, -330.1, -384.2, -457.1, -694.2, -858.5, -553.8, -318.3))         # Thermodynamic data - (0) Zn, (1) Zn2+, (2) ZnOH+, (3) ZnO22-, (4) HZnO2-, (5) Zn(OH)3-, (6) Zn(OH)4-, (7) Zn(OH)2, (8) ZnO

dG_Metal_Species[1] = (dG_Metal_Species[1] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Zn2+

dG_Metal_Species[2] = (dG_Metal_Species[2] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment ZnOH+

dG_Metal_Species[3] = (dG_Metal_Species[3] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment ZnO22-

dG_Metal_Species[4] = (dG_Metal_Species[4] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment HZnO2-

dG_Metal_Species[5] = (dG_Metal_Species[5] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Zn(OH)3-

dG_Metal_Species[6] = (dG_Metal_Species[6] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Zn(OH)4-


dG_Chloride_Species = np.array((-131.2, -275.3, -540.5, -666.0, -369.4, -472.7))          # Thermodynamic data - (0) Cl-, (1) ZnCl+, (2) ZnCl3-, (3) ZnCl42-, (4) ZnCl2, (5) Zn(OH)Cl

dG_Chloride_Species[0] = (dG_Chloride_Species[0] + (C[1] * T * np.log(Cl) * 1E-03))          # Concentration (activity) adjustment Cl-

dG_Chloride_Species[1] = (dG_Chloride_Species[1] + (C[1] * T * np.log(N) * 1E-03))          # Concentration (activity) adjustment ZnCl+

dG_Chloride_Species[2] = (dG_Chloride_Species[2] + (C[1] * T * np.log(N) * 1E-03))          # Concentration (activity) adjustment ZnCl3-

dG_Chloride_Species[3] = (dG_Chloride_Species[3] + (C[1] * T * np.log(N) * 1E-03))          # Concentration (activity) adjustment ZnCl42-


dG_Sulfate_Species = np.array((-744.5, -755.9, -871.5))          # Thermodynamic data - (0) SO42-, (1) HSO4-, (2) ZnSO4

dG_Sulfate_Species[0] = (dG_Sulfate_Species[0] + (C[1] * T * np.log(SO42) * 1E-03))          # Concentration (activity) adjustment SO42-

dG_Sulfate_Species[1] = (dG_Sulfate_Species[1] + (C[1] * T * np.log(SO42) * 1E-03))          # Concentration (activity) adjustment HSO4-


dG_Carbonate_Species = np.array((-527.8, -586.8, -731.5))          # Thermodynamic data - (0) CO32-, (1) HCO3-, (2) ZnCO3

dG_Carbonate_Species[0] = (dG_Carbonate_Species[0] + (C[1] * T * np.log(CO32) * 1E-03))          # Concentration (activity) adjustment CO32-

dG_Carbonate_Species[1] = (dG_Carbonate_Species[1] + (C[1] * T * np.log(CO32) * 1E-03))          # Concentration (activity) adjustment HCO3-


## Gibbs free energy mimimization

dG_H2O_Grid = np.full_like(dG_Proton_Grid, dG_H2O)

dG_Zn_Grid_0 = np.full_like(dG_Proton_Grid, dG_Metal_Species[0])

dG_Zn_Grid_1 = np.full_like(dG_Proton_Grid, dG_Metal_Species[1])

dG_Zn_Grid_2 = np.full_like(dG_Proton_Grid, dG_Metal_Species[2])

dG_Zn_Grid_3 = np.full_like(dG_Proton_Grid, dG_Metal_Species[3])

dG_Zn_Grid_4 = np.full_like(dG_Proton_Grid, dG_Metal_Species[4])

dG_Zn_Grid_5 = np.full_like(dG_Proton_Grid, dG_Metal_Species[5])

dG_Zn_Grid_6 = np.full_like(dG_Proton_Grid, dG_Metal_Species[6])

dG_Zn_Grid_7 = np.full_like(dG_Proton_Grid, dG_Metal_Species[7])

dG_Zn_Grid_8 = np.full_like(dG_Proton_Grid, dG_Metal_Species[8])


### Chloride species

dG_Cl_Grid_0 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[0])

dG_Cl_Grid_1 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[1])

dG_Cl_Grid_2 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[2])

dG_Cl_Grid_3 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[3])

dG_Cl_Grid_4 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[4])

dG_Cl_Grid_5 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[5])


### Sulfate species

dG_SO4_Grid_0 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[0])

dG_SO4_Grid_1 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[1])

dG_SO4_Grid_2 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[2])


### Carbonate species

dG_CO3_Grid_0 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[0])

dG_CO3_Grid_1 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[1])

dG_CO3_Grid_2 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[2])


### Reaction 0 - Zn

dG_Reaction_0_Grid = (dG_Zn_Grid_0 - dG_Zn_Grid_0)

### Reaction 1 - Zn2+

dG_Reaction_1_Grid = ((dG_Zn_Grid_1 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0))

### Reaction 2 - ZnOH+

dG_Reaction_2_Grid = ((dG_Zn_Grid_2 + dG_Proton_Grid + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_H2O_Grid))

### Reaction 3 - ZnO22-

dG_Reaction_3_Grid = ((dG_Zn_Grid_3 + (4 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 4 - HZnO2-

dG_Reaction_4_Grid = ((dG_Zn_Grid_4 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 5 - Zn(OH)3-

dG_Reaction_5_Grid = ((dG_Zn_Grid_5 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (3 * dG_H2O_Grid)))

### Reaction 6 - Zn(OH)42-

dG_Reaction_6_Grid = ((dG_Zn_Grid_6 + (4 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (4 * dG_H2O_Grid)))

### Reaction 7 - Zn(OH)2

dG_Reaction_7_Grid = ((dG_Zn_Grid_7 + (2 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 8 - ZnO

dG_Reaction_8_Grid = ((dG_Zn_Grid_8 + (2 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_H2O_Grid))

### Reaction 9 - ZnCl+

dG_Reaction_9_Grid = ((dG_Cl_Grid_1 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_Cl_Grid_0))

### Reaction 10 - ZnCl3-

dG_Reaction_10_Grid = ((dG_Cl_Grid_2 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (3 * dG_Cl_Grid_0)))

### Reaction 11 - ZnCl42-

dG_Reaction_11_Grid = ((dG_Cl_Grid_3 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (4 * dG_Cl_Grid_0)))

### Reaction 12 - ZnCl2

dG_Reaction_12_Grid = ((dG_Cl_Grid_4 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + (2 * dG_Cl_Grid_0)))

### Reaction 13 - Zn(OH)Cl

dG_Reaction_13_Grid = ((dG_Cl_Grid_5 + dG_Proton_Grid + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_Cl_Grid_0 + dG_H2O_Grid))

### Reaction 14 - ZnSO4

dG_Reaction_14_Grid = ((dG_SO4_Grid_2 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_SO4_Grid_0))

### Reaction 15 - ZnCO3

dG_Reaction_15_Grid = ((dG_CO3_Grid_2 + (2 * dG_Electron_Grid)) - (dG_Zn_Grid_0 + dG_CO3_Grid_0))


### Minimum Gibbs free energy

dG_Energy_Minimization = np.stack([dG_Reaction_0_Grid, dG_Reaction_1_Grid, dG_Reaction_2_Grid, dG_Reaction_3_Grid, dG_Reaction_4_Grid, dG_Reaction_5_Grid, dG_Reaction_6_Grid, dG_Reaction_7_Grid, dG_Reaction_8_Grid, dG_Reaction_9_Grid, dG_Reaction_10_Grid, dG_Reaction_11_Grid, dG_Reaction_12_Grid, dG_Reaction_13_Grid, dG_Reaction_14_Grid, dG_Reaction_15_Grid], axis = 0)

dG_Pourbaix_Diagram = np.argmin(dG_Energy_Minimization, axis = 0)

np.savetxt(r'C:\\Users\\lukas.boonen\\OneDrive - Proviron Holding nv\\Desktop\\Research Papers\\Paper 1\\E_pH\\E_pH_Zn_V.csv', dG_Pourbaix_Diagram, delimiter = ';', comments='')

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

    ax1.plot(np.repeat(-3, 100), np.linspace(-2, 3, 100), '--', color = '#7E7E7E', alpha = 0.25)
    ax1.plot(np.repeat(1.92, 100), np.linspace(-2, 3, 100), '--', color = '#7E7E7E', alpha = 0.25)

    ax1.plot(np.repeat(6.4, 100), np.linspace(-2, 3, 100), '--', color = '#7E7E7E', alpha = 0.25)
    ax1.plot(np.repeat(10.4, 100), np.linspace(-2, 3, 100), '--', color = '#7E7E7E', alpha = 0.25)

    ax1.plot(np.linspace(8.5, 12, 100), np.repeat(-0.50, 100), '-', color = '#141414')
    ax1.plot(np.linspace(8.5, 14, 100), np.repeat(0.25, 100), '-', color = '#141414')
    ax1.plot(np.linspace(9.0, 10, 100), np.repeat(1.50, 100), '-', color = '#141414')

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

    plt.text(7.15, 1.85, r'N', fontsize = 6)

    plt.text(0.60, -1.90, r'HSO$_{4}^{-}$', fontsize = 6)
    plt.text(2.05, -1.90, r'SO$_{4}^{2-}$', fontsize = 6)

    plt.text(5.45, -1.90, r'CO$_{2}$', fontsize = 6)
    plt.text(6.55, -1.90, r'HCO$_{3}^{-}$', fontsize = 6)

    plt.text(9.05, -1.90, r'HCO$_{3}^{-}$', fontsize = 6)
    plt.text(10.50, -1.90, r'CO$_{3}^{2-}$', fontsize = 6)

    plt.text(0.00, 1.30, r'O$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00, 1.15, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')

    plt.text(0.00, 0.08, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00,-0.08, r'H$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    
    plt.text(4.00, -1.50, r'Zn$_{(s)}$', fontsize = 10)
    plt.text(1.25, 1.50, r'ZnCO$_{3 (s)}$', fontsize = 10)
    plt.text(6.75, 1.50, r'ZnO$_{(s)}$', fontsize = 10)
    plt.text(5.00, -0.50, r'HZnO$_{2 (aq)}^{-}$', fontsize = 10)
    plt.text(4.50, 0.25, r'Zn(OH)$_{4 (aq)}^{2-}$', fontsize = 10)

    fig_Pourbaix_Diagram.tight_layout()
Plot_Pourbaix_Diagram()

def Export_Plot(fig, filename, dpi):
    fig.savefig(f"{filename}.png", dpi = dpi)
    fig.savefig(f"{filename}.jpg", dpi = dpi)
    fig.savefig(f"{filename}.pdf")
    fig.savefig(f"{filename}.svg")

Export_Plot(fig_Pourbaix_Diagram, r'C:\\Users\\lukas.boonen\\OneDrive - Proviron Holding nv\\Desktop\\Research Papers\\Paper 1\\E_pH\\E_pH_Zn_V', 1200)

plt.show()