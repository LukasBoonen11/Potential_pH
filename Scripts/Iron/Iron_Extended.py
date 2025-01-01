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


O = np.array((1E-09))

SO42 = np.array((0.75E-03))


CO32 = np.array((1.15E-03))


dG_H2O = np.array((-237.1))         # Thermodynamic data - (0) H2O

dG_Metal_Species = np.array((0, -78.9, -4.7, -295.3, -377.7, -438.0, -229.4, -614.9, -769.7, -742.2, -1015.4, -486.5))            # Thermodynamic data - (0) Fe, (1) Fe2+, (2) Fe3+, (3) FeO22-, (4) HFeO2-, (5) Fe(OH)2+, (6) FeOH2+, (7) Fe(OH)3-, (8) Fe(OH)42-, (9) Fe2O3, (10) Fe3O4, (11) Fe(OH)2

dG_Metal_Species[1] = (dG_Metal_Species[1] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[2] = (dG_Metal_Species[2] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[3] = (dG_Metal_Species[3] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[4] = (dG_Metal_Species[4] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[5] = (dG_Metal_Species[5] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[6] = (dG_Metal_Species[6] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe2+

dG_Metal_Species[7] = (dG_Metal_Species[7] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe(OH)3-

dG_Metal_Species[8] = (dG_Metal_Species[8] + (C[1] * T * np.log(M) * 1E-03))          # Concentration (activity) adjustment Fe(OH)42-


dG_Chloride_Species = np.array((-131.2, -143.9, -302.3, -334.0))          # Thermodynamic data - (0) Cl-, (1) FeCl2+, (2) FeCl2, (3) FeCl3

dG_Chloride_Species[0] = (dG_Chloride_Species[0] + (C[1] * T * np.log(Cl) * 1E-03))          # Concentration (activity) adjustment Cl-

dG_Chloride_Species[1] = (dG_Chloride_Species[1] + (C[1] * T * np.log(N) * 1E-03))          # Concentration (activity) adjustment FeCl2+


dG_Sulfate_Species = np.array((-744.5, -755.9, -772.7, -1524.5, -820.8))          # Thermodynamic data - (0) SO42-, (1) HSO4-, (2) FeSO4+, (3) Fe(SO4)2-, (4) FeSO4

dG_Sulfate_Species[0] = (dG_Sulfate_Species[0] + (C[1] * T * np.log(SO42) * 1E-03))          # Concentration (activity) adjustment SO42-

dG_Sulfate_Species[1] = (dG_Sulfate_Species[1] + (C[1] * T * np.log(SO42) * 1E-03))          # Concentration (activity) adjustment HSO4-

dG_Sulfate_Species[2] = (dG_Sulfate_Species[2] + (C[1] * T * np.log(O) * 1E-03))          # Concentration (activity) adjustment FeSO4+

dG_Sulfate_Species[3] = (dG_Sulfate_Species[3] + (C[1] * T * np.log(O) * 1E-03))          # Concentration (activity) adjustment Fe(SO4)2-


dG_Carbonate_Species = np.array((-527.8, -586.8, -731.5))          # Thermodynamic data - (0) CO32-, (1) HCO3-, (2) FeCO3

dG_Carbonate_Species[0] = (dG_Carbonate_Species[0] + (C[1] * T * np.log(CO32) * 1E-03))          # Concentration (activity) adjustment CO32-

dG_Carbonate_Species[1] = (dG_Carbonate_Species[1] + (C[1] * T * np.log(CO32) * 1E-03))          # Concentration (activity) adjustment HCO3-


## Gibbs free energy mimimization

dG_H2O_Grid = np.full_like(dG_Proton_Grid, dG_H2O)

dG_Fe_Grid_0 = np.full_like(dG_Proton_Grid, dG_Metal_Species[0])

dG_Fe_Grid_1 = np.full_like(dG_Proton_Grid, dG_Metal_Species[1])

dG_Fe_Grid_2 = np.full_like(dG_Proton_Grid, dG_Metal_Species[2])

dG_Fe_Grid_3 = np.full_like(dG_Proton_Grid, dG_Metal_Species[3])

dG_Fe_Grid_4 = np.full_like(dG_Proton_Grid, dG_Metal_Species[4])

dG_Fe_Grid_5 = np.full_like(dG_Proton_Grid, dG_Metal_Species[5])

dG_Fe_Grid_6 = np.full_like(dG_Proton_Grid, dG_Metal_Species[6])

dG_Fe_Grid_7 = np.full_like(dG_Proton_Grid, dG_Metal_Species[7])

dG_Fe_Grid_8 = np.full_like(dG_Proton_Grid, dG_Metal_Species[8])

dG_Fe_Grid_9 = np.full_like(dG_Proton_Grid, dG_Metal_Species[9])

dG_Fe_Grid_10 = np.full_like(dG_Proton_Grid, dG_Metal_Species[10])

dG_Fe_Grid_11 = np.full_like(dG_Proton_Grid, dG_Metal_Species[11])


### Chloride species

dG_Cl_Grid_0 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[0])

dG_Cl_Grid_1 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[1])

dG_Cl_Grid_2 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[2])

dG_Cl_Grid_3 = np.full_like(dG_Proton_Grid, dG_Chloride_Species[3])


### Sulfate species

dG_SO4_Grid_0 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[0])

dG_SO4_Grid_1 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[1])

dG_SO4_Grid_2 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[2])

dG_SO4_Grid_3 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[3])

dG_SO4_Grid_4 = np.full_like(dG_Proton_Grid, dG_Sulfate_Species[4])


### Carbonate species

dG_CO3_Grid_0 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[0])

dG_CO3_Grid_1 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[1])

dG_CO3_Grid_2 = np.full_like(dG_Proton_Grid, dG_Carbonate_Species[2])


### Reaction 0 - Fe

dG_Reaction_0_Grid = (dG_Fe_Grid_0 - dG_Fe_Grid_0)

### Reaction 1 - Fe2+

dG_Reaction_1_Grid = ((dG_Fe_Grid_1 + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0))

### Reaction 2 - Fe3+

dG_Reaction_2_Grid = ((dG_Fe_Grid_2 + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0))

### Reaction 3 - FeO22-

dG_Reaction_3_Grid = ((dG_Fe_Grid_3 + (4 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 4 - HFeO2-

dG_Reaction_4_Grid = ((dG_Fe_Grid_4 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 5 - Fe(OH)2+

dG_Reaction_5_Grid = ((dG_Fe_Grid_5 + (2 * dG_Proton_Grid) + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 6 - FeOH2+

dG_Reaction_6_Grid = ((dG_Fe_Grid_6 + dG_Proton_Grid + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + dG_H2O_Grid))

### Reaction 7 - Fe(OH)3-

dG_Reaction_7_Grid = ((dG_Fe_Grid_7 + (3 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (3 * dG_H2O_Grid)))

### Reaction 8 - Fe(OH)42-

dG_Reaction_8_Grid = ((dG_Fe_Grid_8 + (4 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (4 * dG_H2O_Grid)))

### Reaction 9 - Fe2O3

dG_Reaction_9_Grid = (((1/2 * dG_Fe_Grid_9) + (3 * dG_Proton_Grid) + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (3/2 * dG_H2O_Grid)))

### Reaction 10 - Fe3O4

dG_Reaction_10_Grid = (((1/3 * dG_Fe_Grid_10) + (8/3 * dG_Proton_Grid) + (8/3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (4/3 * dG_H2O_Grid)))

### Reaction 11 - Fe(OH)2

dG_Reaction_11_Grid = ((dG_Fe_Grid_11 + (2 * dG_Proton_Grid) + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_H2O_Grid)))

### Reaction 12 - FeCl2+

dG_Reaction_12_Grid = ((dG_Cl_Grid_1 + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + dG_Cl_Grid_0))

### Reaction 13 - FeCl2

dG_Reaction_13_Grid = ((dG_Cl_Grid_2 + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_Cl_Grid_0)))

### Reaction 14 - FeCl3

dG_Reaction_14_Grid = ((dG_Cl_Grid_3 + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (3 * dG_Cl_Grid_0)))

### Reaction 15 - FeSO4+

dG_Reaction_15_Grid = ((dG_SO4_Grid_2 + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + dG_SO4_Grid_0))

### Reaction 16 - Fe(SO4)2-

dG_Reaction_16_Grid = ((dG_SO4_Grid_3 + (3 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + (2 * dG_SO4_Grid_0)))

### Reaction 17 - FeSO4

dG_Reaction_17_Grid = ((dG_SO4_Grid_4 + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + dG_SO4_Grid_0))

### Reaction 18 - FeCO3

dG_Reaction_18_Grid = ((dG_CO3_Grid_2 + (2 * dG_Electron_Grid)) - (dG_Fe_Grid_0 + dG_CO3_Grid_0))


### Minimum Gibbs free energy

dG_Energy_Minimization = np.stack([dG_Reaction_0_Grid, dG_Reaction_1_Grid, dG_Reaction_2_Grid, dG_Reaction_3_Grid, dG_Reaction_4_Grid, dG_Reaction_5_Grid, dG_Reaction_6_Grid, dG_Reaction_7_Grid, dG_Reaction_8_Grid, dG_Reaction_9_Grid, dG_Reaction_10_Grid, dG_Reaction_11_Grid, dG_Reaction_12_Grid, dG_Reaction_13_Grid, dG_Reaction_14_Grid, dG_Reaction_15_Grid, dG_Reaction_16_Grid, dG_Reaction_17_Grid, dG_Reaction_18_Grid], axis = 0)

dG_Pourbaix_Diagram = np.argmin(dG_Energy_Minimization, axis = 0)

np.savetxt(r'C:\\Users\\lukas.boonen\\OneDrive - Proviron Holding nv\\Desktop\\Research Papers\\Paper 1\\E_pH\\E_pH_Fe_V.csv', dG_Pourbaix_Diagram, delimiter = ';', comments='')

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

    ax1.plot(np.repeat(14.3, 100), np.linspace(-0.75, -0.25, 100), '-', color = '#141414')
    ax1.plot(np.repeat(15.0, 100), np.linspace(-1.20, -0.90, 100), '-', color = '#141414')

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

    plt.text(0.00, 1.40, r'O$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00, 1.15, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')

    plt.text(0.00, 0.08, r'H$_{2}$O', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    plt.text(0.00,-0.08, r'H$_{2}$', fontsize = 6, rotation = -15, ha = 'center', va = 'center')
    
    plt.text(4.00, -1.50, r'Fe$_{(s)}$', fontsize = 9)
    plt.text(3.00, 0.00, r'FeCO$_{3 (s)}$', fontsize = 9)
    plt.text(-1.00, 1.65, r'FeSO$_{4 (aq)}^{+}$', fontsize = 9)
    plt.text(13.00, -0.15, r'Fe$_{3}$O$_{4 (s)}$', fontsize = 9)
    plt.text(11.00, 1.00, r'Fe$_{2}$O$_{3 (s)}$', fontsize = 9)
    plt.text(12.25, -1.50, r'Fe(OH)$_{4 (aq)}^{2-}$', fontsize = 9)

    fig_Pourbaix_Diagram.tight_layout()
Plot_Pourbaix_Diagram()

def Export_Plot(fig, filename, dpi):
    fig.savefig(f"{filename}.png", dpi = dpi)
    fig.savefig(f"{filename}.jpg", dpi = dpi)
    fig.savefig(f"{filename}.pdf")
    fig.savefig(f"{filename}.svg")

Export_Plot(fig_Pourbaix_Diagram, r'C:\\Users\\lukas.boonen\\OneDrive - Proviron Holding nv\\Desktop\\Research Papers\\Paper 1\\E_pH\\E_pH_Fe_V', 1200)

plt.show()