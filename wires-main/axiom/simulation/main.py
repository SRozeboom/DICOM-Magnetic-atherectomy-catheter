#!/usr/bin/env python3

import subprocess
import sys

def install_and_import(package, import_name=None):
    import_name = import_name or package
    try:
        __import__(import_name)
    except ModuleNotFoundError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    finally:
        globals()[import_name] = __import__(import_name)

install_and_import("numpy")
install_and_import("matplotlib")
install_and_import("scipy")

from segments import MagnetSegment, PassiveSegment
from solver import Solver
import numpy as np
import matplotlib.pyplot as plt

def create_rod(dc1_value, remanence_value=1.42):
    # Segment lengths
    L0, L1, L2, L3, L4, L5 = 10, 80, 0, 3, 12, 3

    # Core diameters
    dc0 = 0.36
    dc1 = dc1_value

    # Sleeve diameters
    ds0, ds3, ds5 = 0.89, 3, 2

    # Magnetic field
    bs = 35  # mT

    # Core properties
    rod_core_sections = 6
    rod_core_section_magnetic = np.zeros(rod_core_sections)
    rod_core_section_magnetization_axis = ['0']*rod_core_sections
    rod_core_section_magnetization_axis_coefficient = np.zeros(rod_core_sections)
    rod_core_radius_start = 0.5 * np.array([dc0, dc0, dc1, dc1, dc1, dc1]) * 1e-3
    rod_core_radius_endings = 0.5 * np.array([dc0, dc1, dc1, dc1, dc1, dc1]) * 1e-3
    rod_core_lengths = np.array([L0, L1, L2, L3, L4, L5]) * 1e-3
    rod_core_precurvature = np.zeros((3, rod_core_sections))
    rod_core_poisson_ratio = 0.3

    rod_core_density = []
    rod_core_young_modulus = []
    for i in range(rod_core_sections):
        if rod_core_section_magnetic[i]:
            rod_core_density.append(8000.0)
            rod_core_young_modulus.append(200e9)
        else:
            rod_core_density.append(6450.0)
            rod_core_young_modulus.append(50e9)

    # Sleeve properties
    rod_sleeve_sections = rod_core_sections
    rod_sleeve_section_magnetic = np.array([0,0,0,1,0,1])
    rod_sleeve_section_magnetization_axis = ['0','0','0','x','0','z']
    rod_sleeve_section_magnetization_axis_coefficient = np.array([0,0,0,1,0,-1])
    rod_sleeve_radius_start = 0.5 * np.array([ds0, ds0, ds0, ds3, ds0, ds5]) * 1e-3
    rod_sleeve_radius_endings = rod_sleeve_radius_start.copy()
    rod_sleeve_lengths = rod_core_lengths
    rod_sleeve_precurvature = np.zeros((3, rod_core_sections))
    rod_sleeve_poisson_ratio = 0.3

    rod_sleeve_density = []
    rod_sleeve_young_modulus = []
    for i in range(rod_sleeve_sections):
        if rod_sleeve_section_magnetic[i]:
            rod_sleeve_density.append(1200)
            rod_sleeve_young_modulus.append(200e9)
        else:
            rod_sleeve_density.append(8000)
            rod_sleeve_young_modulus.append(100e6)

    # Generate core and sleeve
    rod_core = []
    rod_sleeve = []
    for i in range(rod_core_sections):
        # Core
        if rod_core_section_magnetic[i]:
            rod_core.append(MagnetSegment(
                magnetization_axis=rod_core_section_magnetization_axis[i],
                magnetization_coeff=rod_core_section_magnetization_axis_coefficient[i],
                remanence=remanence_value,
                length=rod_core_lengths[i],
                radius_inner_start=0.0,
                radius_inner_end=0.0,
                radius_outer_start=rod_core_radius_start[i],
                radius_outer_end=rod_core_radius_endings[i],
                young_modulus=rod_core_young_modulus[i],
                poisson_ratio=rod_core_poisson_ratio,
                density=rod_core_density[i],
                precurvature=rod_core_precurvature[:,i]))
        else:
            rod_core.append(PassiveSegment(
                length=rod_core_lengths[i],
                radius_inner_start=0.0,
                radius_inner_end=0.0,
                radius_outer_start=rod_core_radius_start[i],
                radius_outer_end=rod_core_radius_endings[i],
                young_modulus=rod_core_young_modulus[i],
                poisson_ratio=rod_core_poisson_ratio,
                density=rod_core_density[i],
                precurvature=rod_core_precurvature[:,i]))
        # Sleeve
        if rod_sleeve_section_magnetic[i]:
            rod_sleeve.append(MagnetSegment(
                magnetization_axis=rod_sleeve_section_magnetization_axis[i],
                magnetization_coeff=rod_sleeve_section_magnetization_axis_coefficient[i],
                remanence=remanence_value,
                length=rod_sleeve_lengths[i],
                radius_inner_start=rod_core_radius_start[i],
                radius_inner_end=rod_core_radius_endings[i],
                radius_outer_start=rod_sleeve_radius_start[i],
                radius_outer_end=rod_sleeve_radius_endings[i],
                young_modulus=rod_sleeve_young_modulus[i],
                poisson_ratio=rod_sleeve_poisson_ratio,
                density=rod_sleeve_density[i],
                precurvature=rod_sleeve_precurvature[:,i]))
        else:
            rod_sleeve.append(PassiveSegment(
                length=rod_sleeve_lengths[i],
                radius_inner_start=rod_core_radius_start[i],
                radius_inner_end=rod_core_radius_endings[i],
                radius_outer_start=rod_sleeve_radius_start[i],
                radius_outer_end=rod_sleeve_radius_endings[i],
                young_modulus=rod_sleeve_young_modulus[i],
                poisson_ratio=rod_sleeve_poisson_ratio,
                density=rod_sleeve_density[i],
                precurvature=rod_sleeve_precurvature[:,i]))
    
    return np.array([rod_core, rod_sleeve]), bs, rod_core_precurvature

def main():
    # Both rods use the magnetic moment of the 0.9mm rod
    rods1, bs1, curv1 = create_rod(0.3556, remanence_value=1.42)
    rods2, bs2, curv2 = create_rod(0.9, remanence_value=1.42)

    # Solve both rods
    solver1 = Solver(rods=rods1, gravity=np.array([0,0,0]), field_T=np.array([0,0,bs1*1e-3]))
    Y1 = solver1.solve(p1=np.array([0,0,0]), q1=np.array([0,0,0,1]), guess=np.zeros(5+len(rods1)), thetas=np.zeros(len(rods1)))

    solver2 = Solver(rods=rods2, gravity=np.array([0,0,0]), field_T=np.array([0,0,bs2*1e-3]))
    Y2 = solver2.solve(p1=np.array([0,0,0]), q1=np.array([0,0,0,1]), guess=np.zeros(5+len(rods2)), thetas=np.zeros(len(rods2)))

    # Plot both rods with curvature
    ax = plt.axes(projection='3d')

    # Rod 1 (guidewire)
    for i in range(len(Y1)):
        ax.plot3D(Y1[i][0,:]*1e3, Y1[i][1,:]*1e3, Y1[i][2,:]*1e3,
                  color='blue', linewidth=4, label='0.014" guidewire' if i==0 else "")
        # plot curvature as small arrows along the rod
        for j in range(0, len(Y1[i][0,:]), 10):
            dx, dy, dz = curv1[:,i] * 1e3
            ax.quiver(Y1[i][0,j]*1e3, Y1[i][1,j]*1e3, Y1[i][2,j]*1e3, dx, dy, dz, color='cyan', length=2, normalize=True)

    # Rod 2 (catheter)
    for i in range(len(Y2)):
        ax.plot3D(Y2[i][0,:]*1e3, Y2[i][1,:]*1e3, Y2[i][2,:]*1e3,
                  color='red', linewidth=4, label='0.9 mm catheter' if i==0 else "")
        for j in range(0, len(Y2[i][0,:]), 10):
            dx, dy, dz = curv2[:,i] * 1e3
            ax.quiver(Y2[i][0,j]*1e3, Y2[i][1,j]*1e3, Y2[i][2,j]*1e3, dx, dy, dz, color='magenta', length=2, normalize=True)

    ax.set_xlim([-100,100])
    ax.set_ylim([-100,100])
    ax.set_zlabel('z (mm)')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
