#!/usr/bin/env python3

import subprocess
import sys

def install_and_import(package, import_name=None):
    """
    Try to import a package. If missing, install it automatically.
    """
    import_name = import_name or package
    try:
        __import__(import_name)
    except ModuleNotFoundError:
        print(f"Module '{import_name}' not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"Installed '{package}' successfully.")
    finally:
        globals()[import_name] = __import__(import_name)

# 🔹 Try to import required packages or install them if missing
install_and_import("numpy")
install_and_import("matplotlib")
install_and_import("scipy")        # if Solver or segments rely on it
# Custom modules (segments, solver) won't install via pip — they must exist in directory

from segments import MagnetSegment, PassiveSegment
from solver import Solver
import numpy as np
import matplotlib.pyplot as plt

from segments import MagnetSegment
from segments import PassiveSegment
from solver import Solver
import numpy as np
import matplotlib.pyplot as plt

def main(args=None):

    # Parameterization length (mm). Magnets: L3 & L5, Taper" L1
    L0 = 10
    L1 = 80
    L2 = 0
    L3 = 3
    L4 = 12
    L5 = 3

    # Parameterization diameter core (mm)
    dc0 = 0.36
    dc1 = 0.09

    # Parameterization diameter sleeve (mm)
    ds0 = 0.89
    ds3 = 3
    ds5 = 2

    # Field strength (mT)
    bs = 10
    
    # ROD CORE PROPERTIES
    rod_core_sections = 6
    rod_core_section_magnetic = np.zeros(rod_core_sections) # 0 = passive, 1 = magnetic
    rod_core_section_magnetization_axis = ['0', '0', '0', '0', '0', '0']
    rod_core_section_magnetization_axis_coefficient = np.array([0,0,0,0,0,0])
    rod_core_radius_start =     0.5 * np.array([dc0, dc0, dc1, dc1, dc1, dc1]) * 1.0e-3 # m
    rod_core_radius_endings =   0.5 * np.array([dc0, dc1, dc1, dc1, dc1, dc1]) * 1.0e-3 #m
    rod_core_lengths = np.array([L0, L1, L2, L3, L4, L5]) * 1.0e-3 # m
    rod_core_precurvature = np.zeros((3,rod_core_sections)) # Straight tip
    rod_core_poisson_ratio = 0.3
    rod_core_density = []
    rod_core_young_modulus = []
    for i in range(rod_core_sections):
        if rod_core_section_magnetic[i]:
            rod_core_density.append(8000.0) # kg / m3
            rod_core_young_modulus.append(200.0e9) # Pa
        else:
            rod_core_density.append(6450.0) # kg / m3
            rod_core_young_modulus.append(50.0e9) # Pa

    # ROD SLEEVE PROPERTIES
    rod_sleeve_sections = rod_core_sections
    rod_sleeve_section_magnetic = np.array([0, 0, 0, 1, 0, 1]) # 0 = passive, 1 = magnetic
    rod_sleeve_section_magnetization_axis = ['0', '0', '0', 'x', '0', 'z']
    rod_sleeve_section_magnetization_axis_coefficient = np.array([0,0,0, 1, 0,-1])
    rod_sleeve_radius_start = 0.5 * np.array([ds0, ds0, ds0, ds3, ds0, ds5]) * 1e-3
    rod_sleeve_radius_endings = 0.5 * np.array([ds0, ds0, ds0, ds3, ds0, ds5]) * 1e-3
    rod_sleeve_lengths = rod_core_lengths
    rod_sleeve_precurvature = np.zeros((3,rod_core_sections)) # Straight tip
    rod_sleeve_poisson_ratio = 0.3
    rod_sleeve_density = []
    rod_sleeve_young_modulus = []
    for i in range(rod_sleeve_sections):
        if rod_sleeve_section_magnetic[i]:
            rod_sleeve_density.append(1200) # kg / m3
            rod_sleeve_young_modulus.append(200e9) # Pa
        else:
            rod_sleeve_density.append(8000) # kg / m3
            rod_sleeve_young_modulus.append(100e6) # Pa


    # CHECK
    if (len(rod_core_section_magnetic) != rod_core_sections):
        raise Exception("len(rod_core_section_magnetic) !=  rod_core_sections")
    elif (len(rod_core_radius_start) != rod_core_sections):
        raise Exception("len(rod_core_radius_start) != rod_core_sections")
    elif (len(rod_core_radius_endings) != rod_core_sections):
        raise Exception("len(rod_core_radius_endings) != rod_core_sections")
    elif (len(rod_core_lengths) != rod_core_sections):
        raise Exception("len(rod_core_lengths) != rod_core_sections")
    elif (len(rod_sleeve_section_magnetic) != rod_core_sections):
        raise Exception("len(rod_sleeve_section_magnetic) != rod_core_sections")
    elif (len(rod_sleeve_radius_start) != rod_core_sections):
        raise Exception("len(rod_sleeve_radius_start) != rod_core_sections")
    elif (len(rod_sleeve_radius_endings) != rod_core_sections):
        raise Exception("len(rod_sleeve_radius_endings) != rod_core_sections")
    

    # GENERATE RODS
    rod_core = []
    rod_sleeve = []
    for i in range(rod_core_sections):
        # Core
        if rod_core_section_magnetic[i]:
            rod_core.append( MagnetSegment(magnetization_axis   =   rod_core_section_magnetization_axis[i],
                                           magnetization_coeff  =   rod_core_section_magnetization_axis_coefficient[i], 
                                           remanence            =   1.42, 
                                           length               =   rod_core_lengths[i], 
                                           radius_inner_start   =   0.0, 
                                           radius_inner_end     =   0.0,
                                           radius_outer_start   =   rod_core_radius_start[i], 
                                           radius_outer_end     =   rod_core_radius_endings[i],
                                           young_modulus        =   rod_core_young_modulus[i], 
                                           poisson_ratio        =   rod_core_poisson_ratio, 
                                           density              =   rod_core_density[i],
                                           precurvature         =   rod_core_precurvature[:,i]) )
        else:
            rod_core.append( PassiveSegment(length               =   rod_core_lengths[i], 
                                            radius_inner_start   =   0.0, 
                                            radius_inner_end     =   0.0,
                                            radius_outer_start   =   rod_core_radius_start[i], 
                                            radius_outer_end     =   rod_core_radius_endings[i],
                                            young_modulus        =   rod_core_young_modulus[i], 
                                            poisson_ratio        =   rod_core_poisson_ratio, 
                                            density              =   rod_core_density[i],
                                            precurvature         =   rod_core_precurvature[:,i]) ) 
        
        # Sleeve
        if rod_sleeve_section_magnetic[i]:
            rod_sleeve.append( MagnetSegment(magnetization_axis   =   rod_sleeve_section_magnetization_axis[i],
                                             magnetization_coeff  =   rod_sleeve_section_magnetization_axis_coefficient[i], 
                                             remanence            =   1.42, 
                                             length               =   rod_sleeve_lengths[i], 
                                             radius_inner_start   =   rod_core_radius_start[i], 
                                             radius_inner_end     =   rod_core_radius_endings[i],
                                             radius_outer_start   =   rod_sleeve_radius_start[i], 
                                             radius_outer_end     =   rod_sleeve_radius_endings[i],
                                             young_modulus        =   rod_sleeve_young_modulus[i], 
                                             poisson_ratio        =   rod_sleeve_poisson_ratio, 
                                             density              =   rod_sleeve_density[i],
                                             precurvature         =   rod_sleeve_precurvature[:,i]) )
        else:
            rod_sleeve.append( PassiveSegment(length               =   rod_sleeve_lengths[i], 
                                              radius_inner_start   =   rod_core_radius_start[i], 
                                              radius_inner_end     =   rod_core_radius_endings[i],
                                              radius_outer_start   =   rod_sleeve_radius_start[i], 
                                              radius_outer_end     =   rod_sleeve_radius_endings[i],
                                              young_modulus        =   rod_sleeve_young_modulus[i], 
                                              poisson_ratio        =   rod_sleeve_poisson_ratio, 
                                              density              =   rod_sleeve_density[i],
                                              precurvature         =   rod_sleeve_precurvature[:,i]) ) 
    
    rods = np.array([rod_core, rod_sleeve])



    # INITIALIZE SOLVER
    ctr_solver = Solver(rods    =   rods, 
                        gravity =   np.array([0,0,0]), 
                        field_T =   np.array([0,0,bs*1e-3]))

    # COMPUTE SOLUTION
    Y = ctr_solver.solve(p1     =   np.array([0,0,0]), 
                         q1     =   np.array([0,0,0,1]), 
                         guess  =   np.zeros(5+len(rods)),
                         thetas =   np.zeros(len(rods)))
    

    # Plot
    ax = plt.axes(projection='3d')
    for i in range(len(Y)):
        if rod_sleeve_section_magnetic[i] == 0: # Plot blue
            ax.plot3D(Y[i][0, :] * 1e3, Y[i][1, :] * 1e3, Y[i][2, :] * 1e3, color='blue', linewidth=rod_sleeve_radius_endings[i]*4*1e3)
        else: # Plot red
            ax.plot3D(Y[i][0, :] * 1e3, Y[i][1, :] * 1e3, Y[i][2, :] * 1e3, color='red', linewidth=rod_sleeve_radius_endings[i]*4*1e3)

    ax.set_xlim([-100, 100])
    ax.set_ylim([-100, 100])
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_zlabel('z (mm)')
    plt.show()


if __name__ == "__main__":
    main()

