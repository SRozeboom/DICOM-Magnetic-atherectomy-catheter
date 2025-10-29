#!/usr/bin/env python3

import numpy as np

# Cylindrical rod
class Segment:
    # Units in meter, Pascal, kg/m3
    def __init__(self, length, 
                 radius_inner_start, radius_inner_end,
                 radius_outer_start, radius_outer_end,
                 young_modulus, poisson_ratio, density, 
                 precurvature):
        # Assign member variables
        self.length = length
        self.radius_inner_start = radius_inner_start
        self.radius_inner_end = radius_inner_end
        self.radius_outer_start = radius_outer_start
        self.radius_outer_end = radius_outer_end
        self.E = young_modulus
        self.G = self.E / (2 * (1 + poisson_ratio))
        self.density = density
        self.precurvature = np.array(precurvature).reshape(3,1)
        # Give zero magnetization
        self.M = np.array([0,0,0]).reshape(3,1)
        self.Br = 0

    def getStiffnessAtPosition(self, s):
        """
        if s > self.length:
            raise Exception("Centerline position exceeds segment length")
        else:
        """
        # Get radii at centerline position
        radius_outer = self.radius_outer_start + s * (self.radius_outer_end - self.radius_outer_start) / self.length
        radius_inner = self.radius_inner_start + s * (self.radius_inner_end - self.radius_inner_start) / self.length
        # Compute area of cross-section
        A = np.pi * (radius_outer**2 - radius_inner**2)
        # Compute moments of area
        I = (np.pi / 4) * (radius_outer**4 - radius_inner**4)
        J = 2 * I
        # Compute stiffness matrices
        Kse = np.diag([self.G*A, self.G*A, self.E*A])
        Kbe = np.diag([self.E*I, self.E*I, self.G*J])
        # Return
        return I, J, A, Kse, Kbe


# Magnetic Segment
class MagnetSegment(Segment):
    def __init__(self, magnetization_axis, magnetization_coeff, remanence, **kwargs):
        # Pass and assign base class parameters
        super().__init__(**kwargs)
        # Construct magnetization vector
        if magnetization_axis == 'x':
            self.M = magnetization_coeff * np.array([1,0,0]).reshape(-1,1)
        elif magnetization_axis == 'y':
            self.M = magnetization_coeff * np.array([0,1,0]).reshape(-1,1)
        elif magnetization_axis == 'z':
            self.M = magnetization_coeff * np.array([0,0,1]).reshape(-1,1)
        
        # Assign remanence
        self.Br = remanence


# Passive segment
class PassiveSegment(Segment):
    pass