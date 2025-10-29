#!/usr/bin/env python3

import numpy as np
import scipy.optimize
from scipy import integrate
from scipy.spatial.transform import Rotation as R

class Solver:
    
    def __init__(self, rods, gravity, field_T):
        self.rods = rods
        self.gravity = 9.81 * np.array(gravity).reshape(-1,1)
        self.b_T_target = field_T
    

    def solve(self, p1, q1, guess, thetas):   

        # Check vector compatibility
        if len(thetas) != len(self.rods):
            raise Exception("Incompatible number of rod body rotations")
        
        
        # Extract and convert to np.array
        self.p      = p1
        self.q      = q1
        init_guess  = guess # n, mbxy, mb1z, mb2z, ...
        self.thetas = thetas

        # Initialize solution matrix
        self.Y      = []

        # Prepare to loop/solve over field magnitudes
        db_T            = 0.5e-3
        b_T_target_str  = np.linalg.norm(self.b_T_target)
        b_T_strs        = np.linspace(0, b_T_target_str, int(np.floor(b_T_target_str/db_T))+1)
        b_T_dir         = self.b_T_target / b_T_target_str

        # Solve
        print("Progress:")
        for b_T_str in b_T_strs:
            # Report
            print(str(np.where(b_T_strs == b_T_str)[0]+1) + " / " + str(len(b_T_strs)))
            # Assign field
            self.b_T = b_T_str * b_T_dir
            # Compute solution
            sol = scipy.optimize.root(self.obj_fun, init_guess, method="lm").x
            # Pass solution forward
            init_guess = sol

        # Return
        return self.Y
    

    def obj_fun(self, guess):

        # Compose proximal state vector
        y0 = np.concatenate((self.p, self.q, guess, self.thetas))

        # Integrate along rods (rods[0][:] inner rod, rods[1][:] first layer, rods[2][:] second layer, etc.)
        Y = []
        for i in range(len(self.rods[0,:])):
            # Get i-th set of concentric tube robot (ctr) segments
            ctr_set = self.rods[:,i]
            # Get length over which to integrate
            L = ctr_set[0].length
            # Get centerline positions
            self.ds=0.1e-3
            if L == 0:
                Yi = y0.reshape(-1,1)
            else:
                if L < self.ds:
                    self.ds = L
                s = np.linspace(0, L, int(np.floor(L/self.ds))+1)
                # Forward integrate
                Yi = integrate.solve_ivp(self.solver_ode, t_span=(0,L), y0=y0, args=(ctr_set,), t_eval=s, method='RK45').y
                # Extract updated initial state vector
                y0 = Yi[:,-1].flatten()
            # Append integration to solution
            Y.append(Yi)

        self.Y = Y

        # Compute distal boundary constraints
        guessL = Y[-1][7:-len(self.thetas),-1]
        target = np.zeros(len(guess))
        error = target - guessL
        return error
        
    

    def solver_ode(self, s, y, ctr_set):
        # Get number of concentric tube robots
        N = len(ctr_set)

        # Unpack state vector and reshape to column vectors
        p1 = y[0:3]
        q1 = y[3:7]
        #q1 = q1 / np.linalg.norm(q1)
        Rq1 = R.from_quat(q1.flatten()).as_matrix()
        n = y[7:10]
        mbxy = y[10:12]
        mbz = y[12:12+N]
        thetas = y[-N:]
        
        # --------- Compute curvature vectors ---------
        u1xy_num = mbxy.reshape(-1,1)
        u1xy_den = 0
        uiz = np.zeros(N)
        dthetas = np.zeros((N,1))
        ui = np.zeros((3,N))
        
        for i in range(N):
            # Get stiffness values at centerline position s
            Is, Js, *_ = ctr_set[i].getStiffnessAtPosition(s)
            # Get body rotation
            Rz = R.from_euler('z', thetas[i], degrees=False).as_matrix()
            # Inner rod curvature numerator and denominator (eqn. 5a)
            #u1xy_num += + ctr_set[i].E * Is * self.rotation_matrix_z(thetas[i])[:2, :2] * ctr_set[i].precurvature[:2]
            u1xy_num = u1xy_num + ctr_set[i].E * Is * Rz[:2,:2] @ ctr_set[i].precurvature[:2]
            u1xy_den = u1xy_den + ctr_set[i].E * Is
            # i-th rod polar curvature (eqn. 5b)
            uiz[i] = (mbz[i] + ctr_set[i].G * Js * ctr_set[i].precurvature[-1:]) / (ctr_set[i].G * Js)
            # i-th rod body rotation derivative (eqn. 4f)
            dthetas[i] = uiz[i] - uiz[0]


        # Compose curvature vectors of rod (eqn. 5c)
        u1xy = u1xy_num / u1xy_den
        for i in range(N):
            if i == 0:
                ui[:,i] = np.concatenate([u1xy, [[uiz[0]]]])[:,i]
            else:
                Rz = R.from_euler('z', thetas[i], degrees=False).as_matrix()
                ui[:,i] = (Rz @ ui[:,0].reshape(-1,1) + dthetas[i] * np.array([0,0,1]).reshape(3,1)).flatten()

        
        # --------- Compute loads ---------
        ti = np.zeros((3,N))
        fbi = np.zeros((3,N))
        mbi = np.zeros((3,N))
        mb = np.zeros((3,1))
        for i in range(N):
            # Rotation matrix
            Rz = R.from_euler('z', thetas[i], degrees=False).as_matrix()
            _, _, A, _, Kbe = ctr_set[i].getStiffnessAtPosition(s)
            # i-th rod gravity force
            fbi[:,i] = (ctr_set[i].density * A * self.gravity).flatten()
            # i-th rod magnetic torque
            mu = Rq1 @ Rz @ ((ctr_set[i].Br * A / (np.pi * 4e-7)) * ctr_set[i].M)
            ti[:,i] = (self.skew(mu.flatten()) @ self.b_T).flatten() 
            # i-th rod internal moment (eqn. 5d)
            mbi[:,i] = (Kbe @ (ui[:,i].reshape(-1,1) - ctr_set[i].precurvature)).flatten()
            # total moment (eqn. 5e)
            mb = mb + Rz @ mbi[:,i].reshape(-1,1) #- ti[:,i].reshape(-1,1)

        # Sum forces
        fb = np.sum(fbi, axis=1).reshape(-1,1)

        # --------- Arc-length derivatives ---------
        # Position (eqn. 4a)
        dp1 = Rq1 @ np.array([0,0,1]).reshape(-1,1)
        # Quaternion (eqn. 4b) --> rows are switched up compared to paper because here, real part of q1 is last (instead of first)
        Q = np.zeros((4,3))
        Q[3,:] = -q1[:3].T
        Q[0:3:,:] = q1[-1] * np.identity(3) - self.skew(q1[:3].flatten())
        dq1 = (0.5 * Q @ Rq1 @ ui[:,0]).reshape(-1,1)
        # Force (eqn. 4c)
        dn = -fb
        # Moment (eqn. 4d-e)
        dmbxy = ( self.skew(ui[:,0]) @ mb - self.skew([0,0,1]) @ Rq1.T @ (n.reshape(-1,1)) - np.sum(ti, axis=1).reshape(-1,1))[:2,:2]
        dmbz = np.zeros((N,1))
        for i in range(N):
            dmbz[i] = -np.array([0,0,1]).reshape(1,3) @ self.skew(ui[:,i]) @ mbi[:,i] - ti[-1,i]
        # Gradient state vector
        dy = np.concatenate((dp1, dq1, dn, dmbxy, dmbz, dthetas)).flatten()
        
        return dy
        

            




    def skew(self, v):
        return np.array([[    0, -v[2],  v[1]],
                        [ v[2],     0, -v[0]],
                        [-v[1],  v[0],     0]])

