from cProfile import label
from re import T
from this import d
from tkinter import E
from turtle import color
from xml.sax.handler import DTDHandler
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cvxopt
from cvxopt import matrix
from cvxopt import solvers
import time

import osqp
from scipy import sparse

import cvxpy as cp

import csv

# VERSION 1.5
# CHANGE THE DYNAMIC EQUATION WITH HUMAN TORQUE AND ROBOT TORQUE
# MODIFY THE QP SOLVER WITH Control input u1(for acceleration limit) and delt(for slack)
class SafetyControl():
    def __init__(self, max_joint_limit, max_joint_acceleration, min_joint_acceleration, H_acc):
        
        ### Admittance Controller I,B,K SET
        self.IBK = dict() # order : x, y, z, roll, pitch, yaw
        self.IBK['orientation'] = dict()
        self.IBK['orientation']['I'] = np.diag(np.array([0.01, 0.01, 0.01]))
        #self.IBK['orientation']['B'] = np.diag(np.array([0.002, 0.002, 0.002]))
        self.IBK['orientation']['K'] = np.diag(np.array([0.0, 0.0, 0.0]))
        I = self.IBK['orientation']['I']
        #B = self.IBK['orientation']['B']
        K = self.IBK['orientation']['K']
        
        # Human Damping Now Constant but after experiment it should be equation
        #self.B_H = np.array([0, 0, 0])
        
        ### Control Affine Sys Matrix
        ### Constraint Parameter
        self.max_joint_limit = max_joint_limit
        self.min_joint_acceleration = min_joint_acceleration
        self.max_joint_acceleration = max_joint_acceleration
        self.max_b_r = 100.0 # for set min max input robot damping b_r

        # for change FT sensor direction 
        self.FTdirectionFlag = True
        self.wtConflictValue = True   # Turn on the conflict values 
        self.ConstantForceFlag = True # Constant Force
        self.variable_damping_flag = False # Using variable damping ??
        self.UsingOSQP_flag = True # using OSQP Solver
        self.UsingClosedQP_flag = False # Using Closed form QP.
        
        
        
        self.b_r = 0.002
        self.b_h = 0.0
        if(self.UsingOSQP_flag == True):
            self.P = sparse.csc_matrix([[H_acc, 0, 0], [0, H_acc, 0], [0, 0, H_acc]])
            self.q_matrix = np.zeros((3))
        else:
            self.P_matrix = np.diag(np.array([H_acc, H_acc, H_acc]))
            self.q_matrix = np.zeros((3))
        
        self.b_r_x = 0
        self.b_r_y = 0
        self.b_r_z = 0
        
        
    # This function have impedance control eq
    def calControlAffineSysMatrix(self,x,I,B,K,x_eq, force_applied): 
        self.force_applied = force_applied   # from F/T sensor force translated to torques on human shoulder
        self.x = x
        self.position = np.array([x[0],x[1],x[2]])   # 3 by 1 array
        self.velocity = np.array([x[3],x[4],x[5]])   # 3 by 1 array

        A_temp = np.dot(np.linalg.pinv(-I), K) 
        self.A_temp = A_temp
        a = np.zeros((3,3))
        b = np.diag(np.array([1,1,1]))
        a_temp = np.column_stack((a,b))
        b_temp = np.column_stack((A_temp,a))
        self.F = np.vstack((a_temp,b_temp))    # 6 by 6
        
        # Gamma need to change this 'Lambda'
        Lambda_temp = np.matmul(np.linalg.pinv(-I) ,self.velocity)  # 3 by 1
        B_mat_temp = np.matmul(np.linalg.pinv(-I), self.velocity) # 3 by 1
        #B_mat_temp = np.matmul(np.linalg.pinv(-I), self.velocity) # 3 by 1
        # this for the b_h part lambda*b_h in tarun code F
        
        self.Lambda = np.array([[0],
                                [0],
                                [0],
                                Lambda_temp[0],
                                Lambda_temp[1],
                                Lambda_temp[2]])   # 6 by 1
        
        # this for the input lambda*b_r  in tarun code B
        self.B_mat = np.array([[0],
                              [0],
                              [0],
                              B_mat_temp[0],
                              B_mat_temp[1],
                              B_mat_temp[2]]) # 6 by 1
        
        # TO DO LIST about     af=[0;Id\(Kd*Th_r)-Id\(P*b_r)]; this term 
        # E matrix  6 by 1   in tarun code af matrix
        #################################### TO DO LIST - np.dot(B,self.velocity)!!! CHECK
        #e_temp = np.dot(np.linalg.pinv(I), (np.dot(K,x_eq) + self.force_applied) - np.dot(B,self.velocity))
        e_temp = np.dot(np.linalg.pinv(I), (np.dot(K,x_eq) + self.force_applied)) # TODO : Check this part sign -> 03.23.2022
        boldmat = np.diag(B).reshape((3,1))
        f_temp = np.dot(np.linalg.pinv(I),np.array([self.velocity[0]*boldmat[0],self.velocity[1]*boldmat[1],self.velocity[2]*boldmat[2]]))
        
        self.E_mat = np.array([[0],
                              [0],
                              [0],
                              e_temp[0],
                              e_temp[1],
                              e_temp[2]])
        """
        self.E_mat = np.array([[0],
                              [0],
                              [0],
                              e_temp[0]-f_temp[0],
                              e_temp[1]-f_temp[1],
                              e_temp[2]-f_temp[2]])
        """
    
    def MakeCbf(self):
        # this function make reciprocal control barrier function h(x)>=0
        #h=(y_m-sign(x(2))*x(1))-((1/2)*(x(2)^2)/(a_m));
        
        print("self.max_joint_acceleration",self.max_joint_acceleration)
        h_x = (self.max_joint_limit - np.dot(np.sign(self.velocity[0]), self.position[0])) - 1/2*(self.velocity[0]**2/self.max_joint_acceleration)
        h_y = (self.max_joint_limit - np.dot(np.sign(self.velocity[1]), self.position[1])) - 1/2*(self.velocity[1]**2/self.max_joint_acceleration)
        h_z = (self.max_joint_limit - np.dot(np.sign(self.velocity[2]), self.position[2])) - 1/2*(self.velocity[2]**2/self.max_joint_acceleration)
        
        #dh_dx_x2 = -self.position[0] + self.position[0]*np.tanh(self.velocity[0])**2 - self.velocity[0]/self.max_joint_acceleration
        
        
        dh_dx_x = np.array([-np.sign(self.velocity[0]),
                            -self.velocity[0]/self.max_joint_acceleration]) 
        dh_dx_y = np.array([-np.sign(self.velocity[1]),
                            -self.velocity[1]/self.max_joint_acceleration])
        dh_dx_z = np.array([-np.sign(self.velocity[2]),
                            -self.velocity[2]/self.max_joint_acceleration])
        self.h_x = h_x
        self.h_y = h_y
        self.h_z = h_z
        
        self.dh_dx_x = dh_dx_x        # 2 by 1
        self.dh_dx_y = dh_dx_y        # 2 by 1
        self.dh_dx_z = dh_dx_z        # 2 by 1
        
        total_dh_dx = np.hstack((np.hstack((np.transpose(self.dh_dx_x), np.transpose(self.dh_dx_y))), np.transpose(self.dh_dx_z)))
        self.total_dh_dx = total_dh_dx  ## 1 by 6 array
        
        ## CBF Constraints
        # L_fh = dh_dx * (A*x); L_Fh = dh_dx * F; L_gh = dh_dx * B;  L_haf = dh_dx * af;

        A_x_mat = np.array([[0, 1],
                            [self.A_temp[0][0], 0]])   # 2by2
        A_y_mat = np.array([[0, 1],
                    [self.A_temp[1][1], 0]])           # 2by2
        A_z_mat = np.array([[0, 1],
                    [self.A_temp[2][2], 0]])           # 2by2
        
        # L_fh = dh_dx * (A*x);
        L_fh_x = np.matmul(self.dh_dx_x.T, np.matmul(A_x_mat, np.array([self.position[0],
                                                                          self.velocity[0]]))) # 1 by 1
        L_fh_y = np.matmul(self.dh_dx_y.T, np.matmul(A_y_mat, np.array([self.position[1],
                                                                          self.velocity[1]]))) # 1 by 1
        L_fh_z = np.matmul(self.dh_dx_z.T, np.matmul(A_z_mat, np.array([self.position[2],
                                                                          self.velocity[2]]))) # 1 by 1
        
        
        
        ## This term need * b_h?
        L_Fhh_x = np.matmul(self.dh_dx_x.T, np.array([self.Lambda[0],self.Lambda[3]]))  # 1 by 1
        L_Fhh_y = np.matmul(self.dh_dx_y.T, np.array([self.Lambda[1],self.Lambda[4]]))  # 1 by 1
        L_Fhh_z = np.matmul(self.dh_dx_z.T, np.array([self.Lambda[2],self.Lambda[5]]))  # 1 by 1
        
        
        # dh/dx*f(x) = L_fh_x + L_Fhh_x + L_haf_x
        # dh/dx*g(x) = L_gh_x
        
        ### This term g(x) vector lie derivative 
        L_gh_x = np.matmul(self.dh_dx_x.T, np.array([self.B_mat[0],self.B_mat[3]])) # 1 by 1
        L_gh_y = np.matmul(self.dh_dx_y.T, np.array([self.B_mat[1],self.B_mat[4]])) # 1 by 1
        L_gh_z = np.matmul(self.dh_dx_z.T, np.array([self.B_mat[2],self.B_mat[5]])) # 1 by 1
        
        L_haf_x = np.matmul(self.dh_dx_x.T, np.array([self.E_mat[0],self.E_mat[3]])) # 1 by 1
        L_haf_y = np.matmul(self.dh_dx_y.T, np.array([self.E_mat[1],self.E_mat[4]])) # 1 by 1
        L_haf_z = np.matmul(self.dh_dx_z.T, np.array([self.E_mat[2],self.E_mat[5]])) # 1 by 1
        
        L_ghu1_x = np.matmul(self.dh_dx_x.T, np.array([0, 1]).T) # 1 by 1
        L_ghu1_y = np.matmul(self.dh_dx_y.T, np.array([0, 1]).T) # 1 by 1
        L_ghu1_z = np.matmul(self.dh_dx_z.T, np.array([0, 1]).T) # 1 by 1
        
        self.L_fh = np.array([[L_fh_x.item(0)], [L_fh_y.item(0)], [L_fh_z.item(0)]]) # 3 by 1          dh_dx * (A*x);
        self.L_Fhh = np.array([[L_Fhh_x.item(0)], [L_Fhh_y.item(0)], [L_Fhh_z.item(0)]])   # 3 by 1    dh_dx * F
        self.L_gh = np.array([[L_gh_x.item(0)], [L_gh_y.item(0)], [L_gh_z.item(0)]]) # 3 by 1          dh_dx * B
        self.L_haf = np.array([[L_haf_x.item(0)], [L_haf_y.item(0)], [L_haf_z.item(0)]]) # 3 by 1      dh_dx * af
        self.L_ghu1 = np.array([[L_ghu1_x.item(0)], [L_ghu1_y.item(0)], [L_ghu1_z.item(0)]])
        
        
        
        ## For FC Constarints
        self.f_x_x_temp = np.dot(A_x_mat, np.array([self.position[0],self.velocity[0]]))
        self.f_x_y_temp = np.dot(A_y_mat, np.array([self.position[1],self.velocity[1]]))
        self.f_x_z_temp = np.dot(A_z_mat, np.array([self.position[2],self.velocity[2]]))
        self.lambda_x_temp = np.array([self.Lambda[0],self.Lambda[3]])
        self.lambda_y_temp = np.array([self.Lambda[1],self.Lambda[4]])
        self.lambda_z_temp = np.array([self.Lambda[2],self.Lambda[5]])
        self.g_x_temp = np.array([self.B_mat[0],self.B_mat[3]])
        self.g_y_temp = np.array([self.B_mat[1],self.B_mat[4]])
        self.g_z_temp = np.array([self.B_mat[2],self.B_mat[5]])
        self.e_x_temp = np.array([self.E_mat[0],self.E_mat[3]])
        self.e_y_temp = np.array([self.E_mat[1],self.E_mat[4]])
        self.e_z_temp = np.array([self.E_mat[2],self.E_mat[5]])
        
        
    def InitQP(self):
        pass
    
    def QPMatrixFormulization(self,position,velocity,h_x,L_fh_x,L_Fhh_x,L_gh_x,L_haf_x,L_ghb_x,gamma_safe,B_old):
        # Make Matrix for QP Objective and Constraints
        # AA*x <= bb
        # Decision variable vector : 3 by 1 [b_r, delta, u1]
        # Constraint order 1)delta>0, 2)control input constraint, 3)Acceleration contraints, 4)CBF constraints. 
        II = np.array([0,1,0,-1]).reshape(2,2)
        ac_AA_input = np.dot(II,self.Bmat)
        ac_AA_u1 = np.dot(II,self.BB)
        
        AA = np.array([0,-1,0, 
                       1,-1,0, 
                       -1,-1,0, 
                       ac_AA_input[0,0],0,ac_AA_u1[0,0], 
                       ac_AA_input[1,0],0,ac_AA_u1[1,0], 
                       -L_gh_x[0,0],0,-L_ghb_x[0,0]]).reshape(6,3)
        temp_state = np.array([position[0], velocity[0]]).reshape(2,1)
        ac_bb1 = np.dot(self.Amat,temp_state) + self.Fmat*self.b_h + self.AFmat  # 2 by 1
        ac_bb2 = -np.dot(II,ac_bb1) + np.array([self.max_joint_acceleration,self.max_joint_acceleration]).reshape(2,1) # 2 by 1
        
        cbf_bb = L_fh_x+L_Fhh_x*self.b_h+L_haf_x+gamma_safe*h_x
        
        bb = np.array([0,self.max_b_r,self.max_b_r,ac_bb2[0,0],ac_bb2[1,0],cbf_bb[0,0]]).reshape(6,1)
        
        HH = np.diag(np.array([1,1,1]))
        FF = np.array([-B_old[0,0],0,0]).reshape(3,1)
        print("AA",AA)
        print("bb",bb)
        print("HH",HH)
        print("FF",FF)
        
        print("AA",AA.shape)
        print("bb",bb.shape)
        print("HH",HH.shape)
        print("FF",FF.shape)
        
        
        return AA, bb, HH, FF
    
    def QpSolver(self,B_old,x_t_s, v_t_s, a_t_s, gamma_safe = 10, epsilon = 10):
        
        #this term for making constraint   h_dot >=0   ->   dh/dx*[f(x) + g(x)*b_r] >= -alpha*h   ->   -dh/dx*g(x)*b_r =< alpha*h + dh/dx*f(x)
        # h_temp = dh/dx * f(x) + alpha*h = L_fh_x + L_Fhh_x + L_haf_x + alpha*h
        h_temp_x = self.L_fh[0] + self.L_Fhh[0]*self.b_h + self.L_haf[0] + gamma_safe*self.h_x    ## dh/dx * f(x) + alpha*h
        h_temp_y = self.L_fh[1] + self.L_Fhh[1]*self.b_h + self.L_haf[1] + gamma_safe*self.h_y
        h_temp_z = self.L_fh[2] + self.L_Fhh[2]*self.b_h + self.L_haf[2] + gamma_safe*self.h_z
        
        gamma_mat = np.array([[0.0, 1.0],
                              [0.0,-1.0]])
        #just f(x)
        x_dot_x = self.f_x_x_temp + self.lambda_x_temp*self.b_h + self.e_x_temp
        
        fc_x = np.dot(gamma_mat,x_dot_x)
        fc_g = np.dot(gamma_mat,self.g_x_temp)
        fc_g1 = np.dot(gamma_mat, np.array([0,1]).T)
        
        # just f(x)
        f_temp_x = self.f_x_x_temp[1,0] + self.lambda_x_temp[1,0]*self.b_h + self.e_x_temp[1,0]
        f_temp_y = self.f_x_y_temp[1,0] + self.lambda_y_temp[1,0]*self.b_h + self.e_y_temp[1,0]
        f_temp_z = self.f_x_z_temp[1,0] + self.lambda_z_temp[1,0]*self.b_h + self.e_z_temp[1,0]
        # (B_old[0][0]/self.max_joint_acceleration)*self.velocity[0]**2 = affine term
        h_dot_x = self.L_fh[0] + self.L_Fhh[0]*self.b_h + self.L_haf[0] + self.L_gh.item(0)*B_old[0][0] +self.L_ghu1[0]
        h_dot_x = h_dot_x[0]
        #direct_h_dot =  -a_t_s[0]*x_t_s[0] + np.tanh(v_t_s[0])**2*x_t_s[0]*a_t_s[0] - np.tanh(v_t_s[0])*v_t_s[0] - (v_t_s[0]*a_t_s[0])/self.max_joint_acceleration
        
        direct_h_dot =  -np.sign(v_t_s[0,0])*v_t_s[0,0] - (v_t_s[0,0]*a_t_s[0,0])/self.max_joint_acceleration
        print("direct_h_dot",direct_h_dot)
        
        print("a_t_s[0,0]",a_t_s[0,0])
        print("v_t_s[0,0]",v_t_s[0,0])
        print("direct_h_dot",direct_h_dot)
        
        kappa_h_dot = -gamma_safe*self.h_x
        
        print("kappa_h_dot",kappa_h_dot[0])
        
        self.kappa_h_dot = kappa_h_dot[0]
        
        # Quadratic Program
        if(self.UsingClosedQP_flag == True):
            """
            P = np.diag(np.array([1.0, 1.0, 1.0]))
            q = np.array([[-B_old[0][0]],[-B_old[1][1]],[-B_old[2][2]]])
            x_d = np.array([[B_old[0][0]],[B_old[1][1]],[B_old[2][2]]])
            a = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0],[-self.L_gh.item(0), 0.0, 0.0], [0.0, -self.L_gh.item(1), 0.0], [0.0, 0.0, -self.L_gh.item(2)],[fc_g[0,0], 0.0, 0.0],[fc_g[1,0], 0.0, 0.0]]).T
            b = np.array([self.max_b_r,self.max_b_r,self.max_b_r,self.max_b_r,self.max_b_r,self.max_b_r,h_temp_x[0],h_temp_y[0],h_temp_z[0],-fc_x[0,0] + self.max_joint_acceleration,-fc_x[1,0] + self.max_joint_acceleration]).reshape((11,1))
            """
            if(direct_h_dot < self.kappa_h_dot):
                # From https://dev10110.github.io/tech-notes/maths/qp_closed_form.html
                # Our Objective Function J = 1/2*(b_new -b_old)^2 +
                P = np.diag(np.array([0.5, 0.5, 1000]))
                #P = P.T @ P
                q = np.array([-B_old[0][0], 10, 0]).reshape(3,1) 
                
                a = np.array([0,-1,0 ,1,-1,0, -1,-1,0, -self.L_gh.item(0),0,0, fc_g[0,0],0,0, fc_g[1,0],0,0]).reshape((6,3)).T   # a=3 by 6
                b = np.array([0, self.max_b_r,self.max_b_r, h_temp_x, -fc_x[0,0] + self.max_joint_acceleration, -fc_x[1,0] + self.max_joint_acceleration]).reshape((6,1))
                #h_dot_x
                
                
                    
                s1 = -np.linalg.multi_dot([a.T, np.linalg.pinv(P), a])  ## S = -a^T*P^-1*a
                a1 = -np.dot(np.linalg.pinv(P),q)
                a2 = -np.linalg.multi_dot([np.linalg.inv(P), a, np.linalg.inv(s1), a.T, np.linalg.inv(P), q])
                a3 = -np.linalg.multi_dot([np.linalg.inv(P), a, np.linalg.inv(s1), b])
                sol_value = a1+a2+a3
                
                                    
                #sol_value = 
                OpDampingFlag = True
            else:
                # From https://dev10110.github.io/tech-notes/maths/qp_closed_form.html
                P = np.diag(np.array([1, 1, 1000]))
                P = P.T @ P
                q = np.array([-B_old[0][0], 10, 0]).reshape(3,1) 
                
                a = np.array([0,-1,0 ,1,-1,0, -1,-1,0, fc_g[0,0],0,0, fc_g[1,0],0,0]).reshape((5,3)).T   # a=3 by 6
                b = np.array([0, self.max_b_r,self.max_b_r, -fc_x[0,0] + self.max_joint_acceleration, -fc_x[1,0] + self.max_joint_acceleration]).reshape((5,1))
                h_dot_x
                
                
                    
                s1 = -np.linalg.multi_dot([a.T, np.linalg.pinv(P), a])  ## S = -a^T*P^-1*a
                a1 = -np.dot(np.linalg.pinv(P),q)
                a2 = -np.linalg.multi_dot([np.linalg.pinv(P), a, np.linalg.pinv(s1), a.T, np.linalg.pinv(P), q])
                a3 = -np.linalg.multi_dot([np.linalg.pinv(P), a, np.linalg.pinv(s1), b])
                sol_value = a1+a2+a3
                OpDampingFlag = True
            #print("P",P.shape)
            #print("q",q.shape)
            #print("x_d",x_d.shape)
            #print("a",a.shape)
            #print("b",b.shape)
            #print("sol_value",sol_value)
            damping_status = "Closed"
        else:
            ## Modify We have 3 decision variable
            x = cp.Variable((3,1))
            P = np.diag(np.array([1, 1, 1000]))
            P = P.T @ P
            q = np.array([-B_old[0][0], 10, 0]).reshape(3,1)   ### 마이너스 안붙였었음
            #G = np.array([1, -1, -self.L_gh.item(0), fc_g[0,0], fc_g[1,0]])
            G = np.array([0,-1,0 ,1,-1,0, -1,-1,0, -self.L_gh.item(0),0,-self.L_ghu1.item(0), fc_g[0, 0],0,fc_g1[0], fc_g[1,0],0,fc_g1[1]]).reshape((6,3))
            h = np.array([0, self.max_b_r,self.max_b_r, h_temp_x, -fc_x[0,0] + self.max_joint_acceleration, -fc_x[1,0] + self.max_joint_acceleration]).reshape((6,1))
            #G = np.array([1, -1, -self.L_gh.item(0)])
            #h = np.array([self.max_b_r, self.max_b_r, h_temp_x])
            print("P size", P.shape)
            print("x size", x.shape)
            print("q size", q.T.shape)
            print("G size", G.shape)
            print("h size", h.shape)
            
            objective = cp.Minimize((1/2)*cp.quad_form(x,P) + q.T @ x)  ## 1/2 안붙였었음
            constraints = [G @ x <= h]
            
            prob = cp.Problem(objective, constraints)
            

            prob.solve(solver='OSQP',verbose =True)
            sol_value = x.value
        
            if(prob.status == "optimal"):
                OpDampingFlag = True
            else:
                OpDampingFlag = False
            
            damping_status = "Libarary"
        
        #print("prob.status",prob.status)
        #print("sol_value",sol_value)
        
        print("OpDampingFlag", OpDampingFlag)
    
        return sol_value, OpDampingFlag, damping_status, h_dot_x, direct_h_dot
    
    def ClosedQP(self):
        
        
                
        pass
    
    def setStateInit(self, x_init,y_init,z_init,x_dot_init,y_dot_init,z_dot_init,x_ddot_init,y_ddot_init,z_ddot_init): # inpout : in Degree (theta = pitch, phi = yaw)
        
        self.theta_0 = np.array([[x_init, y_init, z_init]]).T
        self.theta_00 = np.array([[x_init, y_init, z_init]]).T
        self.theta_000 = np.array([[x_init, y_init, z_init]]).T
        self.theta_new = np.array([[x_init, y_init, z_init]]).T
        
        self.theta_dot_0 = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_00 = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_000 = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_new = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        
        self.theta_ddot_0 = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_00 = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_000 = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_new = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        
    def setSafetyStateInit(self, x_init,y_init,z_init,x_dot_init,y_dot_init,z_dot_init,x_ddot_init,y_ddot_init,z_ddot_init): # inpout : in Degree (theta = pitch, phi = yaw)
        
        self.theta_0_s = np.array([[x_init, y_init, z_init]]).T
        self.theta_00_s = np.array([[x_init, y_init, z_init]]).T
        self.theta_000_s = np.array([[x_init, y_init, z_init]]).T
        self.theta_new_s = np.array([[x_init, y_init, z_init]]).T
        
        self.theta_dot_0_s = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_00_s = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_000_s = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        self.theta_dot_new_s = np.array([[x_dot_init, y_dot_init, z_dot_init]]).T
        
        self.theta_ddot_0_s = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_00_s = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_000_s = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
        self.theta_ddot_new_s = np.array([[x_ddot_init, y_ddot_init, z_ddot_init]]).T
    
    def getHumanIntent(self,v_t,a_t,v_t_s,a_t_s):
        # human intent = velocity * acceleration
        
        #WO SAFETY
        humanintent_x = v_t[0] * a_t[0]
        humanintent_y = v_t[1] * a_t[1]
        humanintent_z = v_t[2] * a_t[2]
        humanintent = np.array([humanintent_x, humanintent_y, humanintent_z])
        #W SAFETY
        humanintent_s_x = v_t_s[0] * a_t_s[0]
        humanintent_s_y = v_t_s[1] * a_t_s[1]
        humanintent_s_z = v_t_s[2] * a_t_s[2]
        humanintent_s = np.array([humanintent_s_x, humanintent_s_y, humanintent_s_z])
        
        return humanintent, humanintent_s
    
    def VariableDamping(self, b_c, b_low,b_upper,s, humanintent):
        xx_max = 49250
        xx_min = -9859
        k_p = -np.log((1-s)/(1+s)) / xx_max
        k_n = -np.log((1+s)/(1-s)) / xx_min
        
        if(humanintent[0]>xx_max):
            humanintent[0][0]=xx_max
        elif(humanintent[0]<xx_min):
            humanintent[0][0] = xx_min
        else:
            pass
        
        if (humanintent[0][0]>=0):
            self.b_r_x = (2*b_low/(1+np.exp(-k_p*humanintent[0][0]))) -b_low + b_c
        elif(humanintent[0][0]<0):
            self.b_r_x = -(2*b_upper/(1+np.exp(-k_n*humanintent[0][0]))) +b_upper + b_c
        if(humanintent[1][0]>=0):
            self.b_r_y = (2*b_low/(1+np.exp(-k_p*humanintent[1][0]))) -b_low + b_c
        elif(humanintent[1][0]<0):
            self.b_r_y = -(2*b_upper/(1+np.exp(-k_n*humanintent[1][0]))) +b_upper + b_c
        if(humanintent[2][0]>=0):
            self.b_r_z = (2*b_low/(1+np.exp(-k_p*humanintent[2][0]))) -b_low + b_c
        elif(humanintent[2][0]<0):
            self.b_r_z = -(2*b_upper/(1+np.exp(-k_n*humanintent[2][0]))) +b_upper + b_c
        else:
            pass
        b_r_x_save = self.b_r_x
        ## After using human intention
        B_ah = np.diag(np.array([self.b_h + self.b_r_x, self.b_h + self.b_r_y, self.b_h + self.b_r_z]))   # INIT ROBOT DAMPING
        
        
        return B_ah, b_r_x_save
    
    
    def Doadmittance(self, time_step, I, B, K, x_eq, x_t_1, x_t_2, force_applied):
        
        if(self.FTdirectionFlag == True):
            force_applied = - force_applied
        
        
        t = time_step
        A_dyn =  (I / (t * t) + B / (t) + K) # temporary matrix only used for inverting dynamics behavior in discrete time 3 by 3
        B_temp = np.matmul(I, ((2.0 * x_t_1 - x_t_2) / (t * t))) #3 by 1
        C_temp = np.matmul(B, (x_t_1 / (t))) #3 by 1
        D_temp = np.matmul(K, x_eq) #3 by 1
        
        x_t = np.matmul(np.linalg.pinv(A_dyn),  (-force_applied + B_temp + C_temp + D_temp))
        
        v_t = (x_t-x_t_1)/(t)
        
        a_t = (v_t - self.theta_dot_0)/(t)
        
        self.theta_new = x_t
        self.theta_dot_new = v_t
        self.theta_ddot_new = a_t
        
        self.theta_000 = self.theta_00
        self.theta_00 = self.theta_0
        self.theta_0 = self.theta_new
        
        self.theta_dot_000 = self.theta_dot_00
        self.theta_dot_00 = self.theta_dot_0 
        self.theta_dot_0 = self.theta_dot_new
        
        self.theta_ddot_000 = self.theta_ddot_00
        self.theta_ddot_00 = self.theta_ddot_0 
        self.theta_ddot_0 = self.theta_ddot_new
        
        return x_t, v_t, a_t
    
    def SafeDoadmittance(self, time_step,sol_value, I, B, K, x_eq, x_t_1, x_t_2, force_applied, status):
        
        if(self.FTdirectionFlag == True):
            force_applied = - force_applied
        
        u1 = sol_value[2,0]
        t = time_step
        print("B",B)
        A_dyn =  (I / (t * t) + B / (t) + K) # temporary matrix only used for inverting dynamics behavior in discrete time 3 by 3
        print("here")
        B_temp = np.matmul(I, ((2.0 * x_t_1 - x_t_2) / (t * t))) #3 by 1
        C_temp = np.matmul(B, (x_t_1 / (t))) #3 by 1
        D_temp = np.matmul(K, x_eq) #3 by 1
        
        x_t_s = np.matmul(np.linalg.pinv(A_dyn),  (-force_applied + u1 + B_temp + C_temp + D_temp))
        
        v_t_s = (x_t_s-x_t_1)/(t)
        
        a_t_s = (v_t_s - self.theta_dot_0_s)/(t)#(x_t_s-2*x_t_1+x_t_2)/t**2
        
        
        if(status == "UsingSafetyValue"):
            self.theta_new_s = x_t_s
            self.theta_dot_new_s = v_t_s
            self.theta_ddot_new_s = a_t_s
            
            self.theta_000_s = self.theta_00_s
            self.theta_00_s = self.theta_0_s
            self.theta_0_s = self.theta_new_s
            
            self.theta_dot_000_s = self.theta_dot_00_s
            self.theta_dot_00_s = self.theta_dot_0_s 
            self.theta_dot_0_s = self.theta_dot_new_s
            
            self.theta_ddot_000_s = self.theta_ddot_00_s
            self.theta_ddot_00_s = self.theta_ddot_0_s
            self.theta_ddot_0_s = self.theta_ddot_new_s
        else:
            #print("Not Save")
            pass
        return x_t_s, v_t_s, a_t_s
    
    def torque_trajectory(self, a, n,index, torque_type):
        
        if(torque_type == "constant"):
            torque_applied = a
        if(torque_type == "sine"):
            torque_applied = a*np.sin(index*(2*np.pi/n))
        if(torque_type == "minjerk"):
           torque_applied = a*(10*(index/n)**3 - 15*(index/n)**4 + 6*(index/n)**5)
        torques_applied = np.array([[torque_applied, torque_applied, torque_applied]]).T
        if(index == n):
            index = 0
        return torques_applied, index
    
    def DoSafetyController(self, cycle_time,sol_value,I, B_S, K, x_eq, forces_applied, status):
        # 여기서는 state를 저장 안함
        
        x_t_s, v_t_s, a_t_s = self.SafeDoadmittance(cycle_time,sol_value, I, B_S, K, x_eq, self.theta_0_s,
                                                                self.theta_00_s, 
                                                                forces_applied,
                                                                status)
        # 6 by 1 state vector with x(t)
        x = np.array([x_t_s[0],x_t_s[1],x_t_s[2],v_t_s[0],v_t_s[1],v_t_s[2]])  
        acc_save = a_t_s[0]
        # with b_old and state made by b_old control affine construct with x(t)
        self.calControlAffineSysMatrix(x,I,B_S,K,x_eq, forces_applied)
        #self.MakeClf()
        self.MakeCbf()
        
        #To Do fix gamma safe term   constant 10 = gammasafe
        #h_dot_x = self.L_fh[0] + self.L_Fhh[0]*self.b_h + self.L_haf[0] + self.L_gh.item(0)  #10*self.h_x 
        
        
        sol_value, OpDampingFlag,damping_status ,h_dot_x, direct_h_dot = self.QpSolver(B_S,x_t_s,v_t_s,a_t_s, gamma_safe = 10, epsilon = 10)
        
        return sol_value, OpDampingFlag,x , v_t_s,a_t_s, damping_status ,h_dot_x, direct_h_dot
    
    def after_h_dot(self,v_t_s,a_t_s):
        
        after_hdot =  -np.sign(v_t_s[0,0])*v_t_s[0,0] - (v_t_s[0,0]*a_t_s[0,0])/self.max_joint_acceleration
        return after_hdot

    
    def Run(self):

        #Init IBK values
        I = self.IBK['orientation']['I']
        #B = self.IBK['orientation']['B']
        B = np.diag(np.array([self.b_h + self.b_r, self.b_h + self.b_r, self.b_h + self.b_r]))   # INIT ROBOT DAMPING
        
        K = self.IBK['orientation']['K']
        x_eq = np.zeros((3,1))
        count = 0.0
        B_S = np.diag(np.array([self.b_h + self.b_r, self.b_h + self.b_r, self.b_h + self.b_r]))  # INIT ROBOT DAMPING FOR SAFETY CONTROLLER
        sol_value = np.zeros((3,1))
        # init val
        self.setStateInit(x_init = 0.0,y_init = 0.0,z_init = 0.0,x_dot_init = 0.0,y_dot_init = 0.0,z_dot_init = 0.0,x_ddot_init = 0.0,y_ddot_init = 0.0,z_ddot_init = 0.0)
        self.setSafetyStateInit(x_init = 0.0,y_init = 0.0,z_init = 0.0,x_dot_init = 0.0,y_dot_init = 0.0,z_dot_init = 0.0,x_ddot_init = 0.0,y_ddot_init = 0.0,z_ddot_init = 0.0)
        cycle_time = 0.0001#0.1ms
        num = 0.0
        
        
        ## Save Without safetycontroller
        save_x_t = []
        save_v_t = []
        save_a_t = []
        save_b_r = []
        
        ## Save With safetycontroller
        save_x_t_s = []
        save_v_t_s = []
        save_a_t_s = []
        save_b_r_s = []
        
        save_humanitent_wo = []
        save_humanitent_w = []
        save_torque = []
        
        save_h_dot_x = []
        save_direct_h_dot = []
        
        save_kappa_h_dot = []
        save_h_x = []
        
        save_after_h_dot_x = []
        
        humanintent = np.zeros((3,1))
        humanintent_s = np.zeros((3,1))
        ## Now not using this
        """
        if(self.UsingOSQP_flag == True):
            prob = self.InitQP()
        else:
            prob = 0
        """
        #try:
        while( num <= 200):
            
            forces_applied, count = self.torque_trajectory(a = 2, n = 200,index = count, torque_type = "sine")
            save_torque.append(forces_applied[0,0])
            save_humanitent_wo.append(humanintent[0,0])
            save_humanitent_w.append(humanintent_s[0,0])
            ## Variable Damping Part
            if(self.variable_damping_flag == True):
                #humanintent = self.getHumanIntent()
                B_W,b_r_x_save = self.VariableDamping(b_c = 0, b_low = -40,b_upper = 50,s=0.95,humanintent=humanintent_s)
                B_WO,b_r_x_save = self.VariableDamping(b_c = 0, b_low = -40,b_upper = 50,s=0.95,humanintent=humanintent)
                print("B_WO",B_WO)
            else:
                B_W = np.diag(np.array([0.1, 0.1, 0.1]))
                B_WO = np.diag(np.array([0.1, 0.1, 0.1]))
                
            
            
            b_x_old = B_WO[0][0]
            
        
        
            x_t, v_t, a_t = self.Doadmittance(cycle_time, I, B_WO, K, x_eq, self.theta_0, self.theta_00, forces_applied)
            
            ## with out safecontroller data save
            save_x_t.append(x_t[0,0])
            save_v_t.append(v_t[0,0])
            save_a_t.append(a_t[0,0])
            save_b_r.append(b_x_old)
            
            
            sol_value, OpDampingFlag, x_old, v_t_s, a_t_s,damping_status, h_dot_x, direct_h_dot = self.DoSafetyController(cycle_time,sol_value,I, B_W, K, x_eq, forces_applied, status = "SafetyController")
            
            save_kappa_h_dot.append(self.kappa_h_dot)
            save_h_x.append(self.h_x[0])
            # compare h_dot value , this value is before h_dot
            save_h_dot_x.append(h_dot_x)
            save_direct_h_dot.append(direct_h_dot)
            
            # If QP solver solve the SAFE robot damping b_r, using this safe value we do admittance onmore time (safety version admittance)
            if(OpDampingFlag == True):
                ## Now we just solve x orientation b sol I put this value on all B_S

                print("sol_value",sol_value)
                print("sol_value",sol_value.shape)
                if(damping_status == "Libarary"):
                    B_S = np.diag(np.array([self.b_h + sol_value[0][0], self.b_h + sol_value[0][0], self.b_h + sol_value[0][0]]))
                else:
                    a = sol_value[0,0]
                    
                    B_S = np.diag(np.array([a, a, a]))
                    
                    #B_S = np.array([[self.b_h + sol_value.item(0),0,0], [0,self.b_h + sol_value.item(1),0], [0,0,self.b_h + sol_value.item(2)]])
                    print("hihihi")
                    print("B_S",B_S.shape)
            else:
                B_S = sol_value
                print("B_S",sol_value)
            x_t_s, v_t_s, a_t_s= self.SafeDoadmittance(cycle_time,sol_value, I, B_S, K, x_eq, x_t_1 = self.theta_0_s,x_t_2 = self.theta_00_s, force_applied = forces_applied,status = "UsingSafetyValue")
            ## safecontroller data save
            save_x_t_s.append(x_t_s[0,0])
            save_v_t_s.append(v_t_s[0,0])
            save_a_t_s.append(a_t_s[0,0])
            save_b_r_s.append(B_S[0,0])
            after_hdot = self.after_h_dot(v_t_s,a_t_s)
            save_after_h_dot_x.append(after_hdot)
            
            humanintent, humanintent_s = self.getHumanIntent(v_t,a_t,v_t_s,a_t_s)  ## get humanintent after robot moving with admittance controller. 
            

            
            count += 1
            
            num += 1
        #plt.show()
        #print("save_br", save_br)
        
        #### Make CSV
        row_title =["save_x_t", "save_v_t", "save_a_t", "save_b_r","save_x_t_s","save_v_t_s","save_a_t_s","save_b_r_s","save_humanitent_wo","save_humanitent_w","save_torque", "save_h_dot_x", "save_direct_h_dot","save_kappa_h_dot","save_h_x","save_after_h_dot_x"]
        rows = zip(save_x_t, save_v_t, save_a_t, save_b_r, save_x_t_s, save_v_t_s, save_a_t_s, save_b_r_s,save_humanitent_wo,save_humanitent_w,save_torque, save_h_dot_x, save_direct_h_dot,save_kappa_h_dot,save_h_x,save_after_h_dot_x)
        newfilePath = "Data.csv"
        with open(newfilePath, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(row_title)
            for row in rows:
                writer.writerow(row)
        
        
    
        if(self.ConstantForceFlag == True):
            plt.savefig('simul result Constant Force', dpi=300)
        else:
            plt.savefig('simul result Sine Force', dpi=300)
        
        #except KeyboardInterrupt:
        #    print("Exception")    
        #    exit()

def main():
    safe_controller = SafetyControl(max_joint_limit = 0.5, max_joint_acceleration = 100, min_joint_acceleration = 1000, H_acc = 0.5)
    
    #try:
    #print("trying to run...")
    safe_controller.Run()
        
    #except KeyboardInterrupt:
    #    exit()
    
    
if __name__ == '__main__':
    #try:
    main()
    #except:
    #    pass
