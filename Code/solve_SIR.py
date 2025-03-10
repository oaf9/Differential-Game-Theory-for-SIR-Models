from scipy.integrate import odeint, solve_ivp
from scipy.interpolate import interp1d
from scipy.integrate import simpson
import numpy as np

class solveSIR:

    #class variables 
    S_hat = None
    I_hat = None
    R_hat = None
    λ = None


    #constructor
    def __init__(self, S_true, I_true, R_true, N, p_0):
    
        self.tf = len(S_true)
        self.N = N
        self.p = p_0
        self.time_steps = np.arange(0, len(S_true), 1)

        self.I = I_true
        self.S = S_true
        self.R = R_true

    #the forwards equations
    def SIR(self, V, t):

        """defines the SIR model dynamics
        V = (S, I, R)
        args = (N,Β, γ)
        """

        B, g  = self.p[0:2]
        S, I, R = V

        dS = -(B/self.N)*S*I
        dI = (B/self.N)*S*I - g*I
        dR = g*I

        return dS, dI, dR
    
    def forwardSolve(self):

        """solves the forward problem
        """
        return odeint(self.SIR, 
                      y0 = self.p[2:],
                      t = self.time_steps)


    #the backwards equtions
    def lamdaDE(self, t, λ):
        """
        solve the backwards problem for λ(t)
        """
        β, γ  = self.p[0:2]

        # simpler names for readability
        t = int(t)
        I, I_hat, S_hat, N = self.I[t], self.I_hat[t], self.S_hat[t], self.N
 
        df_dx = 2*np.array([0, I - I_hat, 0])
        dh_dx = np.array([[β*I_hat/N  ,  β*S_hat/N, 0],
                          [-β*I_hat/N , -β*S_hat/N + γ, 0],
                          [0          ,         -γ, 0]])

        #return the negative since we are moving backwards
        return  -(df_dx + np.array(λ).T@dh_dx)
    
    
    def backwardsSolve(self):
        """Backwards solve for λ(t). 
        """
        #we have to integrate backwards
        #output must be flipped since odeint will sovle a forward problem.
        return solve_ivp(fun = self.lamdaDE, 
                        y0 = [0,0,0], 
                        t_span = [0, self.time_steps[-1]],
                        t_eval = self.time_steps).y.T[::-1]


    #the integral for dF_dp
    def dh_dp(self,t):

        I, I_hat, S_hat, N = self.I[t], self.I_hat[t], self.S_hat[t], self.N

        return np.array([[ I_hat*S_hat/N,      0,0,0,0],
                         [-I_hat*S_hat/N,  I_hat,0,0,0],
                         [ 0            , -I_hat,0,0,0]])

    def integrand(self, t, λ): 

        df_dp = np.zeros((5))
        
        return df_dp + λ[t].T@self.dh_dp(t)
    
    def dF_dp(self):

        dg_dx0 = np.eye(3)
        dh_dx_dt = np.eye(3)
        dg_dp = -np.array([[0,0,1,0,0],
                           [0,0,0,1,0],
                           [0,0,0,0,1]])
 
        #we get all the timesteps: 
        integrands = [self.integrand(t, self.λ) for t in self.time_steps]

        #integrate with simpson's method and add the gx(0)gp term 
        return (simpson(y = integrands, x = self.time_steps, axis = 0) + 
                self.λ[0].T@dh_dx_dt@ np.linalg.inv(dg_dx0)@dg_dp)
    

    def fit(self, max_iter = 100, η = .01, ε = .01):
        """A gradient descent implementation to find p*
        """
        i = 0
        while(True):

            print(self.p)

            #step 1 is to integrate for (S,I,R) and update the SIR values
            V = self.forwardSolve()

            self.S_hat = V[:, 0]
            self.I_hat = V[:, 1]
            self.R_hat = V[:, 2]

            #step 2 is to solve the backwards problems 
            self.λ = self.backwardsSolve()

            print("lambda ", self.λ)

            #compute the gradient
            dL_dp = np.clip(self.dF_dp(),-1e3, 1e3)

            #perform the gradient update
            self.p = self.p - η*dL_dp

            #check for convergence: 
            if( (np.sqrt(np.dot(self.p, self.p)) < ε)  or i > max_iter):
                break;
            else:
                i += 1
    
        return self.p, self.forwardSolve()