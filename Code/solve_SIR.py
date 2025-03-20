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

        self.S =  interp1d(np.arange(len(S_true)), S_true, 
                           kind = 'quadratic', 
                           fill_value = "extrapolate", axis=-1)
        
        self.I =  interp1d(np.arange(len(I_true)), I_true, 
                           kind = 'quadratic', 
                           fill_value = "extrapolate", axis=-1)

        self.R =  interp1d(np.arange(len(R_true)), R_true, 
                           kind = 'quadratic', 
                           fill_value="extrapolate", axis=-1)

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
    
    def bSIR(self, V, t):
        return - self.SIR(V,t)

    def bSIRSolve(self, y0):
        """solves the forward problem
        """
        return odeint(self.SIR, 
                      y0 = y0,
                      t = self.time_steps[::-1])[::-1]

    #the backwards eqautions
    def lamdaDE(self, λ, t):
        """
        solve the backwards problem for λ(t)
        """
        β, γ  = self.p[0:2]

        # simpler names for readability

        I, I_hat, S_hat, N = self.I(t), self.I_hat(t), self.S_hat(t), self.N
 
        df_dx = -2*np.array([0, I - I_hat, 0])
        dh_dx = np.array([[β*I_hat/N  ,      β*S_hat/N, 0],
                          [-β*I_hat/N , -β*S_hat/N + γ, 0],
                          [0          ,             -γ, 0]])

        #return the negative since we are integrating backwards
        return  -(df_dx + np.array(λ).T@dh_dx)


    def backwardsSolve(self):
        """Backwards solve for λ(t). 
        """
        #output must be flipped ...the solver will only solve a forward problem.
        return solve_ivp(self.lamdaDE, 
                         y0 =  [0,0,0],
                         t_span = [self.time_steps[-1], self.time_steps[0]],
                         t_eval = self.time_steps[::-1],
                         method='RK45')
    
    def backwardsSolve2(self):
        """Backwards solve for λ(t) using solve_ivp."""
        def lambdaDE_wrapper(t, λ):
            return self.lamdaDE(λ, t)

        # Use solve_ivp for backward integration
        sol = solve_ivp(
            fun=lambdaDE_wrapper,
            t_span=[self.time_steps[-1], self.time_steps[0]],  # Backward in time
            y0=[0, 0, 0],
            t_eval=self.time_steps[::-1],  # Ensure time points match
            method='RK45',  # Adaptive step-size method
        )

        return sol.y.T[::-1]
    
    def forwardLambdaSolve(self, y0):
        """Solve the adjoint λ(t) forward in time using solve_ivp."""
        def lambdaDE_wrapper(t, λ):
            return self.lamdaDE(λ, t)

        sol = solve_ivp(
            fun=lambdaDE_wrapper,
            t_span=[self.time_steps[0], self.time_steps[-1]],
            y0=y0,
            t_eval=self.time_steps,
            method='RK45',
        )

        return sol.y.T



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
        integrands = np.array([self.integrand(t, self.λ) for t in self.time_steps])

        #integrate with simpson's method and add the gx(0)gp term 
        return (simpson(y = integrands, x = self.time_steps, axis = 0) + 
                self.λ[0].T@dh_dx_dt@ np.linalg.inv(dg_dx0)@dg_dp)
    

    def fit(self, max_iter = 200, η = .0001, ε = .01):
        """A gradient descent implementation to find p*
        """
        i = 1
        while(True):

            print(f"round: {i}")

            #step 1 is to integrate for (S,I,R) and update the SIR values
            V = self.forwardSolve()

            self.S_hat =  interp1d(np.arange(len(V)), V[:, 0], 
                    kind = 'quadratic', fill_value="extrapolate", axis=-1)
            self.I_hat =  interp1d(np.arange(len(V)), V[:, 1], 
                    kind = 'quadratic', fill_value="extrapolate", axis=-1)
            self.R_hat =  interp1d(np.arange(len(V)), V[:, 2], 
                    kind = 'quadratic', fill_value="extrapolate", axis=-1)
        
            #step 2 is to solve the backwards problems 
            self.λ = self.backwardsSolve()

            #compute the gradient
            dL_dp = np.clip(self.dF_dp(),-100, 100)
            #dL_dp = self.dF_dp()
            #perform the gradient update
            self.p = self.p - (η/(np.sqrt(i)))*dL_dp


            #check for convergence: 
            if( (np.sqrt(np.dot(dL_dp,  dL_dp)) < ε)  or i > max_iter):
                break;
            else:
                i += 1

    
        return self.p, self.forwardSolve()
    

    #some functions that are usefull for debugging: 
    # this checks that the forward solvers are actually solving the equations 
    # correctly by checking the finite
    def checkForwardEquations(self):

        #It should be the case that 
        #x[i+1] - x[i] = dx
        X = self.forwardSolve()

        S, I, R =  X[:, 0], X[:, 1], X[:, 2]

        B, g  = self.p[0:2]
        N = self.N
        dt = self.time_steps[1]-self.time_steps[0]

        dX = np.array([-B*I*S/N, B*I*S/N - g*I, g*I]).T

        dF = (X[1:,] - X[0:-1,])/dt

        #print(dX[0]) 
        diffs = np.abs((dX[1:,] - dF)/dF)

        return np.mean(diffs)
    
    def checkBackwardEquations(self):
        #we do not update any class variables for the check function.

        #solve the forward problem
        X = self.forwardSolve()

        self.S_hat, self.I_hat, self.R_hat =  X[:, 0], X[:, 1], X[:, 2]
        print(np.array([ self.S_hat, self.I_hat, self.R_hat]).T)


        λ = self.backwardsSolve()

        β, γ  = self.p[0:2]
        N = self.N

        dλ = []

        for t in self.time_steps:

            I  = self.I[t]
            I_hat = self.I_hat[t]
            S_hat = self.S_hat[t]
            N = self.N

            λ_t = np.array(λ[t])
    
            df_dx = -2*np.array([0, I - I_hat, 0])
            dh_dx = np.array([[β*I_hat/N  ,      β*S_hat/N, 0],
                            [-β*I_hat/N  , -β*S_hat/N + γ, 0],
                            [0          ,             -γ, 0]])
            
            dλ.append((df_dx + λ_t.T@dh_dx)) 

        dλ = np.array(dλ)
        dt = self.time_steps[1]-self.time_steps[0]

   
        #the solution to the last coordinate is 0, so we add some noise to avoid
        #division errors. 
        df = (λ[1:,] - λ[0:-1,])/dt + .0001
        print(λ)
        print(df)
        print(dλ)
        #print(dX[0])

        diffs = np.abs((dλ[1:,] - df)/df)

        return np.mean(diffs), df,  λ, dλ

