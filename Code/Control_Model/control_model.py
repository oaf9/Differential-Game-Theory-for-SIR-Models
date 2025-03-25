import torch
import numpy as np

 #ode_adjoint expects a forward function that inherits from nn.module

def solveQuadratic(a,b,c):
    """ Solves ax^2 + bx +c = 0
        Note: this solver assumes the existence of two real roots
    """
    #a,b,c = a.detach(), b.detach(), c.detach()
    x1 = -b + torch.sqrt(b*b - 4*a*c)
    x2 = -b - torch.sqrt(b*b - 4*a*c)

    #print(torch.sqrt(b*b - 4*a*c))

    return torch.stack([x1,x2])/(2*a)

class control_model(torch.nn.Module):
    def __init__(self, θ, m):
        super(control_model, self).__init__()

        """Constructor: θ = (Β_0, γ_0,  m,  S_0, I_0, V_0)
        """
        self.θ = torch.nn.Parameter(θ)
        self.m = torch.nn.Parameter(m)

        self.cost = []

    def σ(self, c):
        """ models the returns the risk of infection in terms of control"""
        return 1/(self.m*c + 1)

    def forward(self, t, X):

        # defines the SIR equations
        Β, γ  = self.θ[0:2]

        S, I, R, Vs = X
        N = S + I + R

        #update the control values
        a,b,c = (self.m*self.m), 2*self.m, (1-self.m*I*(1+Vs))

        c_optimal = 0

        if c >= 0: 
            c_optimal = solveQuadratic(a,b,c)[0]

        σ_c = self.σ(c_optimal)

        
        dS = -σ_c*(Β/N)*S*I
        dI = σ_c*(Β/N)*S*I - γ*I
        dR = γ*I
        dVs = (1+Vs)* σ_c * I + c_optimal

        return torch.stack([dS, dI, dR, -dVs])
    

    
    def get_initial_conditions(self):

        #this function is just a helpful getter ... could be removed without 
        #affecting anything except some syntactical stuffs.
        return self.θ[2:]
    

    

