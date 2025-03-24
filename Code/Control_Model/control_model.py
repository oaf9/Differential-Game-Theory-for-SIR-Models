import torch
import numpy as np

 #ode_adjoint expects a forward function that inherits from nn.module

def solveQuadratic(a,b,c):
    """ Solves ax^2 + bx +c = 0
        Note: this solver assumes the existence of two real roots
    """
    x1 = -b + np.sqrt(b*b - 4*a*c)
    x2 = -b - np.sqrt(b*b - 4*a*c)

    return 2*a*np.array(x1,x2)


class control_model(torch.nn.Module):
    def __init__(self, θ):
        super(control_model, self).__init__()

        """Constructor: θ = (Β_0, γ_0, m,  S_0, I_0, V_0)
        """
        self.θ = torch.nn.Parameter(θ) 

    def forward(self, t, X):
        # defines the SIR equations
        Β, γ  = self.θ[0:2]
        S, I, R = X

        N = S + I + R

        dS = -(Β/N)*S*I
        dI = (Β/N)*S*I - γ*I
        dR = γ*I

        return torch.stack([dS, dI, dR])
    
    def get_initial_conditions(self):

        #this function is just a helpful getter ... could be removed without 
        #affecting anything except some syntactical stuffs.

        return self.θ[2:]
    

    

