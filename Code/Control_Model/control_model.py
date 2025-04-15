import torch
import numpy as np


def solveQuadratic(a,b,c):
    """ Solves ax^2 + bx +c = 0
        Note: this solver assumes the existence of two real roots
    """
    #a,b,c = a.detach(), b.detach(), c.detach()
    x1 = -b + torch.sqrt(b*b - 4*a*c)
    x2 = -b - torch.sqrt(b*b - 4*a*c)

    #print(torch.sqrt(b*b - 4*a*c))
    return torch.stack([x1,x2])/(2*a)


def Linear_interpolation(t, x_prev, x_next, y_prev, y_next):
    """Simple Linear Interpolater
    """
    a = y_prev
    b = (y_next - y_prev)/(x_next - x_prev)

    return a + b*(t - x_prev)


#ode_adjoint expects a forward function that inherits from nn.module
class control_model(torch.nn.Module):
    def __init__(self, θ, m, time_steps):
        super(control_model, self).__init__()

        """ Constructor: θ = (Β_0, γ_0, m, S_0, I_0, V_0) """

        self.θ = torch.nn.Parameter(θ)
        self.m = torch.nn.Parameter(m)

        #these are the control paramaters
        self.c = torch.nn.Parameter(torch.zeros(time_steps)+.1)

    def σ(self, c):
        """ models the returns the risk of infection in terms of control"""
        return 1/(self.m*c + 1)
    

    def interpolate(self, t):

        t_max = self.c.shape[0] - 1
        if t <= 0:
            return self.c[0]
        elif t >= t_max:
            return self.c[-1]
        else: 
            t_prev = torch.floor(t).long()
            t_next = torch.ceil(t).long()
            y_prev = self.c[t_prev]
            y_next = self.c[t_next]

            return Linear_interpolation(t = t, 
                                        x_prev = t_prev,
                                        y_prev = y_prev,
                                        x_next = t_next, 
                                        y_next = y_next)



    def forward(self, t, X):
        """defines the SIR equations"""

        #unpack model params
        Β, γ  = self.θ[0:2]

        #unpack the system state
        S, I, R, Vs = X
        N = S + I + R

        #update the control values
        c_optimal = torch.relu(self.interpolate(t))

        #pass c(t) through σ()
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
    

    

