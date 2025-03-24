import torch

 #ode_adjoint expects a forward function that inherits from nn.module
class SIR(torch.nn.Module):
    def __init__(self, θ):
        super(SIR, self).__init__()

        """Constructor: θ = (Β_0, γ_0, S_0, I_0, R_0)
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
    

    

