from SIR_Adjoint_model import SIR
from torchdiffeq import odeint_adjoint as odeint
import torch

#this just returns the MSE
def lossF(pred, true): 
        return torch.mean((pred - torch.tensor(true, dtype=torch.float32))**2)

def fit(p_0, X, epochs = 100, lr = .001):

    #first we ensure that the initial paramaters are compatible with torch
    # p_0 = (Β_0, γ_0, S_0, I_0, R_0 ) -- matches params for SIR()
    p_0 = torch.tensor(p_0, dtype = torch.float32)
    N = p_0[2:].sum()

    #setting the timesteps to fit at ... this will depend on your data
    t = torch.linspace(0, len(X), 152)

    #initilize the model
    model = SIR(p_0)

    #using stochastic gradient descent 
    optimizer = torch.optim.SGD(model.parameters(), lr = lr)

    for i in range(epochs): 

        #make the forward pass: 
        # retrieve the initil conditions from the model
        
        X_0 = model.get_initial_conditions()
        X_hat = odeint(model, X_0, t, method = 'dopri5')[:,0:2]
        
        X_true = X[:,0:2]
        
        #calculate loss
        loss = 0.1 * lossF(X_hat[:, 0], X_true[:, 0]) + 0.9 * lossF(X_hat[:, 1], X_true[:, 1])



        #compute and update the gradients
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        #project the values back into valid numerical ranges
        with torch.no_grad():
            model.θ.data[0] = model.θ.data[0].clamp(0.001,5)
            model.θ.data[1] = model.θ.data[1].clamp(0.001,1)
            model.θ.data[2:] = model.θ.data[2:].clamp(0.001, N)
            model.θ.data[2:] = model.θ.data[2:] / model.θ.data[2:].sum() * N

        #print values
        print(f"p = {model.θ}")
        print(f"Loss {i} = {loss}")


    # make a forward pass and return those values
    return odeint(model, model.θ[2:], t, method='dopri5'), model
