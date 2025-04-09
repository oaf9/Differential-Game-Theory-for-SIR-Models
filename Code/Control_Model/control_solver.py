from control_model import control_model
from torchdiffeq import odeint_adjoint as odeint
import torch

#this just returns the MSE
def lossF(pred, true): 
        return torch.mean((pred - torch.tensor(true, dtype=torch.float32))**2)

def fit(p_0, m, X, epochs = 100, lr = .0001):
    """
    X is the data = [S,I,R]
    """

    #first we ensure that the initial paramaters are compatible with torch
    # p_0 = (Β_0, γ_0, S_0, I_0, R_0, V_0 ) -- matches params for SIR()
    p_0 = torch.tensor(p_0, dtype = torch.float32)
    m = torch.tensor(m, dtype = torch.float32)
    
    N = p_0[2:].sum()

    #setting the timesteps to fit at ... this will depend on your data
    t = torch.linspace(0, len(X), 152)

    #initilize the model
    model = control_model(p_0, m, len(X))

    #using stochastic gradient descent 
    optimizer = torch.optim.Adam(model.parameters(), lr = lr)


    for i in range(epochs): 

        print(model.m)

        #make the forward pass: 

        # retrieve the initil conditions from the model
        X_0 = model.get_initial_conditions()

        #we are optimizing to minimize error in I and S 
        X_hat = odeint(model, X_0, t, method = 'dopri5', 
                       options={'step_size': 1e-3, 'min_step': 1e-6})[:,0:2]
        X_true = X[:,0:2]
        
        #calculate loss
        loss = lossF(X_hat[:, 0:2], X_true[:, 0:2])
        if(loss < .001): #if we get RSE low enough, we break the loop
             break

        #compute and update the gradients
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        #project the values back into valid numerical ranges
        with torch.no_grad():

            model.θ.data[0] = model.θ.data[0].clamp(0.0001,5)
            model.θ.data[1] = model.θ.data[1].clamp(0.0001,1)
            model.θ.data[2:] = model.θ.data[2:].clamp(0.0001, N)
            model.θ.data[2:] = model.θ.data[2:] / model.θ.data[2:].sum() * N

        #print values
        print(f"p = {model.θ}, {model.m}")
        print(f"Loss {i} = {loss}")


    # make a forward pass and return those values
    return odeint(model, model.θ[2:], t, method='dopri5'), model
