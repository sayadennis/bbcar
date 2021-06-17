import os
import sys
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from torch import nn

class BBCarMLP(torch.nn.Module):
        def __init__(self, input_size, hidden_size):
            super(BBCarMLP, self).__init__()
            self.input_size = input_size
            self.hidden_size  = hidden_size
            self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
            self.relu = torch.nn.ReLU()
            self.fc2 = torch.nn.Linear(self.hidden_size, 1)
            self.sigmoid = torch.nn.Sigmoid()
            ## alternative A 
            # self.layers = nn.Sequential(
            #     nn.Flatten(),
            #     nn.Linear(32 * 32 * 3, 64),
            #     nn.ReLU(),
            #     nn.Linear(64, 32),
            #     nn.ReLU(),
            #     nn.Linear(32, 10)
            #     )
            ## alternative B
            # self.lin1 = nn.Linear(784, 512, bias=True) 
            # self.lin2 = nn.Linear(512, 256, bias=True)
            # self.lin3 = nn.Linear(256, 10, bias=True)
            ## alternative C
            # self.dropout = nn.Dropout(0.2)


        
        def forward(self, x):
            hidden = self.fc1(x)
            relu = self.relu(hidden)
            output = self.fc2(relu)
            output = self.sigmoid(output)
            ## alternative A 
            # return self.layers(x)
            ## alternative B
            # x = xb.view(-1,784) 
            # x = F.relu(self.lin1(x))
            # x = F.relu(self.lin2(x))
            # return self.lin3(x)
            return output


from sklearn.datasets import make_blobs

def blob_label(y, label, loc): # assign labels
    target = np.copy(y)
    for l in loc:
        target[y == l] = label
    return target

x_train, y_train = make_blobs(n_samples=40, n_features=2, cluster_std=1.5, shuffle=True)
x_train = torch.FloatTensor(x_train)
y_train = torch.FloatTensor(blob_label(y_train, 0, [0]))
y_train = torch.FloatTensor(blob_label(y_train, 1, [1,2,3]))
x_test, y_test = make_blobs(n_samples=10, n_features=2, cluster_std=1.5, shuffle=True)
x_test = torch.FloatTensor(x_test)
y_test = torch.FloatTensor(blob_label(y_test, 0, [0]))
y_test = torch.FloatTensor(blob_label(y_test, 1, [1,2,3]))

model = BBCarMLP(2, 10)
criterion = torch.nn.BCELoss()
optimizer = torch.optim.SGD(model.parameters(), lr = 0.01)

model.eval()
y_pred = model(x_test)
before_train = criterion(y_pred.squeeze(), y_test)
print('Test loss before training' , before_train.item())

model.train()
epoch = 20
for epoch in range(epoch):
    optimizer.zero_grad()
    # Forward pass
    y_pred = model(x_train)
    # Compute Loss
    loss = criterion(y_pred.squeeze(), y_train)
   
    print('Epoch {}: train loss: {}'.format(epoch, loss.item()))
    # Backward pass
    loss.backward()
    optimizer.step()

model.eval()
y_pred = model(x_test)
after_train = criterion(y_pred.squeeze(), y_test) 
print('Test loss after Training' , after_train.item())

####

from torchvision.datasets import CIFAR10
from torchvision import transforms

if __name__ == '__main__':
    # Set fixed random number seed
    torch.manual_seed(42)
    # Prepare CIFAR-10 dataset
    dataset = CIFAR10(os.getcwd(), download=True, transform=transforms.ToTensor())
    trainloader = torch.utils.data.DataLoader(dataset, batch_size=10, shuffle=True, num_workers=1)
    # Initialize the MLP
    mlp = BBCarMLP()
    # Define the loss function and optimizer
    loss_function = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(mlp.parameters(), lr=1e-4)
    # Run the training loop
    for epoch in range(0, 5): # 5 epochs at maximum
        # Print epoch
        print(f'Starting epoch {epoch+1}')
        # Set current loss value
        current_loss = 0.0
        # Iterate over the DataLoader for training data
        for i, data in enumerate(trainloader, 0):        
            # Get inputs
            inputs, targets = data
            # Zero the gradients
            optimizer.zero_grad()
            # Perform forward pass
            outputs = mlp(inputs)
            # Compute loss
            loss = loss_function(outputs, targets)
            # Perform backward pass
            loss.backward()
            # Perform optimization
            optimizer.step()
            # Print statistics
            current_loss += loss.item()
            if i % 500 == 499:
                print('Loss after mini-batch %5d: %.3f' %
                        (i + 1, current_loss / 500))
                current_loss = 0.0
    # Process is complete.
    print('Training process has finished.')