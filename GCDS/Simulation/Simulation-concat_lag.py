#!/usr/bin/env python
# coding: utf-8

# In[3]:


import torch
import torch.nn as nn
import torch.optim as optim
import torchvision.utils as utils
import torchvision.datasets as dsets
from torchvision import datasets
from torch.utils.data import Dataset
import torchvision.transforms as transforms
from torchvision.utils import save_image

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ## load the dataset 

# In[4]:


data = pd.read_csv("Simulation_data.csv")


# In[5]:


data.head()


# In[6]:


data.shape


# In[7]:


def imshow(img):
    img = (img+1) / 2
    img = img.squeeze() # 차원 중 사이즈 1 을 제거
    np_img = img.numpy() # 이미지 픽셀을 넘파이 배열로 변환
    plt.imshow(np_img)
    plt.show()


# In[8]:


# for i in range(30):
#     plt.imshow(data.iloc[i,].to_numpy().reshape((30,30)))
#     plt.show()


# ## data preprocessing

# In[9]:


# place holder
dataset = torch.zeros([435,2,900])
lag_lst = torch.zeros([435])


# In[10]:


# all possible pairs of idices
idx_list = range(30)
pair = [(a, b) for idx, a in enumerate(idx_list) for b in idx_list[idx + 1:]]


# In[11]:


for i in range(len(pair)): 
    idx1 = pair[i][0]
    idx2 = pair[i][1]
    lag = idx2-idx1
    
    dataset[i,0,:] = torch.tensor(data.iloc[idx1,].values)
    dataset[i,1,:] = torch.tensor(data.iloc[idx2,].values)
    lag_lst[i] = (lag)
    


# done; dataset and lag_lst

# In[12]:


print(dataset[0,0,:].min())
print(dataset[0,0,:].max())


# In[ ]:





# In[13]:


dataset = (dataset - dataset.min())/(dataset.max()- dataset.min())

print(dataset.min())
print(dataset.max())


# In[14]:


dataset = ((dataset - 0.5)/0.5)

print(dataset.min())
print(dataset.max())


# scaled dataset in [-1,1]

# ## Data loader

# 나중에 모델만 짜는 파일을 만들어서 import 하는 형태로 구성하자 -> training 하는 과정을 함수로 만들어야 의미가 있을 듯 하다. 
# 

# In[15]:


device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"device = {device}")
sample_dir = 'samples'
if not os.path.exists(sample_dir):
    os.makedirs(sample_dir)


# X를 이미지 2개를 받아서 만드는 것으로, y는 label을 받는 대신 lag를 받는 것으로 변경

# In[16]:


class CustomImageDataset(Dataset):
    def __init__(self,X,Y):
        self.X = X
        self.Y = Y
        
    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        X = self.X
        Y = self.Y
        return X[idx], Y[idx]


# In[17]:


dataset = dataset.resize(435,2,30,30)


# In[18]:


trainset = CustomImageDataset(dataset,lag_lst)


# In[19]:


print(trainset.__getitem__(1)[0][0].size())
print(trainset.__getitem__(1)[1])


# In[20]:


plt.imshow(trainset.__getitem__(1)[0][0].numpy())
plt.show()


# In[21]:


plt.imshow(trainset.__getitem__(1)[0][1].numpy())
plt.show()


# In[22]:


batch_size= 50


# In[23]:


# data loader
data_loader = torch.utils.data.DataLoader(dataset=trainset, 
                                          batch_size=batch_size,
                                          shuffle=True) 


# In[24]:


for X,Y in data_loader:
    print(X.size())
    print(Y.size())
    


# ### Model construction

# In[25]:


# Channel changed
class Discriminator(nn.Module):
    def __init__(self):
        super(Discriminator, self).__init__()
        self.D1 = nn.Sequential(
            nn.Conv2d(1,64,5,2,padding=1), 
            nn.BatchNorm2d(64),
            nn.LeakyReLU(0.2),
            nn.Conv2d(64,128,4,1,padding=1), 
            nn.BatchNorm2d(128),
            nn.LeakyReLU(0.2),
            nn.Flatten(1)
            )
        self.D2 = nn.Sequential(
            nn.Conv2d(1,32,5,2,padding=1), 
            nn.BatchNorm2d(32),
            nn.LeakyReLU(0.2),
            nn.Conv2d(32,64,5,2,padding=1),
            nn.BatchNorm2d(64),
            nn.LeakyReLU(0.2),
            nn.Flatten(1) 
            )
        self.D3 = nn.Sequential(
            nn.Linear(23965,128), 
            nn.BatchNorm1d(128),
            nn.LeakyReLU(0.2),
            nn.Linear(128,64),
            nn.BatchNorm1d(64),
            nn.LeakyReLU(0.2),
            nn.Linear(64,1) # [100, 1]
            )
    def forward(self,X,Y,conditional,device=device):
        #one hot vector
        lag_lst = conditional.int().tolist()
        max_lag = 29 # int(lag.max().item())
        lag_idx = [lag_lst[i]-1 for i in range(len(lag_lst))]
        one_hot = nn.functional.one_hot(torch.arange(0, max_lag) % max_lag)[lag_idx]
        
        disc_X = self.D1(X)
        disc_Y = self.D2(Y)
        Dout1 = torch.concat((disc_X,disc_Y),1)
        Dout2 = torch.concat((Dout1,one_hot.to(device)),1)
        Dout3 = self.D3(Dout2)
        
        return Dout3


# In[26]:


for X,Y in data_loader:
    disc_X = X[:,0,:].reshape(batch_size,1,30,30)
    disc_Y = X[:,1,:].reshape(batch_size,1,30,30)
    print(Y.size())
    lag = Y
    break
    


# In[27]:


disc_X.size()


# In[28]:


disc_Y.size()


# In[29]:


# Discriminator check
D = Discriminator().to(device)
out = D(disc_X.to(device),disc_Y.to(device),conditional=lag)

print(out.size())
print(out.device)


# In[30]:


plt.imshow(disc_X[1].squeeze().numpy())
plt.show()


# In[31]:


# lag_lst : Y in data_loader 


# In[32]:


class Generator(nn.Module):
    def __init__(self):
        super(Generator, self).__init__()
        self.G1 = nn.Sequential(
            nn.Flatten(1), # 
            nn.Linear(900, 128), # 100*128
            nn.BatchNorm1d(128), # 100*128
            nn.LeakyReLU(0.2)) # 100*128
        self.G2 = nn.Sequential(
            nn.ConvTranspose2d(257,256,kernel_size=(7,7),stride=(1,1),padding=0), # [100, 256, 7, 7]
            nn.BatchNorm2d(256), # 
            nn.LeakyReLU(0.2), # [100, 256, 7, 7]
            nn.ConvTranspose2d(256,128,kernel_size=(7,7),stride=(1,1)), # [100, 128, 14, 14]
            nn.BatchNorm2d(128), # 
            nn.LeakyReLU(0.2),#
            nn.ConvTranspose2d(128,1,kernel_size=(6,6),stride=(2,2)), # [100, 1, 14, 14] -> should be matched with size of y.
            nn.Tanh()# [100, 1, 14, 14]
            )

    def forward(self,X,eta_dist="Normal",eta_dim=100,conditional = lag,device=device):
        # random error part
        if eta_dist == "Normal":
            eta = torch.randn(X.size()[0], eta_dim, device=device) # X.size()[0] :batch_size

        elif eta_dist == "Unif":
            eta = torch.rand(X.size()[0], eta_dim, device=device)

        else :
            raise Exception("Please specify the correct distribution. Normal or Unif")
        
        #one hot vector
        lag_lst = conditional.int().tolist()
        max_lag = 29 # int(lag.max().item())
        lag_idx = [lag_lst[i]-1 for i in range(len(lag_lst))]
        one_hot = nn.functional.one_hot(torch.arange(0, max_lag) % max_lag)[lag_idx]

        # construcnt the full model
        Gout1 = torch.concat((self.G1(X),one_hot.to(device)),1)
        Gout2 = torch.concat((Gout1,eta),1).resize(X.size()[0],(128 + max_lag + eta_dim),1,1)
        Gout3 = self.G2(Gout2)
    
        return Gout3


# In[33]:


int(lag.max().item())


# In[34]:


# Generator check

print(disc_X.size())
G = Generator().to(device)
out = G(X=disc_X.to(device),eta_dist="Normal",eta_dim=100,conditional=lag)
print(out.size())
print(out.device)


# In[35]:


# visualizing the generator image

rand = torch.randn(2,1, 30,30, device=device)
G = Generator().to(device)
img_1 = G(rand,conditional=torch.tensor([1,2]))[0] # 범위가 [-1,1]임을 감안
imshow(img_1.squeeze().cpu().detach())


# ### Loss function

# In[36]:


# Defince the fenchel conjugacy loss 

# def disc_loss(output,label): # label : 1 for real, -1 for fake
#     loss = - torch.mean( output * (label - 1.) * (-1.0/2.0) - torch.exp(output) * (label + 1.) * (1.0/2.0) )
#     return loss

def disc_loss(output,label): # label : 1 for real, -1 for fake
    loss = - torch.mean( output * (label + 1.) * (1.0/2.0) - torch.exp(output) * (label - 1.) * (-1.0/2.0) )
    return loss

def gen_loss(output,label): # all labels are -1.(all fake data)
    loss = torch.mean(output)
    return loss


# ### Visualizing tool for the results

# In[37]:


train_X = dataset[:,0,:].unsqueeze(1)
train_Y = dataset[:,1,:].unsqueeze(1)
print(train_X.size())


# In[40]:


# making plots 

    
def imshow_sub(img,ax,i,j):
    img = (img+1) / 2
    img = img.squeeze() # 차원 중 사이즈 1 을 제거
    np_img = img.numpy() # 이미지 픽셀을 넘파이 배열로 변환
    ax[i,j].imshow(np_img)
    ax[i,j].axis('off')

def making_plots(X=train_X,Y=train_Y,lag=lag_lst,device=device,G=G,sample_idx = [1,10,8,30,20,100],rep=5):
    sample_batch = X[sample_idx].clone()
    sample_batch_Y = Y[sample_idx]
    sample_batch = (sample_batch*0.5 + 0.5)
    sample_lag = lag[sample_idx]
    # X,Y를 같이 나열 
    fig, ax = plt.subplots(len(sample_idx),rep + 2, figsize=(15,15))

    for i in range(len(sample_idx)):
        for j in range(rep + 2):
            if j == 0:
                imshow_sub(sample_batch[i].detach(),ax,i,j)
                #ax[i,j].set_title(str(i)+str(" ")+str(j))
                ax[i,j].set_title("lag= " + str(int(lag[sample_idx[i]].item())))
                #ax[i,j].set_title(sample_idx[i])
            elif j==1:
                imshow_sub(sample_batch_Y[i].detach(),ax,i,j)
                ax[i,j].axis('off')
            else:
                res_y_tilde = G(sample_batch.to(device),conditional=sample_lag)
                imshow_sub(res_y_tilde[i].detach().cpu(),ax,i,j)
                ax[i,j].axis('off')
    plt.show()


# In[41]:


making_plots(X=train_X,Y=train_Y,device=device,G=G,sample_idx = [1,390,28,250,30,56,65,90,434])


# ### Training

# In[45]:


num_epochs = 3500


# In[ ]:


# training

real_scores = []
fake_scores = []
gloss_lst = []
dloss_lst = [] 

total_step = len(data_loader)
D = Discriminator().to(device)
G = Generator().to(device)

d_optimizer = torch.optim.Adam(D.parameters(), lr=0.00001)
g_optimizer = torch.optim.Adam(G.parameters(), lr=0.0002)

for epoch in range(num_epochs):
    for i, (XY,lag_lst) in enumerate(data_loader):
        X = XY[:,0,:].to(device).unsqueeze(1)
        Y = XY[:,1,:].to(device).unsqueeze(1)
        
        # Create the labels which are later used as input for loss function
        real_labels = -torch.ones(len(XY), 1).to(device)
        fake_labels = torch.ones(len(XY), 1).to(device)

        # ================================================================== #
        #                      Train the discriminator                       #
        # ================================================================== #

        # Compute Loss using real images where custom_loss(x, y)
        
        Dout = D(X,Y,conditional=lag_lst)
        d_loss_real = disc_loss(Dout,real_labels)
        real_score = Dout
        
        # Compute Loss using fake images
        
        Gout = G(X,eta_dist="Normal",eta_dim=100,conditional=lag_lst)
        DGout = D(X,Gout,conditional=lag_lst)
        d_loss_fake = disc_loss(DGout, fake_labels)
        
        # Backprop and optimize
        d_loss = d_loss_real + d_loss_fake
        d_optimizer.zero_grad()
        d_loss.backward(retain_graph=True)
        d_optimizer.step()
        
        # ================================================================== #
        #                        Train the generator                         #
        # ================================================================== #
    
        # Compute loss with fake images
        D.requires_grad = False
        DGout2 = D(X,Gout,conditional=lag_lst)  # G(X,eta_dist="Normal",eta_dim=200)
        fake_score = DGout2

        # Backprop and optimize
        g_loss = gen_loss(DGout2, -1)
        g_optimizer.zero_grad()
        g_loss.backward(retain_graph=True)
        g_optimizer.step()
        
        
    if (epoch) % 50 == 0:
        print('Epoch [{}/{}], Step [{}/{}], d_loss: {:.4f}, g_loss: {:.4f}, D(x): {:.2f}, D(G(z)): {:.2f}' 
              .format(epoch, num_epochs, i+1, total_step, d_loss.item(), g_loss.item(), 
                      real_score.mean().item(), fake_score.mean().item()))
        making_plots(X=train_X,Y=train_Y,device=device,G=G,sample_idx = [1,27,28,29,30,56,65,90,103])
            
    real_scores.append(real_score.mean().item())            
    fake_scores.append(fake_score.mean().item())
    gloss_lst.append(g_loss.item())
    dloss_lst.append(d_loss.item())
    
    
#     # real image 저장
#     if (epoch+1) == 1:
#         images = images.reshape(images.size(0), 1, 28, 28)
#         save_image(denorm(images), os.path.join(sample_dir, 'real_images.png'))
    
#     # 생성된 이미지 저장
#     fake_images = fake_images.reshape(fake_images.size(0), 1, 28, 28)
#     save_image(denorm(fake_images), os.path.join(sample_dir, 'fake_images-{}.png'.format(epoch+1)))

# 생성자, 판별자 각각 모델 저장
torch.save(G.state_dict(), 'G.ckpt')
torch.save(D.state_dict(), 'D.ckpt')


# In[ ]:


mdloss = [dloss_lst[i]*(-1) for i in range(len(dloss_lst))]


# In[ ]:


# plot    
plt.figure(figsize = (12, 8))
plt.xlabel('epoch')
plt.ylabel('loss')
x = np.arange(num_epochs)  # num_epochs
plt.plot(x, gloss_lst, 'g', label='generator loss')
plt.plot(x, mdloss, 'b', label='discriminator loss')
plt.title("loss plot")
plt.legend()
plt.show()


# In[ ]:


# plot    
plt.figure(figsize = (12, 8))
plt.xlabel('epoch')
plt.ylabel('score')
x = np.arange(num_epochs)
plt.plot(x, real_scores, 'g', label='D(x,y)')
plt.plot(x, fake_scores, 'b', label='D(x,G(x,eta))')
plt.legend()
plt.show()


# In[ ]:




