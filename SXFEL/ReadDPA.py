import numpy as np
import matplotlib.pylab as plt

# Change parameter according to your .in file
nslice = 1100  # set to 1 if run in steady-state mode
npart = 8192  # number fo particles in each slice
harmonic = 1
# lammbda = 267e-9
# R56 = 115e-6

# Creat empty matrix
tau = np.empty(npart*nslice)  # "tau" = c*t
p = np.empty(npart*nslice)  # 'p' is canonical momentum E/(c*p0)
x = np.empty(npart*nslice)  # x is horizontal cartesian coordinate
xp = np.empty(npart*nslice)  # 'xp' is conjugate momentum canonical momentum px/p0
y = np.empty(npart*nslice)
yp = np.empty(npart*nslice)

with open('real_current.out.dpa', 'rb') as f:
    for i in range(nslice):
        # gamma = np.empty(nslice)
        # dgamma = np.empty(nslice)
        # sigmax = np.empty(nslice)
        # sigmay = np.empty(nslice)

        single_slice = np.fromfile(f, dtype=np.double, count=npart*6)
        single_slice = single_slice.reshape(6,-1)

        # gamma[i] = np.mean(single_slice[0,:])
        # dgamma[i] = np.sqrt(np.sum((single_slice[0,:]-gamma[i])**2)/npart)
        # sigmax[i] = np.sqrt(np.mean((single_slice[2,:] - np.mean(single_slice[4,:]))**2))
        # sigmay[i] = np.sqrt(np.mean((single_slice[3,:] - np.mean(single_slice[5,:]))**2))

        tau[i*npart: (i+1)*npart] = single_slice[1,:] / harmonic + 2*np.pi*i
        p[i*npart: (i+1)*npart] = single_slice[0,:]
        x[i*npart: (i+1)*npart] = single_slice[2,:]
        xp[i*npart: (i+1)*npart] = single_slice[4,:] / single_slice[0,:]
        y[i*npart: (i+1)*npart] = single_slice[3,:]
        yp[i*npart: (i+1)*npart] = single_slice[5,:] / single_slice[0,:]

# Visualization
plt.figure()
plt.plot(tau/2/np.pi, p, 'r.', markersize=0.01)
plt.xlabel(r'$z/\lambda$')
plt.ylabel('p')
plt.title('phase space')
plt.show()

# tau = tau / 2 / np.pi * lammbda
# tau = R56 * (p - np.nanmean(p)) / np.nanmean(p)

# plt.figure()
# plt.plot(tau, p, 'r.', markersize=0.1)
# plt.xlabel('tau/m')
# plt.ylabel('p')
# plt.show()

# plt.figure()
# plt.hist(tau)
# plt.show()

# plt.figure()
# plt.plot(gamma, 'r-', markersize=0.1)
# plt.show()

# plt.figure()
# plt.plot(dgamma, 'r-', markersize=0.1)
# plt.show()

# plt.figure()
# plt.plot(sigmax, 'r-', markersize=0.1)
# plt.show()

# plt.figure()
# plt.plot(sigmay, 'r-', markersize=0.1)
# plt.show()

# Calcualte Twiss
gamma = np.mean(p)
E = gamma*0.511
xp1 = xp - np.mean(xp)
yp1 = yp - np.mean(yp)
emit_x = (np.mean(x**2) * np.mean(xp1**2) - np.mean(x * xp1)**2)**0.5 * gamma
emit_y = (np.mean(y**2) * np.mean(yp1**2) - np.mean(y * yp1)**2)**0.5 * gamma
beta_x = np.mean(x**2) * gamma / emit_x
beta_y = np.mean(y**2) * gamma / emit_y
gamma_x = np.mean(xp**2) * gamma / emit_x
gamma_y = np.mean(yp**2) * gamma / emit_y
alpha_x = (beta_x * gamma_x - 1)**0.5
alpha_y = -(beta_y * gamma_y - 1)**0.5

print('gamma: %f' %gamma)
print('E: %f' %E)
print('emit_x: %e' %emit_x, 'emit_y: %e' %emit_y)
print('beta_x: %f' %beta_x, 'beta_y: %f' %beta_y)
print('gamma_x: %f' %gamma_x, 'gamma_x: %f' %gamma_y)
print('alpha_x: %f' %alpha_x, 'alpha_y: %f' %alpha_y)



'''
yp1=yp-mean(yp);
emy0=(mean(y.^2)*mean(yp1.^2)-mean(y.*yp1)^2)^0.5*gamma;
xp1=xp-mean(xp);
emx0=(mean(x.^2)*mean(xp1.^2)-mean(x.*xp1)^2)^0.5*gamma;
gamma=mean(g); 
beta0x=mean(x'.^2).*gamma./emx0
beta0y=mean(y'.^2).*gamma./emy0
% 
gamma0x=mean(xp.^2).*gamma./emx0;
alpha0x=(beta0x*gamma0x-1).^0.5
gamma0y=mean(yp.^2).*gamma./emy0; 
alpha0y=-(beta0y*gamma0y-1).^0.5
'''
