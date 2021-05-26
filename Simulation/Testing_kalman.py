import numpy as np

At = np.array([[0, 1], [0, 0]])
Bt = np.array([0], [-1])
Ts = 0.01
Ak = np.array([[1, Ts], [0, 1]])
Bk = np.array([-(Ts**2)/2],[-Ts])
Ck = np.array([1, 0])

altitude = 2000.0

TN = 20

def change_x(x):
    return x + np.array([dx1, dx2])

def GPS(x):
    x += wind
    x += np.random.normal(0, sigma_x, x.shape)
    return x

x1_0 = 0.0
x2_0 = 0.0

if __name__ == '__main__':
    x = np.array([x1_0, x2_0])
    x_meas = x
    for i in range(100):
        x = change_x(x)
        x_meas = GPS(x)
