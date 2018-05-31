import numpy as np


gamma = 2
a = 1


# @profile
def dehnen(n):
    out = []
    fac1 = (3 - gamma) / (4 * np.pi)/a**3/7705
    rd = np.random.rand(n)
    while len(out) < n:
        r = 10*np.random.rand(n)
        r_1 = r+1
        p = fac1 / (np.power(r/r_1, gamma) * np.power(r_1, 4))
        k = rd-p
        p = np.where(k < 0)[0]
        out.extend(r[p])
    out = np.array(out)
    out = out[:n]
    return out


def center_distance(n):
    # r = np.random.randn(n)
    r = dehnen(n)
    return np.abs(r)


# @profile
def positions(n):
    pos = np.zeros((n, 3))
    pos[:, 0] = center_distance(n)
    pos[:, 1] = 2*np.pi*np.random.rand(n)
    pos[:, 2] = np.pi*np.random.rand(n)
    return pos


# @profile
def positions_2d(n, scale):
    pos3d = positions(n)
    pos2d = np.zeros((pos3d.shape[0], 2))
    print(np.max(pos3d[:, 0]))
    x = scale*pos3d[:, 0]*np.cos(pos3d[:, 1])
    y = scale*pos3d[:, 0]*np.sin(pos3d[:, 1])
    print(np.min(x), np.min(y))
    print(np.max(x), np.max(y))
    pos2d[:, 0] = x
    pos2d[:, 1] = y
    return pos2d


if __name__ == '__main__':
    positions_2d(1000, 5)
