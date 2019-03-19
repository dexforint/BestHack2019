import pandas as pd
import numpy as np
from math import pi, sin, cos, sqrt
from scipy.interpolate import make_interp_spline
import argparse

F = pd.read_csv('F.csv')
Wind = pd.read_csv('Wind.csv')

parser = argparse.ArgumentParser()
parser.add_argument(
    '-V0',
    type=float,
    required=True,
    help='Начальная скорость'
)
parser.add_argument(
    '-H0',
    type=float,
    required=True,
    help='Начальная высота'
)
parser.add_argument(
    '-m',
    type=float,
    required=True,
    help='Масса'
)
parser.add_argument(
    '-R',
    type=float,
    default=1,
    help='Радиус'
)
parser.add_argument(
    '-dt',
    type=float,
    default=0.01,
    help='Временной шаг'
)
parser.add_argument(
    '-a',
    type=str,
    default='0',
    help='Углы'
)
args = parser.parse_args()

V0 = args.V0
H0 = args.H0
m = args.m
R = args.R
angles = [float(i) for i in args.a.split(',')]
dt = args.dt


class simulation():
    def __init__(self, H0, V0, m, F, Wind, R=1):
        self.h0 = H0
        self.v0 = V0
        self.m = m
        self.R = R
        S = pi * R**2
        C = 0.47
        p = 1.225
        self.k = C * S * p / 2
        
        V = F['V']
        Fa = F['Fa']
        Y = Wind['Y']
        Wx = Wind['Wx']
        Wz = Wind['Wz']
        self.f_interp = make_interp_spline(V, Fa)
        self.Wx_interp = make_interp_spline(Y, Wx)
        self.Wz_interp = make_interp_spline(Y, Wz)
    
    def emulate(self, angles=[0], dt=0.01):
        g = 9.81
        history = []
        for i, angle in enumerate(angles):
            Coords = np.array([0, 0, self.h0], dtype=np.float32)
            V = np.array([self.v0 * cos(angle), self.v0 * sin(angle), 0])
            a0 = np.array([0, 0, -g])
            F = self.aerodynamic_power(V) + self.wind_power(V, self.h0)
            a = a0 + F / self.m
            t = 0.
            dots = [[t, Coords.copy(), sqrt(np.sum(np.square(V)))]]
            while (Coords[2] - self.R) > 0:
                Coords += V * dt
                V += a * dt
                F = self.aerodynamic_power(V) + self.wind_power(V, Coords[2])
                a = a0 + F / self.m
                dots.append([t, Coords.copy(), sqrt(np.sum(np.square(V)))])
                t += dt
            for j in range(len(dots)):
                dots[j][1][0] -= dots[-1][1][0]
                dots[j][1][1] -= dots[-1][1][1]
            history.append(dots)
            print('angle =', round(angle, 3), ', X =', -Coords[0], ', Z =', -Coords[1])
        self.save_history(history, angles)
    
    def save_history(self, history, angles):
        for i in range(len(angles)):
            file = open('log_angle_' + str(round(angles[i], 2)) + '.txt', 'w')
            for j in range(len(history[i])):
                file.write('t = ' + str(round(history[i][j][0], 3)) + ' , X = ' + str(round(history[i][j][1][0], 2)) + ', Z = ' + 
                           str(round(history[i][j][1][1], 2)) + ' , Y = ' + str(round(history[i][j][1][2], 2)) + 
                           ' , V = ' + str(round(history[i][j][2], 2)) + '\n')
            file.close()
    def aerodynamic_power(self, V):
        v = sqrt(np.square(V).sum())
        f = self.f_interp(v)
        F = - f * V / v 
        return F
    
    def wind_power(self, V, h): ##!!!!!!!!
        V_relative = np.zeros(3)
        V_relative[0] = V[0] - self.Wx_interp(h)
        V_relative[1] = V[1] - self.Wz_interp(h)
        v = sqrt(np.square(V_relative).sum())
        if v == 0:
            return np.zeros(3)
        f = self.k * v ** 2
        F = - f * V_relative / v
        return F

obj = simulation(H0=H0, V0=V0, m=m, F=F, Wind=Wind, R=R)
obj.emulate(angles=angles, dt=dt)