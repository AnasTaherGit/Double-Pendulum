import numpy as np
import matplotlib.pyplot as plt
from math import *
from tkinter import *


#   Définition des Constantes -----------------
solution = 'v'  # e for euler and v for verlet
simulation = 1  # if set to 1 it will make an animated simulation
g = 9.80665
m1, m2 = 1, 1
l1, l2 = 1, 1
TRAJECTORY = False  # Si True affiche la trajectoire
ENERGY = False  # Si True affiche la representation graphique de l'energie au cours du temps
ANGLE = False  # Si True affiche la representation graphique des angles au cours du temps
PHASE = False

# N.B : Chaque plot doit être fermé pour que le prochain s'ouvre

theta = [(pi / 2, pi / 2), (pi / 6, pi / 6), (10**-8, sqrt(2) * 10**-8), (10**-8, 10**-8), (5 * 10**-3, 5 * 10**-3), (1, 1), (1.2, 1.2), (1.7, 1.7), (1.9, 1.9)]

# La constante c represente les conditions initiales suivant les valeurs suivantes
# 0 for Pi/2
# 1 for Pi/6
# 2 for 10**-8, sqrt(2)*10**-8
# 3 for 10**-8
# 4 for 5*10**-3
# 5 for 1
# 6 for 1.2
# 7 for 1.7
# 8 for 1.9

c = 8

y0 = np.array([theta[c][0], theta[c][1], 0, 0])
P = 1000
T = 20

# ------------------------------------------


def f(y, l1, l2, m1, m2):

    inter_11 = -m2 * l1 * (y[3]**2) * sin(y[0] - y[1]) * cos(y[0] - y[1]) + g * m2 * sin(y[1]) * cos(y[0] - y[1])

    inter_12 = l1 * (m1 + m2) - m2 * l1 * (cos(y[0] - y[1])**2)

    inter_1 = inter_11 / inter_12

    inter_21 = -m2 * l2 * (y[3]**2) * sin(y[0] - y[1]) - (m1 + m2) * g * sin(y[0])
    inter_22 = l1 * (m1 + m2) - m2 * l1 * (cos(y[0] - y[1])**2)

    inter_2 = inter_21 / inter_22

    s1 = inter_1 + inter_2

    inter_31 = m2 * l2 * (y[3]**2) * sin(y[0] - y[1]) * cos(y[0] - y[1]) + g * (m1 + m2) * sin(y[0]) * cos(y[0] - y[1])
    inter_32 = l2 * (m1 + m2) - m2 * l2 * (cos(y[0] - y[1])**2)

    inter_3 = inter_31 / inter_32

    inter_41 = (m1 + m2) * l1 * (y[3]**2) * sin(y[0] - y[1]) - (m1 + m2) * g * sin(y[1])
    inter_42 = l2 * (m1 + m2) - m2 * l2 * (cos(y[0] - y[1])**2)

    inter_4 = inter_41 / inter_42

    s2 = inter_3 + inter_4

    return np.array([y[2], y[3], s1, s2])


def stepEuler(y, h, l1, l2, m1, m2):

    w = y + h * f(y, l1, l2, m1, m2)

    return w


def trajectoire(y, l1, l2, m1, m2, P, T, sol='v'):
    w = np.array(y)
    print("y1 = {0:.5f}  y2 = {1:.5f}  y3 = {2:.5f}  y4 = {3:.5f}".format(w[0], w[1], w[2], w[3]))
    path = []
    h = T / P
    for i in range(P):
        if sol == 'e':
            w = np.array(stepEuler(w, h, l1, l2, m1, m2))
        else:
            w = np.array(stepVerlet(w, h, l1, l2, m1, m2))
        path.append(w)
    return np.array(path)


def plot_angles(TH, L):
    plt.title(r'$\theta$1 = {0:.8f}     $\theta$2 = {1:.8f} '.format(theta[c][0], theta[c][1]))
    plt.plot(TH, L[:, 0], 'r')
    plt.plot(TH, L[:, 1], 'b')
    plt.show()


def evol_vitesse_angulaires(TH, L):
    plt.title(r'Vitesse angulaire pour  : $\theta$1 = {0:.8f}     $\theta$2 = {1:.8f} '.format(theta[c][0], theta[c][1]))
    plt.plot(TH, L[:, 2], 'r')
    plt.plot(TH, L[:, 3], 'b')
    plt.show()


def plot_trajectoire(x1, y1, x2, y2):
    plt.title(r'$\theta$1 = {0:.8f}     $\theta$2 = {1:.8f} '.format(theta[c][0], theta[c][1]))
    plt.plot(x1, y1, 'r')
    plt.plot(x2, y2, 'b')
    plt.show()


def energie(y):

    # Cette fonction calcul l'energie pour chaque vecteur y càd à chaque instant tp

    Ec = (0.5 * (m1 + m2) * (y[2]**2) * l1) + (0.5 * m2 * (y[3]**2) * l2) + (m2 * y[2] * y[3] * l1 * l2 * cos(y[0] - y[1]))
    Ep = -(m1 + m2) * g * l1 * cos(y[0]) - m2 * l2 * g * cos(y[1])

    E = Ec + Ep

    return E

#   Calcule de l'energie


def energie_m(L):

    # Cette fonction  renvoie une vecteur contenant tout les valeurs de l'energie a chaque instant tp

    array = []
    for i in range(P):
        array.append(energie(L[i]))
    return np.array(array)


def plot_energy(TH, EH):
    plt.title(r'$\theta$1 = {0:.4f}     $\theta$2 = {1:.4f} '.format(theta[c][0], theta[c][1]))
    plt.ylim((np.min(EH) - 10, np.max(EH) + 10))
    plt.plot(TH, EH)
    plt.show()


# On remarque que l'énergie augmente or le pendule subit des forces conservatifs donc l'explosion de l'energie mécanique ...
# est du au accumulation des erreurs donc l'algorithme d'euler est inapropririé pour simuler le mouvement d'un pendule double

def stepVerlet(y, h, l1, l2, m1, m2):

    v = f(y, l1, l2, m1, m2)

    y1 = y[0] + v[0] * h + (h**2) * v[2] * 0.5
    y2 = y[1] + v[1] * h + (h**2) * v[3] * 0.5

    w = np.array([y1, y2, y[2], y[3]])

    inter = f(w, l1, l2, m1, m2)

    w1 = y[2] + (h / 2) * (inter[2] + v[2])
    w2 = y[3] + (h / 2) * (inter[3] + v[3])
    print("y1 = {0:.5f}  y2 = {1:.5f}  y3 = {2:.5f}  y4 = {3:.5f} s1={4:.5f} s2 ={5:.5f}".format(y1, y2, w1, w2, v[2], v[3]))
    return np.array([y1, y2, w1, w2])


#


# Partie 5 - Mouvement Chaotique et sensibilté au conditions initiales

# Question (1)

def plot_espaces_phases(y):

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot(y[:, 0], y[:, 2], color='b')
    ax1.set_title('Espace des phases 1')
    ax2.plot(y[:, 1], y[:, 3], color='orange')
    ax2.set_title('Espace des phases 2')

    plt.show()

#  Question (2) : Régime des petites oscillations

# (a) la forme de l'attracteur est une ellipse

# (b) la forme de l'attracteur est devenu très complexe malgré un petit changement dans les conditions initiales donc le systémes et très sensible

# au conditions initiales

# Question 3 : Transition vers le mouvement chaotique

# (a) Dans le régime des petite oscillations le pendule double est equivalent a un pendule simple qui n'a pas de caractère chaotique et un attracteur elliptique

#  mais au fur et à mesure qu'on augmente alpha on constate un mouvement complètement différent du pendule et qui est très sensible

#  au conditions initiales


def main(simulation, c=0):

    L = trajectoire(y0, l1, l2, m1, m2, P, T, sol=solution)
    # for x in L[:, 0]:
    #    print('{0:.8f}'.format(x))
    x1 = l1 * np.sin(L[:, 0])
    y1 = -l1 * np.cos(L[:, 0])

    x2 = x1 + l2 * np.sin(L[:, 1])
    y2 = y1 - l2 * np.cos(L[:, 1])

    if simulation:
        p = Simulation(x1, y1, x2, y2, T, P)
        p.mainloop()

    else:

        TH = np.linspace(0, T, P)

        if ANGLE:  # Representation graphique de l'evolution des angles du double pendule

            plot_angles(TH, L)
            evol_vitesse_angulaires(TH, L)

        #   Representation graphique de la trajectoire du double Pendule

        if TRAJECTORY:

            plot_trajectoire(x1, y1, x2, y2)

        # Represenation graphique de l'evolution de l'énergie
        if ENERGY:

            EH = energie_m(L)

            plot_energy(TH, EH)

        # Affichage de l'espace de phases
        if PHASE:

            plot_espaces_phases(L)


class Simulation(Tk):

    def __init__(self, x1, y1, x2, y2, T, P):

        Tk.__init__(self)
        self.geometry('600x600+400+40')
        self.title('Pendule')
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.T = T
        self.P = P
        self.i = 1

        # ----- Initialisation du canvas ------------
        self.field = Canvas(self, width=600, height=600, background='white')
        x01 = int(self.x1[0] * 100 + 300)
        y01 = abs(int(self.y1[0] * 100 + 300))
        x02 = int(self.x2[0] * 100 + 300)
        y02 = abs(int(self.y2[0] * 100 + 300))
        # print(x01, y01, x02, y02)
        self.line1 = self.field.create_line(300, 300, x01, y01, fill='red', width=5)
        self.line2 = self.field.create_line(x01, y01, x02, y02, fill='blue', width=5)
        self.field.pack()

        self.update()

    def update(self):
        h = self.T / self.P
        i = self.i
        x1 = int(self.x1[i] * 100 + 300)
        y1 = -int(self.y1[i] * 100 - 300)
        x2 = int(self.x2[i] * 100 + 300)
        y2 = -int(self.y2[i] * 100 - 300)
        self.field.coords(self.line1, 300, 300, x1, y1)
        self.field.coords(self.line2, x1, y1, x2, y2)
        # sleep(h)

        self.i += 1
        # print(self.i)
        if self.i < P:
            self.after(int(h * 10**3), self.update)
        else:
            pass


main(simulation, c)
