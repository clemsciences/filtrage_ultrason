# -*- coding: utf-8 -*-
import pickle

from conversion_mesure_etat import *
from math import sqrt
import collections
import scipy.linalg
from scipy.stats import norm
import generateur_chemin
from structures import Point, Velocity
___author__ = "spraakforskaren"

"""
class Simulator:
    def __init__(self, sizeTableX, sizeTableY):
        self.convertisseur = Convertisseur()

    def add_noise_to_real_position(self, x_reel, y_reel, var):
        x_reel +=  np.random.randn()*var
        y_reel += np.random.randn()*var
        m1, m2, m3 = self.convertisseur.obtenir_mesures(x_reel, y_reel)
        return np.matrix([m1, m2, m3])

    def add_noise_to_sensor_measures(self, x_reel, y_reel, var):
        m1, m2, m3 = self.convertisseur.obtenir_mesures(x_reel, y_reel)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        return np.matrix([m1, m2, m3])

    def simuler1(self, x_reel, y_reel, var):

        # Cette simulation prend la position réelle, bruite les mesures et les convertit pour obtenir les positions possibles

        m1, m2, m3 = self.add_noise_to_sensor_measures(x_reel, y_reel, var)
        res = []
        res.append(self.convertisseur.equation1_1_2(m1, m2))
        res.append(self.convertisseur.equation2_1_2(m1, m2))
        res.append(self.convertisseur.equation1_2_3(m2, m3))
        res.append(self.convertisseur.equation2_2_3(m2, m3))
        res.append(self.convertisseur.equation1_3_1(m3, m1))
        res.append(self.convertisseur.equation2_3_1(m3, m1))
        return res

    def simuler2(self, x_reel, y_reel, var):

        # Cette simulation prend la position réelle, bruite les mesures et les convertit pour obtenir les positions possibles

        m1, m2, m3 = self.convertisseur.obtenir_mesures(x_reel, y_reel)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        res = []
        res.append(self.convertisseur.equation1_1_2(m1, m2))
        res.append(self.convertisseur.equation2_1_2(m1, m2))
        res.append(self.convertisseur.equation1_2_3(m2, m3))
        res.append(self.convertisseur.equation2_2_3(m2, m3))
        res.append(self.convertisseur.equation1_3_1(m3, m1))
        res.append(self.convertisseur.equation2_3_1(m3, m1))
        return res
    def simuler2(self, x_reel, y_reel, var):

        #Cette simulation prend la position réelle, convertit les mesures, et bruite les positions

        m1, m2, m3 = self.convertisseur.obtenir_mesures(x_reel, y_reel)
        res = []
        res.append(self.convertisseur.equation1_1_2(m1, m2))
        res.append(self.convertisseur.equation2_1_2(m1, m2))
        res.append(self.convertisseur.equation1_2_3(m2, m3))
        res.append(self.convertisseur.equation2_2_3(m2, m3))
        res.append(self.convertisseur.equation1_3_1(m3, m1))
        res.append(self.convertisseur.equation2_3_1(m3, m1))
        return [(conv[0]+ np.dot(np.random.randn(),var), conv[1] + np.dot(np.random.randn(),var)) for conv in res]

    def __compute_velocity(self, x1, y1, x2, y2, duree):
        return np.matrix([(x2-x1)/duree, (y2-y1)/duree])

    def main(self, positions, duree):
        for posi in positions:
            #Vec2 p_bruit = laser.position_balise(balise.id);
            #prévoir le cas où il manque une donnée
            filtrage_lu = FiltrageLaserUnscented(posi)

            #vitesse = self.__compute_velocity(x1, y1, x2, y2, duree)
            # Mise à jour du modèle de filtrage
            vecteur_mesure = np.concatenate([posi, ])
            filtrage_lu.update(vecteur_mesure)
            # Récupération des valeurs filtrées
            p_filtre = filtrage_lu.position()

            #Vérification si l'obstacle est sur la table
            if p_filtre.x > -self.sizeTableX/2 and p_filtre.y > 0 and p_filtre.x < self.sizeTableX/2 and p_filtre.y < self.sizeTableY:
                print "sur la table"
"""

# norm.cpf
class IvanKalman:
    """
       Cette classe implémente l'algorithme du Filitrage de Kalman Unscented
    Source : Machine Learning a probabilistic perspective p. 651-652
    """
    def __init__(self, mu0, SIGMA0, si, dt, phi, d=2):
        """
        g est la fonction qui associe la position en t à la position en t+1
        h est la fonction qui associe la position réelle aux mesures correspondantes
        mu0 est la position initiale du robot
        SIGMA est la matrice de covariance initiale

        """
        self.dt = dt
        self.si = si
        self.phi = phi
        self.d = d # c'est la dimension de l'espace cachée
        # print "mu0", mu0
        self.mu = mu0  # position initialle
        self.sigma = SIGMA0  #variance initiale
        self.gamma = np.sqrt(1-phi**2)  # à nouveau bizarre !
        # print 1./(1+2*self.d)

    def _first_step(self, mesure):
        self.d = 2
        self.mu_entre = self.phi * self.mu
        self.sigma_entre = self.phi**2*self.sigma+self.gamma ** 2
        racine_sigma = np.asmatrix(scipy.linalg.sqrtm(self.sigma_entre))

        print "racine carrée", racine_sigma
        # print "mu", self.mu, "racine_sigma", racine_sigma.shape, racine_sigma[:, 1].shape
        # print "self.mu - self.gamma*racine_sigma[:, i]", self.mu - self.gamma*racine_sigma[:, 1]
        self.points_sigma = [self.mu_entre]
        self.points_sigma.extend([self.mu_entre - racine_sigma[:, i] for i in range(0, self.d)])  # de -1 à - self.d
        self.points_sigma.extend([self.mu_entre + racine_sigma[:, i] for i in range(0, self.d)])  # de 1 à self.d

        self.S = np.asmatrix(np.zeros((3, 3)))
        for sigma_point in self.points_sigma:
            # print sigma_point.shape
            transformation = self.hm(sigma_point)
            self.S += np.dot(transformation, transformation.T)
        print "S", self.S

        self.S *= 1./(1+2*self.d)
        self.S += self.si**2*np.asmatrix(np.eye(3))
        # print self.S.shape
        # print "sigma_point", sigma_point[0].shape, self.hm(sigma_point)).shape,
        self.K = np.asmatrix(np.zeros((2, 3)))
        for sigma_point in self.points_sigma:
            # print sigma_point.shape
            self.K += np.dot(np.dot(sigma_point, self.hm(sigma_point).T), scipy.linalg.inv(self.S))
        self.K *= 1./(1+2*self.d)
        # print self.K.shape, mesure.shape, self.hm(self.points_sigma[0]).shape,

        h_sum = np.asmatrix(np.zeros((3, 1)))
        for sigma_point in self.points_sigma:
            h_sum += self.hm(sigma_point)
        print h_sum
        print "mesures", mesure, 1./(1+2*self.d)*h_sum
        self.mu_prochain = self.mu_entre + np.dot(self.K,  (mesure - 1./(1+2*self.d)*h_sum))
        self.sigma_prochain = np.asmatrix(np.zeros((2, 2)))
        # print self.K.shape, self.hm(sigma_point).shape, self.points_sigma[0].T.shape
        for sigma_point in self.points_sigma:
            self.sigma_prochain += np.dot(np.dot(self.K, self.hm(sigma_point)), sigma_point.T)
        self.sigma_prochain = self.sigma_entre - 1./(1+2*self.d) * self.sigma_prochain
        self.mu = self.mu_prochain
        self.sigma = self.sigma_prochain

    def filter(self, y):
        """
        On sépare les deux pas parce que l'un a besoin d'une mesure tandis que l'autre non
        :param y: c'est la mesure provenant du capteur
        """

        if y is None:  # aucune mesure
            self._first_step(self.hm(self.mu))
        else:
            print "mesure !"
            self._first_step(y)

    def h(self, x):
        """
        La fonction renvoie les données comme si elles ont été mesurées.
        """
        # print "x shape", x.shape
        conv = Converter()
        measures = conv.get_measures_from_state(x[0], x[1])
        return np.asmatrix(measures).T#[:,np.newaxis]

    def m(self, gauss_x):
        # print gauss_x
        print gauss_x.shape
        print gauss_x
        # print gauss_x[0, 0]
        # print gauss_x[1, 0]
        x = -1500 + 3000 * norm.cdf(gauss_x[0, 0])
        y = 2000 * norm.cdf(gauss_x[1, 0])
        return np.matrix([[x], [y]])

    def hm(self, gauss_x):
        return self.h(self.m(gauss_x))

    def get_state(self):
        """
        self.mu est la position moyenne
        """
        return self.mu


class IvanFilter:

    def __init__(self, x0, dt):
        """
        x0 est un array(x,y) ou array(x,y,x point, y point)
        """
        self.dt = dt
        x0 = np.array(self.real_position_to_abstract(x0))
        # x0 = np.array([0, 0])
        mu0 = np.asmatrix(x0)
        mu0 = mu0.T
        # print "mu", mu0

        d = 2
        #x = np.array([1400,100,0.,0.])[:, np.newaxis] # vecteur d'état au départ
        SIGMA0 = np.matrix([[1, 0.], [0., 1]])
        si = 0.01
        phi = 0.99
        self.ukf = IvanKalman(mu0, SIGMA0, si, dt, phi, d=d)
        self.historique = collections.deque(maxlen=3)
        self.valeurs_rejetees = 0
        self.acceleration = None

    def real_position_to_abstract(self, p0):
        """

        :param p0:
        :return:
        """
        return [norm.ppf((p0[0] + 1500) / 3000.), norm.ppf(p0[1] / 2000.)]

    def abstract_position_to_real(self, x0):
        """

        :param x0:
        :return:
        """
        return [norm.cdf(x0[0])*3000 - 1500, norm.cdf(x0[1]) * 2000.]

    def get_state(self):
        """

        :return: the state of the robot
        """
        return self.ukf.mu

    def get_state_position(self):
        """

        :return: a Point
        """
        state = self.ukf.mu
        # x, y = state[0, 0], state[1, 0]
        x, y = self.abstract_position_to_real([state[0, 0], state[1, 0]])
        # print state
        return Point(float(x), float(y))

    def get_state_velocity(self):
        """

        :return: a Velocity
        """
        state = self.ukf.mu
        return Velocity(float(state[2]), float(state[3]))

    def update(self, y):
        """
        Je ne sait pas si c'est vraiement utilse
        fonction qui est utilisé à chaque mesure
        y est un vecteur de mesure de dimension 4 : (x, y, x point, y point)
        """

        #if self._filtrage_acceleration(Point(x, y)):
        #    self.last_point = Point(x, y)
        self.ukf.filter(y)
        pos_filtre = self.ukf.get_state()
        #    self.filtre_kalman.measurement(np.array([x,y])[:, np.newaxis])
        self.historique.append(self.ukf.get_state())
        # print y, pos_filtre[0], pos_filtre[1]
        #else:
        #    self.last_point = None
        return pos_filtre
        #    self.filtre_kalman.prediction()

    def _acceleration_filtering(self, pointm0):
        """
        Vérifie si le point est cohérent avec la position actuelle, en basant sur l'accélération
        """
        # Pas suffisamment de valeurs précédentes pour calculer l'accélération
        if len(self.historique) != 3:
            return True

        # 3 derniers points valides
        pointm1 = self.historique[2]
        pointm2 = self.historique[1]
        pointm3 = self.historique[0]

        # Vecteurs vitesses et accélération
        vitesse_actuelle = pointm0 - pointm1
        vitesse_m1 = pointm1 - pointm2
        vitesse_m2 = pointm2 - pointm3
        acceleration_actuelle = vitesse_actuelle - vitesse_m1
        acceleration_precedente = vitesse_m1 - vitesse_m2
        jerk = acceleration_actuelle - acceleration_precedente

        # Produit scalaire pour savoir s'il y a accélération ou décélération
        produit = acceleration_actuelle.x * vitesse_m1.x + acceleration_actuelle.y * vitesse_m1.y

        # Rejette les accélérations brutales
        if acceleration_actuelle.norme() / self.dt**2 > 50000 and self.valeurs_rejetees < 3:
            #~ print("accélération = {0}, produit = {1}, jerk = {2}".format(acceleration_actuelle.norme() / self.dt**2, produit, jerk.norme() / self.dt**3))
            self.valeurs_rejetees += 1
            return False
        else:
            self.valeurs_rejetees = 0
            return True


class UnscentedKalman:
    """
       Cette classe implémente l'algorithme du Filitrage de Kalman Unscented
    Source : Machine Learning a probabilistic perspective p. 651-652
    """
    def __init__(self, mu0, SIGMA0, Q, R, d, alpha, beta, kappa, dt):
        """
        g est la fonction qui associe la position en t à la position en t+1
        h est la fonction qui associe la position réelle aux mesures correspondantes
        mu0 est la position initiale du robot
        SIGMA est la matrice de covariance initiale
        
        """
        self.dt = dt
        self.d = d  # c'est la dimension de l'espace cachée
        self.alpha = float(alpha) # paramètre bizarre
        self.beta = float(beta)  # paramètre bizarre
        self.kappa = kappa  # paramètre bizarre
        self.mu = mu0  # position initiale
        self.SIGMA = SIGMA0  #variance initiale
        self.lam = alpha**2*(d+kappa) - d  # un paramètre un peu étrange aussi
        self.gamma = d+self.lam  # à nouveau bizarre !
        self.Q = Q # matrice de covariance de la gaussienne de l'équation d'évolution
        self.R = R  # matrice de covariance de la gaussienne de l'équation d'observation

        #calculs de coefficients
        self.wm0 = float(self.lam) /(self.d+self.lam)
        self.wm = 1/(2.*(self.d+self.lam))  # wm = wc (wc non implémenté)
        self.wc0 = self.wm0+(1-self.alpha**2+self.beta)

    def _first_step(self):
        """
        Estimation de mu_barre et de sigma_barre à partir de l'itération précédente
        """
        # first Unscented transform
        racine_sigma = np.asmatrix(scipy.linalg.sqrtm(self.gamma*self.SIGMA))
        # self.mu = self.mu.T
        # print self.mu
        # print racine_sigma[:, 0]
        # print "racine carrée", racine_sigma
        # print "mu", self.mu, "racine_sigma", racine_sigma.shape, racine_sigma[:, 1].shape
        # print "self.mu - self.gamma*racine_sigma[:, i]", self.mu - self.gamma*racine_sigma[:, 1]
        self.points_sigma = [self.mu]
        self.points_sigma.extend([self.mu - racine_sigma[:, i] for i in range(0, self.d)])  # de -1 à - self.d
        self.points_sigma.extend([self.mu + racine_sigma[:, i] for i in range(0, self.d)])  # de 1 à self.d

        # print "self.points_sigma", len(self.points_sigma), np.array(self.points_sigma).shape
        # self.points_sigma = np.array(self.points_sigma)
        # self.z_etoile_barre = np.array([self.g(self.points_sigma[i]).reshape((4,)) for i in range(0, 2*self.d)])

        self.z_etoile_barre = []  # np.asmatrix(np.zeros(self.points_sigma.shape))
        # print len(self.points_sigma)
        for i in range(len(self.points_sigma)):
            self.z_etoile_barre.append(self.g(self.points_sigma[i]))
            # self.z_etoile_barre[i,:] = self.g(self.points_sigma[i])

        self.mu_barre = self.wm0*self.z_etoile_barre[0]

        # print "mu_barre", self.mu_barre.shape
        # print self.mu_barre
        # print self.z_etoile_barre[2]
        # print self.wm
        for i in range(1, len(self.z_etoile_barre)):
            # print i
            # print "self.wm*self.z_etoile_barre[i]", (self.wm*self.z_etoile_barre[i]).shape
            self.mu_barre += self.wm*self.z_etoile_barre[i]

        # calcul de SIGMA_barre
        # print "mu barre", self.mu_barre.shape
        # print "z_etoile_barre", len(self.z_etoile_barre), len(self.z_etoile_barre[0]), self.z_etoile_barre[0], self.mu_barre
        self.SIGMA_barre = self.wc0*np.dot((self.z_etoile_barre[0] - self.mu_barre),
                                           (self.z_etoile_barre[0] - self.mu_barre).T)
        for i in range(1, len(self.z_etoile_barre)):
            self.SIGMA_barre += self.wm*np.dot((self.z_etoile_barre[i] - self.mu_barre),
                                               (self.z_etoile_barre[i] - self.mu_barre).T)
        # print "SIGMA_barre", self.SIGMA_barre.shape
        self.SIGMA_barre += self.Q

    def _second_step(self, y):
        """
        y est la mesure à l'instant t
        Estimation de mu et de sigma présent
        """
        #second unscented transform
        racine_sigma_barre = np.asmatrix(scipy.linalg.sqrtm(self.SIGMA_barre))
        self.points_sigma_barre = [self.mu_barre] + \
                                  [self.mu_barre-sqrt(self.gamma)*racine_sigma_barre[:, i] for i in range(0, self.d)] +\
                                  [self.mu_barre+sqrt(self.gamma)*racine_sigma_barre[:, i] for i in range(0, self.d)]
        # print "points_sigma_barre", self.points_sigma_barre
        self.y_etoile_barre = [self.h(self.points_sigma_barre[i]) for i in range(0, len(self.points_sigma_barre))]

        self.y_chapeau = self.wm0*self.y_etoile_barre[0]
        for i in range(1, len(self.y_etoile_barre)):
            self.y_chapeau += self.wm*self.y_etoile_barre[i]
        #Calcul de S
        self.S = self.wc0*np.dot(self.y_etoile_barre[0] - self.y_chapeau,
                                 (self.y_etoile_barre[0] - self.y_chapeau).T)
        for i in range(1, len(self.y_etoile_barre)):
            self.S += self.wm*np.dot(self.y_etoile_barre[i] - self.y_chapeau,
                                     (self.y_etoile_barre[i] - self.y_chapeau).T)
        # print "S", self.S.shape
        self.S += self.R
        #calcul de SIGMA_z_y
        self.SIGMA_z_y_barre = self.wc0*np.dot((self.z_etoile_barre[0] - self.mu_barre),
                                          (self.y_etoile_barre[0] - self.y_chapeau).T)
        for i in range(1, len(self.z_etoile_barre)):
            self.SIGMA_z_y_barre += self.wm*np.dot(self.z_etoile_barre[i] - self.mu_barre,
                                              (self.y_etoile_barre[i] - self.y_chapeau).T)
        self.K = np.dot(self.SIGMA_z_y_barre, scipy.linalg.inv(self.S))
        #print "K", self.K.shape, "y", y.shape
        #Les valeurs qui nous interessent
        self.mu = self.mu_barre + np.dot(self.K, y - self.y_chapeau)
        self.SIGMA = self.SIGMA_barre - np.dot(np.dot(self.K, self.S), self.K.T)

    def filter(self, y):
        """
        On sépare les deux pas parce que l'un a besoin d'une mesure tandis que l'autre non
        :param y: c'est la mesure provenant du capteur
        """
        self._first_step()
        if y is None:  # aucune mesure
            self._second_step(self.mu)  # pour l'instant, c'est une proposition, on pourra trouver mieux
        else:
            self._second_step(y)

    def g(self, x, dim=2):
        """
        Pour notre modélisation, on choisit comme vecteur x = [abscisse de la position,
        ordonnée de la position] ou [abscisse de la position,
        ordonnée de la position, abscisse de la vitesse, ordonnées de la vitesse]
        """
        # if dim == 2:
        #     F = np.matrix([[1., 0.], [0., 1.]])
        # elif dim == 3:
        #     F = np.matrix([[1., 0, 0], [0., 1., 0.], [0., 0., 1.]])
        # else: #dim == 4
        #     F = np.matrix([[1., 0., self.dt, 0.], [0., 1., 0., self.dt], [0., 0., 1., 0.], [0., 0., 0., 1.]])
        # print F.shape, x.shape
        F = np.matrix([[1., 0.], [0., 1.]])
        res = np.dot(F, x)
        # print res.shape
        return res

    def h(self, x):
        """
        La fonction renvoie les données comme si elles ont été mesurées.
        """
        # print "x shape", x.shape
        conv = Converter()
        measures = conv.get_measures_from_state(x[0], x[1])
        return np.asmatrix(measures).T#[:,np.newaxis]

    def get_state(self):
        """
        self.mu est la position moyenne
        """
        return self.mu


# classe à arranger pour le Kalman Unscented
class UnscentedKalmanFilter:
    
    def __init__(self, x0, dt, coeff_s=100, coeff_q=0.00001, coeff_r=0.1, dime=2):
        """
        x0 est un array(x,y) ou array(x,y,x point, y point)
        """
        self.dt = dt
        #x = np.array(x).T
        mu0 = x0
        d = 2
        kappa = 2
        alpha = 1
        beta = 0


        # x = np.array([1400,100,0.,0.])[:, np.newaxis] # vecteur d'état au départ
        # if dime == 2:
        #     SIGMA0 = np.matrix([[0.001, 0.], [0., 0.001]])
        #     R = np.matrix([[90, 0., 0.], [0., 90, 0.], [0., 0., 90]])  #dimension de la matrice égale au nombre de dimensions des mesures !
        #     Q = np.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2],
        #                   [self.dt**2/2., 0, 4*self.dt, 0], [0, self.dt**2/2, 0, 4*self.dt]])
        # else:
        # SIGMA0 = np.matrix([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]]) # incertitude initiale
        SIGMA0 = np.eye(2)
        SIGMA0 *= coeff_s
        R = np.eye(3)
        # R = np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])  #dimension de la matrice égale au nombre de dimensions des mesures !
        R *= coeff_r

        #Q = np.matrix([[self.dt**3/3., self.dt**2/2., 0, 0],[self.dt**2/2.,self.dt, 0, 0],
        #            [0,0,self.dt**3/3.,self.dt**2/2],[0,0,self.dt**2/2,self.dt]])
        #Q *= 20;
        Q = np.eye(2)
        # Q = np.matrix([[1., 0, 0, 0],[0, 1, 0, 0],
        #                [0, 0, 1, 0],[0, 0, 0, 1]])
        # Q = np.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2],
        #                [self.dt**2/2., 0, 4*self.dt, 0], [0, self.dt**2/2, 0, 4*self.dt]])
        #Q = np.matrix([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 4, 0],[0, 0, 0, 4]])
        Q *= coeff_q # (1-0.98**2)

        self.ukf = UnscentedKalman(mu0, SIGMA0, Q, R, d, alpha=alpha, beta=beta, kappa=kappa, dt=dt)
        self.historique = collections.deque(maxlen=3)
        self.valeurs_rejetees = 0
        self.acceleration = None
        
    def get_state(self):
        """

        :return: the state of the robot
        """
        return self.ukf.mu
        
    def update_dt(self, new_dt):
        """
        Modifie la période d'échantillonage
        """
        self.dt = new_dt
        self.ukf.F[0,2] = new_dt
        self.ukf.F[1,3] = new_dt

    def get_state_position(self):
        """

        :return: a Point
        """
        state = self.ukf.mu
        # print state
        return Point(float(state[0]), float(state[1]))

    def get_state_velocity(self):
        """

        :return: a Velocity
        """
        state = self.ukf.mu
        return Velocity(float(state[2]), float(state[3]))
                
    def update(self, y):
        """
        Je ne sait pas si c'est vraiement utilse
        fonction qui est utilisé à chaque mesure
        y est un vecteur de mesure de dimension 4 : (x, y, x point, y point)
        """

        #if self._filtrage_acceleration(Point(x, y)):
        #    self.last_point = Point(x, y)
        self.ukf.filter(y)
        pos_filtre = self.ukf.get_state()
        #    self.filtre_kalman.measurement(np.array([x,y])[:, np.newaxis])
        self.historique.append(self.ukf.get_state())
        # print y, pos_filtre[0], pos_filtre[1]
        #else:
        #    self.last_point = None
        return pos_filtre
        #    self.filtre_kalman.prediction()
        
    def _acceleration_filtering(self, pointm0):
        """
        Vérifie si le point est cohérent avec la position actuelle, en basant sur l'accélération
        """
        # Pas suffisamment de valeurs précédentes pour calculer l'accélération
        if len(self.historique) != 3:
            return True
            
        # 3 derniers points valides
        pointm1 = self.historique[2]
        pointm2 = self.historique[1]
        pointm3 = self.historique[0]
        
        # Vecteurs vitesses et accélération
        vitesse_actuelle = pointm0 - pointm1
        vitesse_m1 = pointm1 - pointm2
        vitesse_m2 = pointm2 - pointm3
        acceleration_actuelle = vitesse_actuelle - vitesse_m1
        acceleration_precedente = vitesse_m1 - vitesse_m2
        jerk = acceleration_actuelle - acceleration_precedente
        
        # Produit scalaire pour savoir s'il y a accélération ou décélération
        produit = acceleration_actuelle.x * vitesse_m1.x + acceleration_actuelle.y * vitesse_m1.y
        
        # Rejette les accélérations brutales
        if acceleration_actuelle.norme() / self.dt**2 > 50000 and self.valeurs_rejetees < 3:
            #~ print("accélération = {0}, produit = {1}, jerk = {2}".format(acceleration_actuelle.norme() / self.dt**2, produit, jerk.norme() / self.dt**3))
            self.valeurs_rejetees += 1
            return False
        else:
            self.valeurs_rejetees = 0
            return True


class Kalman:
    """
    Classe  qui implémente le filtre de Kalman. Utilisé, il permet de filtrer une suite discrétisé de valeur mesurée
    """
  
    def __init__(self, x, P, F, H, R, Q):
        """
        :param x: Etat initial du truc à suivre (un vecteur position le plus souvent)
        :param P: matrice incertitude sur le modèle d'évolution
        :param F: matrice de transition du modèle d'évolution
        :param H: matrice matrice d'observation
        :param R: matrice de covariance de la mesure
        :param Q: matrice de covariance du modèle

        """
        self.x = x
        self.P = P
        self.F = F
        self.H = H
        self.R = R
        self.Q = Q
  
    def predict(self, u=None):
        """
        C'est la partie prédiction, le modèle imagine comment l'état suivant sera
        :param u: un vecteur "déplacement", si onsait de combien ça a dû bouger, il faut le mettre
        """
        if u is None:
            u = np.zeros(self.x.shape[0])[:, np.newaxis]
        self.x = np.dot(self.F, self.x) + u
        self.P = np.dot(np.dot(self.F, self.P), self.F.T) + self.Q

    def measure(self, mes):
        """
        C'est la partie où on prend en compte la nouvelle mesure
        :param mes: vecteur de même dimension que x
        """
        y = mes - np.dot(self.H, self.x)
        S = np.dot(np.dot(self.H, self.P), self.H.T) + self.R
        K = np.dot(np.dot(self.P, self.H.T), np.linalg.inv(S))
        self.x = self.x + np.dot(K, y)
        self.P = np.dot((np.identity(self.x.shape[0]) - np.dot(K, self.H)), self.P)

    def filter(self, mes, u=None):
        """
        Méthode qui condense une étape du filtrage
        """
        self.predict(u)
        self.measure(mes)


class ExtendedKalman:
    """
    Un jour, ça sera fait
    """
    pass
    

class FiltrageKalman:
    """
    Classe qui utilise Kalman
    """
    def __init__(self, x0, dt=0.025):
        """

        :param x0: état initial
        :param dt: pas de la mesure (période d'échantillonage)
        """
        self.dt = dt
        x = x0  #np.array([1400, 100, 0., 0.])[:, np.newaxis] # vecteur d'état au départ
        P = np.matrix([[30.0, 0., 0., 0.], [0., 30., 0., 0.], [0., 0., 10., 0.], [0., 0., 0., 10.]]) # incertitude initiale
        F = np.matrix([[1., 0., self.dt, 0.], [0., 1., 0., self.dt], [0., 0., 1., 0.], [0., 0., 0., 1.]]) # matrice de transition
        H = np.matrix([[1., 0., 0., 0.], [0., 1., 0., 0.]])# matrice d'observation
        R = np.matrix([[900, 0.], [0., 900]]) # incertitude sur la mesure

        #Q = np.matrix([[self.dt**3/3., self.dt**2/2., 0, 0],[self.dt**2/2.,self.dt, 0, 0],
        #            [0,0,self.dt**3/3.,self.dt**2/2],[0,0,self.dt**2/2,self.dt]])
        #Q *= 20;

        Q = np.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2], \
                       [self.dt**2/2., 0, 4*self.dt, 0],[0, self.dt**2/2, 0, 4*self.dt]])
        #Q = np.matrix([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 4, 0],[0, 0, 0, 4]])
        Q *= 30
        self.kalman_filter = Kalman(x, P, F, H, R, Q)
        self.history = collections.deque(maxlen=3)
        self.rejected_values = 0
        self.acceleration = None
        
    def get_current_state(self):
        return self.kalman_filter.x

    def get_current_position(self):
        state = self.get_current_state()
        return Point(float(state[0]), float(state[1]))
        
    def update_dt(self, new_dt):
        self.dt = new_dt
        self.kalman_filter.F[0,2] = new_dt
        self.kalman_filter.F[1,3] = new_dt
    
    def get_last_state(self):
        #state = self.filtre_kalman.x
        #return Point(int(state[0]), int(state[1]))
        return self.last_point
    
    # def vitesse(self):
    #     state = self.kalman_filter.x
    #     return Vitesse(int(state[2]), int(state[3]))
                
    def update(self, x, y):
        if self.acceleration_filtering(Point(x, y)):
            self.last_point = Point(x, y)
            self.kalman_filter.predict()
            self.kalman_filter.measure(np.array([x, y])[:, np.newaxis])
            self.history.append(self.get_last_state())
        else:
            self.last_point = None
            self.kalman_filter.predict()
        
    def acceleration_filtering(self, pointm0):
        """
        Vérifie si le point est cohérent avec la position actuelle, en se basant sur l'accélération
        """
        # Pas suffisamment de valeurs précédentes pour calculer l'accélération
        if len(self.history) != 3:
            return True
            
        # 3 derniers points valides
        pointm1 = self.history[2]
        pointm2 = self.history[1]
        pointm3 = self.history[0]
        
        # Vecteurs vitesses et accélération
        vitesse_actuelle = pointm0 - pointm1
        vitesse_m1 = pointm1 - pointm2
        vitesse_m2 = pointm2 - pointm3
        acceleration_actuelle = vitesse_actuelle - vitesse_m1
        acceleration_precedente = vitesse_m1 - vitesse_m2
        jerk = acceleration_actuelle - acceleration_precedente
        
        # Produit scalaire pour savoir s'il y a accélération ou décélération
        produit = acceleration_actuelle.x * vitesse_m1.x + acceleration_actuelle.y * vitesse_m1.y
        
        # Rejette les accélérations brutales
        if acceleration_actuelle.norme() / self.dt**2 > 50000 and self.rejected_values < 3:
            #~ print("accélération = {0}, produit = {1}, jerk = {2}".format(acceleration_actuelle.norme() / self.dt**2, produit, jerk.norme() / self.dt**3))
            self.rejected_values += 1
            return False
        else:
            self.rejected_values = 0
            return True


def get_velocity(positions, dt):
    """

    :param positions: matrice des positions
    :param dt: période d'échantillonage
    :return:
    """
    velocities = np.zeros(positions.shape)
    for i in range(1, velocities.shape[0]):
        velocities[i, :] = (positions[i, :] - positions[i-1, :])/float(dt)
    return velocities


def get_distance(point1, point2):
    """

    :param point1: couple (x,y)
    :param point2: couple (x,y)
    :return:
    """
    return sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)


def get_mer(real, estimated):
    """
    Renvoie la moyenne des distances
    :param real: liste de couples (x,y)
    :param estimated: liste de couples (x,y)
    :return:
    """
    try:
        for i in range(len(real)):
            print "réel ", real[i], " estimé ", estimated[i]
            print real[i].distance(estimated[i])
        res = [real[i].distance(estimated[i]) for i in range(len(real))]
        # print "La moyenne des distances entre les estimations et la réalité est : "
        return sum(res)/float(len(res))
    except ValueError:
        print real[i]
        print estimated[i]
    # except IndexError:
    #     print i
    #     print len(real)
    #     print len(estimated)


def squared_error(real_values, filtered_values):
    """
    Renvoie l'erreur quadratique moyenne !
    :param real_values:
    :param filtered_values:
    :return:
    """
    diff = real_values[:, :2] - filtered_values[:, :2]
    res = diff.T.dot(diff)
    "L'erreur quadratique est : "
    return res


def script_unscented_with_real_measures():
    """
    Script utilisant le filtre de Kalman unscented avec des mesures réelles !
    :return:
    """
    print "script_unscented_with_real_measures"
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))
    # vite = get_velocity(measures_pos, 0.025)
    # measures = np.concatenate((measures_pos, vite), axis=1)
    measures = np.asmatrix(measures_pos)
    fichier = open("mesures_simulees_10.txt", "w")
    for coeff_s in range(0, 20, 1):
        for coeff_q in range(0, 20, 1):
            for coeff_r in range(1, 20, 1):
                filtering = UnscentedKalmanFilter(measures[0, :].T, dt=0.025, coeff_s=10**-(coeff_s-10),
                                                  coeff_q=10**-(coeff_q-10), coeff_r=10**-(coeff_r-10), dime=4)
                var = 10
                l_pos_filter = []
                for i in range(1, measures.shape[0]):
                    x = measures[i, 0]
                    y = measures[i, 1]
                    # print "x et y", x, y
                    conv = Converter()
                    m1, m2, m3 = conv.get_measures_from_state(x, y)
                    m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
                    filtering.update(np.asmatrix([m1, m2, m3]).T)
                    pos = filtering.get_state_position()
                    l_pos_filter.append(pos)
                fichier.write(str(get_mer(real_path, l_pos_filter))+"\t"+str(10**-(coeff_s-10))+"\t"+str(10**-(coeff_q-10))+"\t"+str(10**-(coeff_r-10))+"\n")
    fichier.close()


def lire_fic_opti():
    valeurs = []
    coeff_s = []
    coeff_q = []
    coeff_r = []
    with open("mesures_simulees_10.txt", "r") as f:
        res = f.read()
    for r in res.split("\n"):
        print r
        if len(r) != 0:
            ligne = r.split("\t")
            valeurs.append(ligne[0])
            coeff_s.append(ligne[1])
            coeff_q.append(ligne[2])
            coeff_r.append(ligne[3])
    maxi = max(valeurs)
    print maxi
    ind = valeurs.index(maxi)
    print ind
    print "Q", coeff_q[ind]
    print "R", coeff_r[ind]
    print "S", coeff_s[ind]


def script_classic_trajectory_with_real_measures():
    """
    Script utilisant le filtre de Kalman  avec des mesures réelles !
    :return:
    """
    print "script_classic_trajectory_with_real_measures"
    dt=0.025
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))

    l_pos_filtre = [real_path[0]]
    vite = get_velocity(measures_pos, dt)
    measures = np.concatenate((measures_pos, vite), axis=1)
    measures = np.asmatrix(measures)
    filtering = FiltrageKalman(measures[0, :].T, dt=dt)
    var = 10
    for i in range(1, measures.shape[0]):
        x = measures[i, 0]
        y = measures[i, 1]
        # print "x et y", x, y, "i" ,i
        x_bruite, y_bruite = x + np.random.randn()*var, y + np.random.randn()*var
        filtering.update(x_bruite, y_bruite)
        pos = filtering.get_current_position()
        l_pos_filtre.append(pos)
    # print erreur_quadratique(measures, np.asmatrix(np.array(l_pos_filtre)))
    print get_mer(real_path, l_pos_filtre)


def script_classic_trajectory():
    """
    Script utilisant le filtre de Kalman  avec des mesures simulées à partir d'une trajectoire inventée !
    :return:
    """
    print "script_classic_trajectory"
    l_points = [[-1000., 200.], [-1000., 800.], [-400., 1200.], [500., 500.], [1100., 180.]]
    dt = 0.025
    real_path = generateur_chemin.generate_path(l_points=l_points, velocity_translation=25,
                                                            velocity_rotation=0.7, dt=dt)
    measures_pos = np.array(real_path)
    real_path_point = []
    for couple in real_path:
        x, y = couple
        pos = Point(x, y)
        real_path_point.append(pos)

    l_pos_filtre = [real_path_point[0]]
    vite = get_velocity(measures_pos, dt)
    measures = np.concatenate((measures_pos, vite), axis=1)
    measures = np.asmatrix(measures)
    filtering = FiltrageKalman(measures[0, :].T, dt=dt)
    var = 10
    for i in range(1, measures.shape[0]):
        x = measures[i, 0]
        y = measures[i, 1]
        # print "x et y", x, y, "i" ,i
        x_bruite, y_bruite = x + np.random.randn()*var, y + np.random.randn()*var
        filtering.update(x_bruite, y_bruite)
        pos = filtering.get_current_position()
        l_pos_filtre.append(pos)
    print get_mer(real_path_point, l_pos_filtre)


def script_unscented_trajectory_with_file():
    """
    Script utilisant le filtre de Kalman unscented avec des mesures simulées à partir d'une trajectoire inventée !
    :return:
    """
    print "script_unscented_trajectory"
    l_points = [[-1000., 200.], [-1000., 800.], [-400., 1200.], [500., 500.], [1100., 180.]]
    dt = 0.025
    real_path = generateur_chemin.generate_path(l_points=l_points, velocity_translation=1,
                                                            velocity_rotation=0.5, dt=dt)
    measures_pos = np.array(real_path)
    real_path_point = []
    for couple in real_path:
        x, y = couple
        pos = Point(x, y)
        real_path_point.append(pos)
    l_pos_filter = [real_path_point[0]]
    # vite = get_velocity(measures_pos, dt)
    # measures = np.concatenate((measures_pos, vite), axis=1)
    measures = np.asmatrix(measures_pos)
    mesures_us = []
    filtering = UnscentedKalmanFilter(measures[0, :].T, dt=dt, dime=2)
    var = 0
    # with open("mesures_chemin_determinites", "rb") as f:
    #     m = pickle.Unpickler(f)
    #     mesures_us = m.load()
    for i in range(1, measures.shape[0]):
        # print mesures_us[i-1]

        x = measures[i, 0]
        y = measures[i, 1]
        conv = Converter()
        m1, m2, m3 = conv.get_measures_from_state(x, y)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        mesures_us.append([m1, m2, m3])
        # m1, m2, m3 = mesures_us[i-1]
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        l_pos_filter.append(filtering.get_state_position())
    with open("mesures_chemin_determinites", "wb") as f:
        m = pickle.Pickler(f)
        m.dump(mesures_us)
    print get_mer(real_path_point, l_pos_filter)


def script_unscented_trajectory():
    """
    Script utilisant le filtre de Kalman unscented avec des mesures simulées à partir d'une trajectoire inventée !
    :return:
    """
    print "script_unscented_trajectory"
    l_points = [[-1000., 200.], [-1000., 800.], [-400., 1200.], [500., 500.], [1100., 180.]]
    dt = 0.025
    real_path = generateur_chemin.generate_path(l_points=l_points, velocity_translation=50,
                                                            velocity_rotation=0.5, dt=dt)
    measures_pos = np.array(real_path)
    real_path_point = []
    for couple in real_path:
        x, y = couple
        pos = Point(x, y)
        real_path_point.append(pos)
    l_pos_filter = [real_path_point[0]]
    # vite = get_velocity(measures_pos, dt)
    # measures = np.concatenate((measures_pos, vite), axis=1)
    measures = np.asmatrix(measures_pos)
    filtering = UnscentedKalmanFilter(measures[0, :].T, dt=dt, dime=4)
    var = 10
    # mesures_us = []
    for i in range(1, measures.shape[0]):
        x = measures[i, 0]
        y = measures[i, 1]
        conv = Converter()
        m1, m2, m3 = conv.get_measures_from_state(x, y)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        # mesures_us.append([m1, m2, m3])
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        l_pos_filter.append(filtering.get_state_position())
    # with open("mesures_chemin_determinites", "wb") as f:
    #     m = pickle.Pickler(f)
    #     m.dump(mesures_us)
    print get_mer(real_path_point, l_pos_filter)


def script_unscented_with_fixed_trajectory_only():
    print "script_unscented_with_fake_ultrasound_measures"
    dt = 0.025
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))
    l_pos_filter = [Point(measures_pos[0, :][0], measures_pos[0, :][1])]
    # vite = get_velocity(measures_pos, dt)
    # measures_pos = np.concatenate((measures_pos, vite), axis=1)
    measures_pos = np.asmatrix(measures_pos)
    filtering = UnscentedKalmanFilter(measures_pos[0, :].T, dt=dt)
    var = 10
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        conv = Converter()
        m1, m2, m3 = conv.get_measures_from_state(x, y)
        print m1, m2, m3
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        pos = filtering.get_state_position()
        l_pos_filter.append(pos)
    print get_mer(real_path, l_pos_filter)


def script_unscented_with_fake_ultrasound_measures():
    print "script_unscented_with_fake_ultrasound_measures"
    dt = 0.025
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))
    l_pos_filter = [Point(measures_pos[0, :][0], measures_pos[0, :][1])]
    # vite = get_velocity(measures_pos, dt)
    # measures_pos = np.concatenate((measures_pos, vite), axis=1)
    measures_pos = np.asmatrix(measures_pos)
    measures_us = np.genfromtxt("mesures_25_sigma10.txt", delimiter=",")
    filtering = UnscentedKalmanFilter(measures_pos[0, :].T, dt=dt)
    # var = 10
    for i in range(1, measures_us.shape[0]):
        m1 = measures_us[i, 0]
        m2 = measures_us[i, 1]
        m3 = measures_us[i, 2]
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        pos = filtering.get_state_position()
        l_pos_filter.append(pos)
    print get_mer(real_path, l_pos_filter)


def script_unscented_ivan():
    dt = 0.025
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))
    l_pos_filter = [Point(measures_pos[0, :][0], measures_pos[0, :][1])]
    # vite = get_velocity(measures_pos, dt)
    # measures_pos = np.concatenate((measures_pos, vite), axis=1)
    measures_pos = np.asmatrix(measures_pos)
    measures_us = np.genfromtxt("mesures_25_sigma10.txt", delimiter=",")
    print measures_pos[0, :].T
    filtering = IvanFilter(measures_pos[0, :].T, dt=dt)
    # var = 10
    for i in range(1, measures_us.shape[0]):
        m1 = measures_us[i, 0]
        m2 = measures_us[i, 1]
        m3 = measures_us[i, 2]
        filtering.update(np.asmatrix([m3, m2, m1]).T)
        pos = filtering.get_state_position()
        print "pos", pos
        l_pos_filter.append(pos)
    print get_mer(real_path, l_pos_filter)


def script_random_unscented_ivan():
    dt = 0.025
    measures_pos = np.genfromtxt("mesures_25.txt", delimiter="\t")
    real_path = []
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        real_path.append(Point(x, y))
    l_pos_filter = [Point(measures_pos[0, :][0], measures_pos[0, :][1])]
    # vite = get_velocity(measures_pos, dt)
    # measures_pos = np.concatenate((measures_pos, vite), axis=1)
    measures_pos = np.asmatrix(measures_pos)
    measures_us = np.genfromtxt("mesures_25_sigma10.txt", delimiter=",")
    print measures_pos[0, :].T
    filtering = IvanFilter(measures_pos[0, :].T, dt=dt)
    # var = 10
    # for i in range(1, measures_us.shape[0]):
    #     m1 = measures_us[i, 0]
    #     m2 = measures_us[i, 1]
    #     m3 = measures_us[i, 2]
    #     filtering.update(np.asmatrix([m1, m2, m3]).T)
    #     pos = filtering.get_state_position()
    #     print "pos", pos
    #     l_pos_filter.append(pos)
    var = 0.01
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        print 'x', x, 'y', y
        conv = Converter()
        m1, m2, m3 = conv.get_measures_from_state(x, y)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        print m1, m2, m3
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        l_pos_filter.append(filtering.get_state_position())

    print get_mer(real_path, l_pos_filter)


def script_random_trajectory_unscented_ivan():
    dt = 0.025
    print "script_unscented_trajectory"
    l_points = [[-1000., 200.], [-1000., 800.], [-400., 1200.], [500., 500.], [1100., 180.]]
    real_path = generateur_chemin.generate_path(l_points=l_points, velocity_translation=50,
                                                            velocity_rotation=0.5, dt=dt)
    measures_pos = np.array(real_path)
    real_path_point = []
    for couple in real_path:
        x, y = couple
        pos = Point(x, y)
        real_path_point.append(pos)
    l_pos_filter = [real_path_point[0]]
    measures_pos = np.asmatrix(measures_pos)
    print measures_pos[0, :].T
    filtering = IvanFilter(measures_pos[0, :].T, dt=dt)
    var = 0.01
    for i in range(1, measures_pos.shape[0]):
        x = measures_pos[i, 0]
        y = measures_pos[i, 1]
        print 'x', x, 'y', y
        conv = Converter()
        m1, m2, m3 = conv.get_measures_from_state(x, y)
        m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
        print m1, m2, m3
        filtering.update(np.asmatrix([m1, m2, m3]).T)
        l_pos_filter.append(filtering.get_state_position())

    print get_mer(real_path, l_pos_filter)


if __name__ == "__main__":
    # script_unscented_with_real_measures()
    # print """
    #
    #
    #
    # """
    script_unscented_with_fixed_trajectory_only()
    # script_unscented_trajectory()
    # script_unscented_trajectory_with_file()
    # print """
    #
    #
    #
    #
    # """
    # script_classic_trajectory()
    # script_classic_trajectory_with_real_measures()
    # script_unscented_with_fake_ultrasound_measures()
    # script_unscented_ivan()
    # script_random_unscented_ivan()
    # script_random_trajectory_unscented_ivan()
    # lire_fic_opti()
