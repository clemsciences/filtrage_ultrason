from bruitage_conversion import *
from math import sqrt
import scipy

l_test = [(0.3, 0.5), (-1.2, 1.5), (1.4, 1.7)]

def simuler1(x_reel, y_reel):
    """
    Cette simulation prend la position réelle, bruite les mesures et les convertit pour obtenir les positions possible
    """
    m1, m2, m3 = obtenir_mesures(x_reel, y_reel, var)
    m1, m2, m3 = m1 + np.random.randn()*var, m2 + np.random.randn()*var, m3 + np.random.randn()*var
    res = []
    res.append(equation1_1_2(m1, m2))
    res.append(equation2_1_2(m1, m2))
    res.append(equation1_2_3(m2, m3))
    res.append(equation2_2_3(m2, m3))
    res.append(equation1_3_1(m3, m1))
    res.append(equation2_3_1(m3, m1))
    return res
def simuler2(x_reel, y_reel, var):
    """
    Cette simulation prend la position réelle, convertit les mesures, et bruite les positions
    """
    m1, m2, m3 = obtenir_mesures(x_reel, y_reel)
    res = []
    res.append(equation1_1_2(m1, m2))
    res.append(equation2_1_2(m1, m2))
    res.append(equation1_2_3(m2, m3))
    res.append(equation2_2_3(m2, m3))
    res.append(equation1_3_1(m3, m1))
    res.append(equation2_3_1(m3, m1))
    return [(conv[0]+ dot(np.random.randn()),var, conv[1] + np.dot(np.random.randn(),var)) for conv in res]
    


def g(x):
	"""
	Pour notre modélisation, on choisit comme vecteur x = [abscisse de la position, 
	ordonnée de la position, abscisse de la vitesse, ordonnées de la vitesse]
	"""
	F = np.matrix([[1.,0.,self.dt,0.],[0.,1.,0.,self.dt],[0.,0.,1.,0.],[0.,0.,0.,1.]])
	return dot(F,x)

def h(x):
	return array(obtenir_mesures(x[0],x[1]))[:,np.newaxis]
class UnscentedKalman:
    """
       Cette classe implémente l'algorithme du Filitrage de Kalman Unscented
    Source : Machine Learning a probabilistic perspective p. 651-652
    """
    def __init__(self, g, h, mu0, SIGMA0, d, alpha, beta, kappa, Q, R):
        """
        g est la fonction qui associe la position en t à la position en t+1
        h est la fonction qui associe la position réelle aux mesures correspondantes
        mu0 est la position initiale du robot
        SIGMA est la matrice de covariance initiale
        
        """
        self.d = d # c'est la dimension de l'espace cachée 
        self.alpha = float(alpha) #paramètre bizarre
        self.beta = float(beta)  #paramètre bizarre
        self.kappa = kappa  #paramètre bizarre
        self.mu = mu0 #position initiale
        self.SIGMA = SIGMA0 #variance initiale
        self.lam = alpha**2*(d+kappa) - d # un paramètre un peu étrange aussi
        self.gamma = math.sqrt(d+self.lam) #à nouveau bizarre ! 
        self.Q = Q # matrice de covariance de la gaussienne de l'équation d'évolution
        self.R = R #matrice de covariance de la gaussienne de l'équation d'observation

		#calculs de coefficients
        self.wm0 = self.lam /(self.d*self.lam)
        self.wm = 1/(2.*(self.d+self.lam)
		self.wc0 = self.wm0+(1-self.alpha**2+self.beta)

		
    def premier_pas(self):
		"""
		Estimation de mu_barre et de sigma_barre à partir de l'itération précédente
		"""
		#first Unscented transform
        racine_sigma = scipy.linalg.sqrtm(self.SIGMA)
        self.points_sigma = array([self.mu]+[self.mu - gamma*racine_sigma[:,i] for i in range(0,self.d)]+[self.mu + gamma * racine_sigma[:,i] for i in range(0,d)])
        self.z_etoile_barre = array([g(self.points_sigma[i]) for i in range(0,2*self.d)])
		
        self.mu_barre = array([wm0*self.points_sigma[:,0]]+[wm*self.points_sigma[:,i] for i in range(1,2*self.d)])
       	#calcul de SIGMA_barre
        self.SIGMA_barre = self.wc0*dot( (self.z_etoile_barre[0] - self.mu_barre),(self.z_etoile_barre[0] - self.mu_barre).T
        for i in range(1,2*self.d):
            self.SIGMA_barre += self.wm*dot((self.z_etoile_barre[i] - self.mu_barre),(self.z_etoile_barre[i] - self.mu_barre).T
        self.SIGMA_barre += self.Q

    def second_pas(self, y):
        """
        y est la mesure à l'instant t
		Estimation de mu et de sigma présent
        """
		#second unscented transform
        racine_sigma = scipy.linalg.sqrtm(self.SIGMA_barre)
        self.points_sigma_barre = array([self.mu_barre]+[self.mu - gamma*racine_sigma[:,i] for i in range(0,self.d)]+[self.mu_barre + gamma * racine_sigma[:,i] for i in range(0,d)])
		#
        self.y_etoile_barre = array([h(self.points_sigma_barre[i]) for i in range(0,2*self.d)])
        
        self.y_chapeau = array([self.wm0*self.points_sigma_barre[:,0]]+[self.wm*self.points_sigma_barre[:,i] for i in range(1,2*self.d)])
		#Calcul de S
        self.S = self.wc0*dot((self.y_etoile_barre[0] - self.y_chapeau),(self.y_etoile_barre[0] - self.y_chapeau).T)
        for i in range(1,2*self.d):
            self.S += wm*dot((self.y_etoile_barre[i] - self.y_chapeau),(self.y_etoile_barre[i] - self.y_chapeau).T)
        self.S += self.R
		#calcul de SIGMA_z_y
        self.SIGMA_z_y_barre = wc0*dot((self.z_etoile_barre[0] - self.mu_barre),(self.y_etoile_barre[0] - self.y_chapeau).T)
        for i in range(1, 2*self.D):
            self.SIGMA_z_y_barre += wm*dot((self.z_etoile_barre[i] - self.mu_barre),(self.y_etoile_barre[i] - self.y_chapeau).T)
        self.K = dot(self.SIGMA_z_y_barre, linalg.inv(self.S))
		#Les valeurs qui nous interessent
        self.mu = self.mu_barre + dot(self.K, (y - self.y_chapeau))
        self.SIGMA = self.SIGMA_barre - dot(dot(self.K,self.S),T)
	def filtrer(self, y):
		"""
		Pöurquoi séparer les deux étapes? Je l'ignore
		"""
		self.premier_pas()
		self.second_pas()
#classe à arranger pour le Kalman Unscented
class FiltrageLaserUnscented:
    
    def __init__(self, config, x0):
        self.config = config
        self.dt = 0.2
        x = np.array(x).T
        #x = np.array([1400,100,0.,0.])[:, np.newaxis] # vecteur d'état au départ
        P = np.matrix([[30.,0.,0.,0.],[0.,30.,0.,0.],[0.,0.,10.,0.],[0.,0.,0.,10.]]) # incertitude initiale
        F = np.matrix([[1.,0.,self.dt,0.],[0.,1.,0.,self.dt],[0.,0.,1.,0.],[0.,0.,0.,1.]]) # matrice de transition
        H = np.matrix([[1.,0.,0.,0.],[0.,1.,0.,0.]])# matrice d'observation
        R = np.matrix([[900,0.],[0.,900]]) # incertitude sur la mesure
        #Q = np.matrix([[self.dt**3/3., self.dt**2/2., 0, 0],[self.dt**2/2.,self.dt, 0, 0],[0,0,self.dt**3/3.,self.dt**2/2],[0,0,self.dt**2/2,self.dt]])
        #Q *= 20;
        Q = np.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2],[self.dt**2/2., 0, 4*self.dt, 0],[0, self.dt**2/2, 0, 4*self.dt]])
        #Q = np.matrix([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 4, 0],[0, 0, 0, 4]])
        Q *= 30
        self.filtre_kalman = Kalman(x, P, F, H, R, Q)
        self.historique = collections.deque(maxlen=3)
        self.valeurs_rejetees = 0
        self.acceleration = None
        
    def etat_robot_adverse(self):
        return self.filtre_kalman.x
        
    def update_dt(self, new_dt):
        self.dt = new_dt
        self.filtre_kalman.F[0,2] = new_dt
        self.filtre_kalman.F[1,3] = new_dt
    
    def position(self):
        #state = self.filtre_kalman.x
        #return Point(int(state[0]), int(state[1]))
        return self.last_point
    
    def vitesse(self):
        state = self.filtre_kalman.x
        return Vitesse(int(state[2]), int(state[3]))
                
    def update(self, x, y):
        if self._filtrage_acceleration(Point(x, y)):
            self.last_point = Point(x, y)
            #self.filtre_kalman.prediction()
            #self.filtre_kalman.measurement(np.array([x,y])[:, np.newaxis])
            #self.historique.append(self.position())
        #else:
            #self.last_point = None
            #self.filtre_kalman.prediction()
        
    def _filtrage_acceleration(self, pointm0):
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
  
  def __init__(self,x,P,F,H,R,Q):
    self.x = x
    self.P = P
    self.F = F
    self.H = H
    self.R = R
    self.Q = Q
  
  def prediction(self, u = None):
    if u == None:
        u = np.zeros(self.x.shape[0])[:, np.newaxis]
    self.x = np.dot(self.F, self.x) + u
    self.P = np.dot(np.dot(self.F, self.P), self.F.T) + self.Q
  
  def mesure(self, Z):
    y = Z - np.dot(self.H, self.x)
    S = np.dot(np.dot(self.H, self.P), self.H.T) + self.R    
    K = np.dot(np.dot(self.P, self.H.T), np.linalg.inv(S))
    self.x = self.x + np.dot(K, y)
    self.P = np.dot((np.identity(self.x.shape[0]) - np.dot(K, self.H)), self.P)
    
  def filtrer(self, Z, u = None):
    prediction(u)
    measurement(Z)


class FiltrageLaser:
    
    def __init__(self, config, x0):
        self.config = config
        self.dt = 0.2
        x = np.array(x).T
        #x = np.array([1400,100,0.,0.])[:, np.newaxis] # vecteur d'état au départ
        P = np.matrix([[30.,0.,0.,0.],[0.,30.,0.,0.],[0.,0.,10.,0.],[0.,0.,0.,10.]]) # incertitude initiale
        F = np.matrix([[1.,0.,self.dt,0.],[0.,1.,0.,self.dt],[0.,0.,1.,0.],[0.,0.,0.,1.]]) # matrice de transition
        H = np.matrix([[1.,0.,0.,0.],[0.,1.,0.,0.]])# matrice d'observation
        R = np.matrix([[900,0.],[0.,900]]) # incertitude sur la mesure
        #Q = np.matrix([[self.dt**3/3., self.dt**2/2., 0, 0],[self.dt**2/2.,self.dt, 0, 0],[0,0,self.dt**3/3.,self.dt**2/2],[0,0,self.dt**2/2,self.dt]])
        #Q *= 20;
        Q = np.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2],[self.dt**2/2., 0, 4*self.dt, 0],[0, self.dt**2/2, 0, 4*self.dt]])
        #Q = np.matrix([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 4, 0],[0, 0, 0, 4]])
        Q *= 30
        self.filtre_kalman = Kalman(x, P, F, H, R, Q)
        self.historique = collections.deque(maxlen=3)
        self.valeurs_rejetees = 0
        self.acceleration = None
        
    def etat_robot_adverse(self):
        return self.filtre_kalman.x
        
    def update_dt(self, new_dt):
        self.dt = new_dt
        self.filtre_kalman.F[0,2] = new_dt
        self.filtre_kalman.F[1,3] = new_dt
    
    def position(self):
        #state = self.filtre_kalman.x
        #return Point(int(state[0]), int(state[1]))
        return self.last_point
    
    def vitesse(self):
        state = self.filtre_kalman.x
        return Vitesse(int(state[2]), int(state[3]))
                
    def update(self, x, y):
        if self._filtrage_acceleration(Point(x, y)):
            self.last_point = Point(x, y)
            #self.filtre_kalman.prediction()
            #self.filtre_kalman.measurement(np.array([x,y])[:, np.newaxis])
            #self.historique.append(self.position())
        #else:
            #self.last_point = None
            #self.filtre_kalman.prediction()
        
    def _filtrage_acceleration(self, pointm0):
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


class Point:
    """
    Classe permettant de définir un point mathématique dans R^2 et les opérations usuelles sur les points dans R^2
    
    :param x: abscisse
    :type x: float
    
    :param y: ordonnée
    :type y: float
    
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def __repr__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'
        
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)
        
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)
        
    def __rmul__(self, other):
        return Point(self.x*other, self.y*other)
    
    def __str__(self) :
        return "("+str(self.x)+"," + str(self.y) + ")"
        
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
        
    def to_list(self):
        return [self.x, self.y]
        
    def unitaire(self):
        n = self.norme()
        return Point(self.x / n, self.y / n)
        
    def norme(self):
        return self.distance(Point(0, 0))
        
    def distance(self, point):
        v = self - point
        d = v.x ** 2 + v.y ** 2
        return math.sqrt(d)
    
    def copy(self):
        return Point(self.x, self.y)

class Vitesse:
    """
    Classe permettant de définir une vitesse sous forme de vecteur mathématique dans R^2 et les opérations usuelles sur les vecteurs dans R^2
    
    :param vx: vitesse sur x
    :type vx: float
    
    :param vy: vitesse sur y
    :type vy: float
    
    """
    def __init__(self, vx, vy):
        self.vx = vx
        self.vy = vy
    
    def __repr__(self):
        return '(' + str(self.vx) + ',' + str(self.vy) + ')'
        
    def __add__(self,other):
        return Vitesse(self.vx + other.vx, self.vy + other.vy)
        
    def __sub__(self,other):
        return Vitesse(self.vx - other.vx, self.vy - other.vy)
    
    def __mul__(self,other):
        return Vitesse(self.vx*other, self.vy*other)
    
    def __str__(self) :
        return "("+str(self.vx)+"," + str(self.vy) + ")"
        
    def norme(self):
        return math.sqrt(self.vx ** 2 + self.vy ** 2)
        
    def to_list(self):
        return [self.vx, self.vy]
if __name__ == "__main__":
    for var in [0.001, 0.01, 0.1]:
        for i in l_test:
            x_reel = i[0]
            y_reel = i[1]
            res1 = simuler1(x_reel, y_reel)
            res2 = simuler2(x_reel, y_reel)
