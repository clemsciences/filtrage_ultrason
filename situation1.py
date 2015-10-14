from bruitage_conversion import *


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
    return [(conv[0]+ np.random.randn()*var, conv[1] + np.random.randn()*var) for conv in res]
    


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
        u = numpy.zeros(self.x.shape[0])[:, numpy.newaxis]
    self.x = numpy.dot(self.F, self.x) + u
    self.P = numpy.dot(numpy.dot(self.F, self.P), self.F.T) + self.Q
  
  def mesure(self, Z):
    y = Z - numpy.dot(self.H, self.x)
    S = numpy.dot(numpy.dot(self.H, self.P), self.H.T) + self.R    
    K = numpy.dot(numpy.dot(self.P, self.H.T), numpy.linalg.inv(S))
    self.x = self.x + numpy.dot(K, y)
    self.P = numpy.dot((numpy.identity(self.x.shape[0]) - numpy.dot(K, self.H)), self.P)
    
  def filtrer(self, Z, u = None):
    prediction(u)
    measurement(Z)


class ExtendedKalman:
    

class FiltrageLaser:
    
    def __init__(self, config):
        self.config = config
        self.dt = 0.2
        x = numpy.array([1400,100,0.,0.])[:, numpy.newaxis] # vecteur d'état au départ
        P = numpy.matrix([[30.,0.,0.,0.],[0.,30.,0.,0.],[0.,0.,10.,0.],[0.,0.,0.,10.]]) # incertitude initiale
        F = numpy.matrix([[1.,0.,self.dt,0.],[0.,1.,0.,self.dt],[0.,0.,1.,0.],[0.,0.,0.,1.]]) # matrice de transition
        H = numpy.matrix([[1.,0.,0.,0.],[0.,1.,0.,0.]])# matrice d'observation
        R = numpy.matrix([[900,0.],[0.,900]]) # incertitude sur la mesure
        #Q = numpy.matrix([[self.dt**3/3., self.dt**2/2., 0, 0],[self.dt**2/2.,self.dt, 0, 0],[0,0,self.dt**3/3.,self.dt**2/2],[0,0,self.dt**2/2,self.dt]])
        #Q *= 20;
        Q = numpy.matrix([[self.dt**3/3., 0, self.dt**2/2., 0],[0, self.dt**3/3., 0, self.dt**2/2],[self.dt**2/2., 0, 4*self.dt, 0],[0, self.dt**2/2, 0, 4*self.dt]])
        #Q = numpy.matrix([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 4, 0],[0, 0, 0, 4]])
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
            #self.filtre_kalman.measurement(numpy.array([x,y])[:, numpy.newaxis])
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
