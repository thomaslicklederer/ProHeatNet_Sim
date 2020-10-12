## inputs per prosumer:    pi, mu, T_sec_in, dotV_sec_in, kappa, u
## inputs per edge:        V_min, V_max
## general inputs: T_soil

# for random generator
dicenbr = 7

# T_soil
T_soil = 12+273.15

# generation of (random) values for scenario
random.seed(dicenbr)
# pi
pi_vec  =   [1 for i in range(self.M)]

#mu
mu_vec  =   [99 for i in range(self.M)]
b   =   int(math.ceil(self.M/2.0)) # more consumers than producers, b=number consumers
a   =   int(self.M-b) # a = number of producers
mu_vec[0:a]   =   [+1 for i in mu_vec[0:a]]  # first half of prosumers is producer
mu_vec[a:]    =   [-1 for i in mu_vec[a:]]   # second half of prosumers is consumer
random.shuffle(mu_vec)

#T_sec_in
T_min_w = 60+273.15
T_max_w = 80+273.15
T_min_c = 30+273.15
T_max_c = 50+273.15
T_sec_in_vec = list(np.zeros(np.shape(mu_vec)))
for entry in range(len(T_sec_in_vec)):
    if mu_vec[entry] == -1:
        T_sec_in_vec[entry]=45+273.15 # ((T_max_c-T_min_c)*random.random()+T_min_c)
    elif mu_vec[entry] == +1:
        T_sec_in_vec[entry]=65+273.15 # ((T_max_w-T_min_w)*random.random()+T_min_w)
T_sec_in_degree_vec = list(np.array(T_sec_in_vec) - 273.15)
 
#dotV_sec_in
dotV_sec_in_vec = list(np.zeros(np.shape(mu_vec)))
dotV_norm   =   math.pi*((0.022)**2)*(1/4)*1*60*1000 # liters per minute
for entry in range(len(dotV_sec_in_vec)):
    if mu_vec[entry] == -1:
        dotV_sec_in_vec[entry]=+14 # +(0.8+0.25*random.uniform(-1.0, 1.0))*dotV_norm
    elif mu_vec[entry] == +1:
        dotV_sec_in_vec[entry]=-14 #-(0.8+0.25*random.uniform(-1.0, 1.0))*dotV_norm
 
#kappa
kappa_vec = list(np.zeros(np.shape(mu_vec)))
kappa_cons_norm =   0.9
for entry in range(len(kappa_vec)):
    if mu_vec[entry] == -1:
        kappa_vec[entry]=0.9 #random.uniform(-0.2, 0.1)+kappa_cons_norm
    elif mu_vec[entry] == +1:
        kappa_vec[entry]= 0

#u
u_vec = list(np.zeros(np.shape(mu_vec)))
u_prod_norm =   0.8
for entry in range(len(u_vec)):
    if mu_vec[entry] == -1:
        u_vec[entry]=0
    elif mu_vec[entry] == +1:
        u_vec[entry]=0.8  # random.uniform(-0.2, 0.2)+u_prod_norm

##
self.scenario   =   {}
for prosumer in range(self.M):
    self.scenario[prosumer+1]   =   {}
    self.scenario[prosumer+1]['pi']   =   pi_vec[prosumer]
    self.scenario[prosumer+1]['mu']   =   mu_vec[prosumer]
    self.scenario[prosumer+1]['T_sec_in']   =   T_sec_in_vec[prosumer]
    self.scenario[prosumer+1]['T_sec_in_degree']   =   T_sec_in_degree_vec[prosumer]
    self.scenario[prosumer+1]['dotV_sec_in']   =   dotV_sec_in_vec[prosumer]
    self.scenario[prosumer+1]['kappa']   =   kappa_vec[prosumer]
    self.scenario[prosumer+1]['u']   =   u_vec[prosumer]
self.scenario['T_soil'] = T_soil




        