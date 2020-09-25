self.scenario   =   {}

self.scenario['T_soil'] = 12 + 273.15

self.scenario[1]                        =   {}
self.scenario[1]['pi']                  =   1
self.scenario[1]['mu']                  =   1
self.scenario[1]['T_sec_in_degree']     =   65
self.scenario[1]['T_sec_in']            =   self.scenario[1]['T_sec_in_degree'] + 273.15 
self.scenario[1]['dotV_sec_in']         =   -10 # l/min
self.scenario[1]['kappa']               =   0
self.scenario[1]['u']                   =   0.8

self.scenario[2]                        =   {}
self.scenario[2]['pi']                  =   1
self.scenario[2]['mu']                  =   -1
self.scenario[2]['T_sec_in_degree']     =   45
self.scenario[2]['T_sec_in']            =   self.scenario[2]['T_sec_in_degree'] + 273.15 
self.scenario[2]['dotV_sec_in']         =   10 # l / min
self.scenario[2]['kappa']               =   0.7
self.scenario[2]['u']                   =   0

self.scenario[3]                        =   {}
self.scenario[3]['pi']                  =   1
self.scenario[3]['mu']                  =   -1
self.scenario[3]['T_sec_in_degree']     =   45
self.scenario[3]['T_sec_in']            =   self.scenario[3]['T_sec_in_degree'] + 273.15 
self.scenario[3]['dotV_sec_in']         =   10 # l / min
self.scenario[3]['kappa']               =   0.5
self.scenario[3]['u']                   =   0

# self.scenario[4]                        =   {}
# self.scenario[4]['pi']                  =   1
# self.scenario[4]['mu']                  =   -1
# self.scenario[4]['T_sec_in_degree']     =   45
# self.scenario[4]['T_sec_in']            =   self.scenario[4]['T_sec_in_degree'] + 273.15 
# self.scenario[4]['dotV_sec_in']         =   10
# self.scenario[4]['kappa']               =   0.7
# self.scenario[4]['u']                   =   0

# self.scenario[5]                        =   {}
# self.scenario[5]['pi']                  =   1
# self.scenario[5]['mu']                  =   -1
# self.scenario[5]['T_sec_in_degree']     =   45
# self.scenario[5]['T_sec_in']            =   self.scenario[5]['T_sec_in_degree'] + 273.15 
# self.scenario[5]['dotV_sec_in']         =   10
# self.scenario[5]['kappa']               =   0.35
# self.scenario[5]['u']                   =   0



