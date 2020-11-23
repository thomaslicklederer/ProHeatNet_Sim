self.scenario   =   {}

self.scenario['T_soil'] = 20 + 273.15

self.scenario[1]                        =   {}
self.scenario[1]['pi']                  =   1
self.scenario[1]['mu']                  =   -1
self.scenario[1]['T_sec_in_degree']     =   40  # [°C]
self.scenario[1]['T_sec_in']            =   self.scenario[1]['T_sec_in_degree'] + 273.15 # [K]
self.scenario[1]['dotV_sec_in']         =   8 # [l / min]
self.scenario[1]['kappa']               =   1
self.scenario[1]['u']                   =   0

self.scenario[2]                        =   {}
self.scenario[2]['pi']                  =   1
self.scenario[2]['mu']                  =   1
self.scenario[2]['T_sec_in_degree']     =   65  # [°C]
self.scenario[2]['T_sec_in']            =   self.scenario[2]['T_sec_in_degree'] + 273.15 # [K]
self.scenario[2]['dotV_sec_in']         =   -5 # [l/min]
self.scenario[2]['kappa']               =   0
self.scenario[2]['u']                   =   0.7

self.scenario[3]                        =   {}
self.scenario[3]['pi']                  =   1
self.scenario[3]['mu']                  =   1
self.scenario[3]['T_sec_in_degree']     =   65  # [°C]
self.scenario[3]['T_sec_in']            =   self.scenario[3]['T_sec_in_degree'] + 273.15 # [K]
self.scenario[3]['dotV_sec_in']         =   -5 # [l / min]
self.scenario[3]['kappa']               =   0
self.scenario[3]['u']                   =   0.7




