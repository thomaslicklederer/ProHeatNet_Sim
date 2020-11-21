# example configuration

## number of prosumers in network
self.M = 3

## topology between prosumers:
##  (1,2) means a connection between prosumer 1 and 2, with positive direction from 1 to 2
self.topology = [(1,2), (2,3)]

## coordinates of nodes - just for visualization purposes
self.coordinates = {'2h': [35,80], '2c':[35,10], '3h': [60,105], '3c': 
                        [60,35], '1h': [10,105],'1c':[10,35]}