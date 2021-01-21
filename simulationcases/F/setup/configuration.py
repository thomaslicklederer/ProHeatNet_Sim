# example configuration

## number of prosumers in network
self.M = 3

## topology between prosumers:
##  (1,2) means a connection between prosumer 1 and 2, with positive direction from 1 to 2
self.topology = [(1,2), (2,3), (1,3)]

## coordinates of nodes - just for visualization purposes
self.coordinates = {'1h': [10,60], '1c':[10,10], '2h': [60,60], '2c': 
                        [60,10], '3h': [110,60], '3c':[110,10]}