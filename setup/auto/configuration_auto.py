# #keyword == '2':
# self.M = 2
# self.topology = [(2,1)]
# self.coordinates = {'1h': [10,60], '1c':[10,10], '2h': [60,60], '2c': 
                        # [60,10]}
# #keyword == '3_circ':
# self.M = 3
# self.topology = [(1,2), (2,3), (3,1)]
# self.coordinates = {'2h': [35,80], '2c':[35,10], '3h': [60,105], '3c': 
                        # [60,35], '1h': [10,105],'1c':[10,35]}
# keyword == '3_lin':
self.M = 3
self.topology = [(1,2), (2,3)]
self.coordinates = {'1h': [10,60], '1c':[10,10], '2h': [60,60], '2c': 
                        [60,10], '3h': [110,60], '3c':[110,10]}
# # keyword == '4_rhomb':
# self.M = 4
# self.topology = [(1,2), (2,3), (3,4), (4,1), (4,2)]
# self.coordinates = {'1h': [10,105], '1c':[10,35], '2h': [70,120], '2c': 
                        # [70,50], '3h': [100,95], '3c':[100,25],'4h': 
                            # [50,80], '4c':[50,10]}
# keyword == '5_CoSES':
# self.M = 5
# self.topology = [(1,2), (2,5), (5,3), (3,4), (4,1), (4,5)]
# self.coordinates = {'1h': [10,170], '2h': [80,210], '3h': [100,150], '4h': 
                        # [30,110], '5h':[120,190], '1c': [10,70], '2c': 
                            # [80,110], '3c': [100,50], '4c':[30,10], '5c': 
                                # [120,90]}
# # else:
# self.M = 0
# self.topology = [] 
# self.coordinates = {}