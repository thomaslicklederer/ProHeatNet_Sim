# define standard components
mypump     =      cm.pump()
myvalve    =      cm.controlvalve()
myhx       =      cm.heatexchanger()
mypipe1    =      cm.pipe()
mypipe2    =      cm.pipe()


# define fluid that is used
myfluid         =   cm.fluid()
self.myfluid    =   myfluid

# define components for each edge automatically
self.components =      {}
for edge in e_vec:
    self.components[edge] = {}
    self.components[edge]['pipe'] = mypipe1
    # filter for connections that connect warm and cold subnetwork and add 
    # devices
    if (('h' in edge[0] and 'c' in edge[1]) or ('h' in edge[1] and 'c' in 
                                                edge[0])):
        self.components[edge]['heatexchanger'] = myhx
        self.components[edge]['pump'] = mypump
        self.components[edge]['controlvalve'] = myvalve
        self.components[edge]['pipe'] = mypipe1
