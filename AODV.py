from random import randrange, random
import numpy
import copy
__version__ = "0.1"
__release__ = "2011-09-12-r2" # last update
__author__= "lujiasen"


class node:
    '''An object of class which contains all the information of a node. 
    '''
    def __init__(self, Length, Width, Id):
        self.Id = Id
        self.Length = Length
        self.Width = Width
        self.Location = [randrange(Length), randrange(Width)]
        self.NLocation = self.Location
        self.Speed = 0
        self.MaxID = 0
        self.Buffer = 30

        #----------------parameter--------------------- 
        self.PktRate = 3
        self.MaxSpeed = 5
        self.PauseTime = 20
        self.PauseTimer = 0
        self.ExpireTime = 20
        #----------------parameter--------------------
        #-------------------------
        self.RREQ_State = []
        self.RREP_State = []
        self.Packet_State = []
        #------------------------
        self.RREQ_RoutingTable = []
        self.RREP_RoutingTable = []
        self.RREQ_Table_Time = {}
        self.RREP_Table_Time = {}
        # the routingtable contains the destination, nextcount and the trail.
        self.RREP_recieve = []
        self.PktRemain = {}
        # the PktRemain contains the destination, RREQ.Id and the remain to forward packet.
        self.Recieve_Pkt = []
        
        self.Timeout = {}

        # the Timeout recorder when sending the RREQ, if the RREQ wasn't recieve at a certain time, the node would 
        # regenerate the RREQ and reset the Time 


    def gen_Speed(self):
        self.Speed = randrange(0,self.MaxSpeed)

    def gen_NLocation(self):
        self.NLocation = [randrange(self.Length), randrange(self.Width)]
    
    def add_PauseTimer(self):
        self.PauseTimer += 1
	
    def clear_PauseTimer(self):
        self.PauseTimer = 0
    
    def update_RREQ_State(self, Pkt):
        if self.Id != Pkt.Dst:
            self.RREQ_State.append([Pkt.Id, Pkt.Src])

    def update_Packet_State(self, i, j, k):
        '''
        i:Packet id
        j:Packet Num
        k:Packet Src
        '''
        self.Packet_State.append([i,j,k])

    def update_Location(self):
        if self.Location == self.NLocation: 
        # if the node reach the next location
            if self.PauseTime > self.PauseTimer:
            # the node still in pause time
                self.add_PauseTimer()
                # add it's pause time
            else:
                # the node out of the pause time, the node choose Next location.
                self.gen_NLocation()
                self.gen_Speed()
                self.clear_PauseTimer()
        else:
        # if the node is not reach the next location.
        # end of class node.
            X = self.NLocation[0] - self.Location[0]
            Y = self.NLocation[1] - self.Location[1]
            Distance = numpy.sqrt(X ** 2 + Y ** 2)

            if Distance > self.Speed:
                self.Location = [(self.Location[0] + X * self.Speed / Distance),
                                (self.Location[1] + Y * self.Speed / Distance)]
            else:
                self.Location = self.NLocation

    def recieve_RREP(self, Pkt):
        self.RREP_recieve.append([Pkt.Id, Pkt.Dst, Pkt.Path])

    def update_RREQ_RoutingTable(self, i, j, k):
        '''
        update the Routing Table of a node, the Routing table has the Destination, Last Hop
        
        Parameters
        ----------
        i: integer
            source 
        j: integer
            RREQ id
        k: integer
            Reverse node number
                
        '''
        Flag = 0
        for m in xrange(0, len(self.RREQ_RoutingTable)):
            if i == self.RREQ_RoutingTable[m][0]:
                if j > self.RREQ_RoutingTable[m][1]:
                    self.RREQ_RoutingTable[m][2] = k
                Flag = 1
        # update the routingTable when new RREQ updated.
        if Flag == 0:
            self.RREQ_RoutingTable.append([i,j,k])

    def update_RREP_RoutingTable(self, i,j,k):
        Flag = 0
        for m in xrange(0, len(self.RREP_RoutingTable)):
            if i == self.RREP_RoutingTable[m][0]:
                if j >= self.RREP_RoutingTable[m][1]:
                    self.RREP_RoutingTable[m][2] = k
                Flag = 1
        if Flag == 0:
            self.RREP_RoutingTable.append([i,j,k])

    def set_RoutingTable_Time(self):
        for i in xrange(0, len(self.RREQ_RoutingTable)):
            if tuple(self.RREQ_RoutingTable[i]) not in self.RREQ_Table_Time:
                self.RREQ_Table_Time[tuple(self.RREQ_RoutingTable[i])] = self.ExpireTime
            
        for i in xrange(0, len(self.RREP_RoutingTable)):
            if tuple(self.RREP_RoutingTable[i]) not in self.RREP_Table_Time:
                self.RREP_Table_Time[tuple(self.RREP_RoutingTable[i])] = self.ExpireTime

    def update_RoutingTable_Time(self):
        

        for i in self.RREQ_Table_Time.iteritems():
            self.RREQ_Table_Time[i[0]] -= 1
            
        j = 0
        Temp_RREQ = self.RREQ_Table_Time.items()

        while True:
            if j < len(Temp_RREQ):
                if Temp_RREQ[j][1] <= 0:
                    if list(Temp_RREQ[j][0]) in self.RREQ_RoutingTable:
                        self.RREQ_RoutingTable.remove(list(Temp_RREQ[j][0]))

                    del self.RREQ_Table_Time[tuple(Temp_RREQ[j][0])]
                    
                    Temp_RREQ.pop(j)
                    j -= 1
                j += 1
            else:
                break            
               
                
        for i in self.RREP_Table_Time.iteritems():
            self.RREP_Table_Time[i[0]] -= 1


        j = 0
        Temp_RREP = self.RREP_Table_Time.items()

        while True:
            if j < len(Temp_RREP):
                if Temp_RREP[j][1] <= 0:
                    if list(Temp_RREP[j][0]) in self.RREP_RoutingTable:
                        self.RREP_RoutingTable.remove(list(Temp_RREP[j][0]))

                    del self.RREP_Table_Time[tuple(Temp_RREP[j][0])]
                    
                    Temp_RREP.pop(j)
                    j -= 1
                j += 1
            else:
                break            

   #-------------------------------------------------------------------     
    def set_PktRemain(self, Pkt):
        if (Pkt.Src, Pkt.Id) not in self.PktRemain:
            self.PktRemain[Pkt.Src, Pkt.Id] = 5
   #--------------------------------------------------------------------
    def set_Timeout(self, i, j):
   
        if (i, j) not in self.Timeout:
            self.Timeout[i, j] = 40
    
    def del_Timeout(self, i, j):
        if (i, j) in self.Timeout:
            del self.Timeout[i, j]     
    
    def update_Timeout(self, i, j):
        if (i, j) in self.Timeout:
            self.Timeout[i, j] -= 1


class RREQ:
    '''An object of class RREQ is a route discovery request packet from a node
    to another.[Pk

    An RREQ package has the follow
        for Nbring variables: Id, Src, Dst, TTL and Prb. 
    They can be instantiated by the constructor when intializing an object.

    Parameters
    ----------.RREQ_State
    Id : integer
        package id of an RREQ

    Src : integer
        the source node of an RREQ

    Dst : integer
        the destination node of an RREQ

    TTL : integer
        Time to live of an RREQ

    Prb : float
        gossiping probability of an RREQ

    Methods
    -------
    reduce_TTL()
        reduce the TTL of an object by 1
    
    '''
    
    def __init__(self, Id, Src, Dst, TTL, Prb):
        self.Id = Id
        self.Src = Src
        self.Dst = Dst
        self.TTL = TTL
        self.Prb = Prb
        self.HopCount = 0
       
         
        self.Flag = 0
        
    def reduce_TTL(self):
        '''this member function is unused right now.     
        
        '''
        self.TTL -= 1


    def add_HopCount(self):
        self.HopCount += 1


class RREP:
    def __init__(self, Id, Src, Dst, TTL):
        self.Id = Id
        self.Src = Src
        self.Dst = Dst
        self.Path = []
        self.Flag = 1
        self.TTL = TTL
        self.HopCount = 0
        self.Path = []

    def add_HopCount(self):
        self.HopCount += 1

        # the TTL is just in case for convinient, no other specific usage.tps://wlan.ct10000.com/
class Packet:
    def __init__(self, Id, Src, Dst, Prb, TTL, Num):
        self.Id = Id
        self.Src = Src
        self.Dst = Dst
        self.Prb = Prb
        self.Flag = 2
        self.TTL = TTL
        self.Num = Num
        self.RouteRecovery = 4
        # same with the RREP's TTL
        # packet Num is to identify the packet, the packet with the same Id and Num the packet is the same.

class RRER:
    def __init__(self, Id, Src, Dst, TTL):
        self.Id = Id
        self.Src = Src
        self.Dst = Dst
        self.TTL = TTL
        self.Flag = 3

def gen_nodes(Length = 1000, Width = 1000, Number = 100):
    '''Initialize all the nodes.
    '''
	
    Nodes = [0] * Number
    for i in xrange(0, Number):
        Nodes[i] = node(Length, Width, i)
    
    return Nodes

def genRQ(Number):
    
    return [[] for i in xrange(0, Number)]

def init(Length = 1000, Width = 1000, Number = 100):
    
    Nodes = gen_nodes(Length, Width, Number)
    RQ = genRQ(Number)
    return Nodes, RQ

#--------------------------------basic setup--------------------------
#------------------------------------------------------------------------------

def dist(X, Y):
    '''Calculate Euclidean distance between two 2-D coordinates.

    Parameters
    ----------

    X : list of two integers
        one coordinate 

    Y : list of two integers
        the other coordinate 

    '''

    return numpy.linalg.norm(numpy.array(X) - numpy.array(Y))

def gen_Adj(Nodes, TxRng):
    '''Generate adjcent matrix for *Nodes* in transmission range *TxRng*.

    Parameters
    -----------

    Nodes : list of lists (2-tuples) of integers
        coordinates of nodes

    TxRng : float
        transmission range of nodes

    Num : integer
        number of nodes in the network

    Returns
    --------
    Adj: list of lists of Bools
        adjacent matrix of nodes, Adj[i][j] = True if nodes i and j are neighbors

    '''

    Num = len(Nodes) # number of nodes
    Adj = [[(dist(Nodes[i].Location, Nodes[j].Location) <= TxRng and i != j) for j in xrange(0, Num)] for i in xrange(0, Num)]

    return Adj

def update_Nodes(Nodes):
    '''update the information of Nodes, include the NLocation, Grid ect.
    '''
    for i in xrange(0,len(Nodes)):
        Nodes[i].update_Location()
        Nodes[i].set_RoutingTable_Time()
        Nodes[i].update_RoutingTable_Time()
    return Nodes

def prepare(Nodes, RQ, TxRng, InitTTL):
    '''
    prepare Adj, RQ, RQp for the simulation.
    '''
    Nodes = update_Nodes(Nodes)
    Adj = gen_Adj(Nodes, TxRng)
    Nodes, RQ = check_RQ(Nodes, RQ)
    RQp = prepare_RQp(RQ)
    return Nodes, Adj, RQ, RQp    

def check_RQ(Nodes, RQ):
    '''check all the packet in RQ, if the RREQ packet is out the TTL, if the RREQ packet arrive the destination 
    node, then the destination node become the Source node and the Source node become the destination node and generate the same ID RREP in RQ,
    the RREP was sent along the RREQ path, when the RREP arrive the desination(original source) node, the nodes 
    start to sending the the packet to the destination node.
    
    '''
    for i in xrange(0, len(RQ)):
        if RQ[i] != []:
            j = 0
            while True:
                if j < len(RQ[i]):  
                    Pkt = RQ[i][j]
                    if Pkt.Dst == i or Pkt.TTL <= 0: 
                        RQ[i].pop(j)
                        if Pkt.Dst == i:
                            if Pkt.Flag == 0:
                            # suppose it's RREQ, then it should foward the RREP
                                NPkt = RREP(Pkt.Id, Pkt.Dst, Pkt.Src, 30)
                                Nodes[i].RREP_RoutingTable.append([Pkt.Dst, Pkt.Id, 'Null'])
                                RQ[i].insert(0,NPkt)
                                
                            if Pkt.Flag == 1:
                            # it's RREP, the node set the PktRemain to forward the packet,
                                Nodes[i].set_PktRemain(Pkt)
                                Nodes[i].del_Timeout(Pkt.Id, Pkt.Src)
                            if Pkt.Flag == 2:
                                if [Pkt.Src, Pkt.Id, Pkt.Num] not in Nodes[i].Recieve_Pkt:
                                    Nodes[i].Recieve_Pkt.append([Pkt.Src, Pkt.Id, Pkt.Num])
                              #  print Pkt.Id
                             #   print Pkt.Num
                        j -= 1
                    j += 1
                else:
                    break
    return Nodes, RQ
                                
def prepare_RQp(RQ):
    RQp = copy.deepcopy(RQ)
    for i in xrange(0, len(RQ)):
        if RQ[i] != []:
            RQp[i] = RQp[i][:-1]
    return RQp
    
#-----------------------------------------------------------------------------------------------------------

def check_Timeout(Nodes, RQ, RQNum, InitTTL):
    '''
    if the Node first send the RREQ, it set an timeout table for it.
    if the Node didn't recieve the RREP in that time, it reforward the 
    RREQ, with different RREQ Id.
    '''
    for i in xrange(0, len(RQ)):
        if Nodes[i].Timeout != {}:
            Temp = Nodes[i].Timeout.items() 
            for j in xrange(len(Temp)):
                if Temp[j][1] > 0:
                    Nodes[i].update_Timeout(Temp[j][0][0], Temp[j][0][1])
                else:
                    Pkt = RREQ(Temp[j][0][0] + 1, i, Temp[j][0][1], InitTTL, 1)
                    Nodes[i].RREP_RoutingTable.append([Pkt.Src, Pkt.Id, 'Null'])
                    Nodes[i].update_RREQ_State(Pkt)
                    Nodes[i].set_Timeout(Pkt.Id, Pkt.Dst)
                    RQ[i].insert(0, Pkt)

    return Nodes, RQ
def run(Nodes, RQ, RQNum, MaxIter, TxRng, InitTTL):
    '''
    Simulate the running of the network.
    '''
    Iter = 0
    Nodes, RQ = gen_RREQ(Nodes, RQ, RQNum, InitTTL)
    while Iter < MaxIter:
        Iter += 1
        Nodes, Adj, RQ, RQp = prepare(Nodes, RQ, TxRng, InitTTL)
        Nodes, RQ = gen_Pkt(Nodes, RQ)
        Nodes, RQ = check_Timeout(Nodes, RQ, RQNum, InitTTL)
        Nodes, RQ = broadcast(Nodes, Adj, RQ, RQp)
        
        print 'Iter= %d' %Iter


    return Nodes, RQ

def gen_RREQ(Nodes, RQ, RQNum, InitTTL):
    InitTTL = 30
    for i in xrange(RQNum):
    
        Src = randrange(len(RQ))
        Dst = randrange(len(RQ))
        while Dst == Src:
            Dst = randrange(len(RQ))
            if Dst != Src:
                break
        Nodes[Src].MaxID += 1
        Pkt = RREQ(Nodes[Src].MaxID, Src, Dst, InitTTL, 1)
        Nodes[Src].RREQ_RoutingTable.append([Src, Pkt.Id, 'Null'])
        Nodes[Src].update_RREQ_State(Pkt)
        Nodes[Src].set_Timeout(Pkt.Id, Pkt.Dst)
        RQ[Src].insert(0, Pkt)

    return Nodes, RQ

def gen_Pkt(Nodes, RQ):
    for i in xrange(0, len(RQ)):
        if Nodes[i].PktRemain != {}:
            Temp = Nodes[i].PktRemain.items() 
            for j in xrange(len(Temp)):
                if Temp[j][1] > 0:
                    NPkt = Packet(Temp[j][0][1], i, Temp[j][0][0], 1, 30, 5-Nodes[i].PktRemain[Temp[j][0]])
                    # here the 10 is the total packet to send.
                    RQ[i].insert(0, NPkt)
                    Nodes[i].update_Packet_State(Temp[j][0][1], 5-Nodes[i].PktRemain[Temp[j][0]], i)
                    Nodes[i].PktRemain[Temp[j][0]] -= 1

    return Nodes, RQ

def nbrs(Adj, i):
    '''Determine the neighbor of a node i
    '''
    Nbrs = []

    for j in xrange(0, len(Adj[i])):
        if Adj[i][j]:
            Nbrs.append(j)
    return Nbrs

def broadcast(Nodes, Adj, RQ, RQp):
    '''broadcast the RREQ from the source node to the destination node using the different routing algrithms
    '''
    for i in xrange(0, len(RQ)):
        if True:
            if RQ[i] != []:
                Pkt = RQ[i][-1]

                Nbrs = nbrs(Adj, i)
                #-----------------for different algrothm---------------#

                # with probability algrithm.
                if Pkt.Flag == 0:
                # the packet is RREQ, use the RREQ broadcast algrithm.
                    Nodes, RQp = forward_RREQ(Nodes, Adj, RQp, Pkt, Nbrs, i)

                # if the node hasn't recieve the same RREP.
                if Pkt.Flag == 1:
                # the packet is RREP, broadcast along the RREQ way.
                    Nodes, RQp = forward_RREP(Nodes, Adj, RQp, Pkt, Nbrs, i)
                    
                if Pkt.Flag == 2:
                # the packet is Packet, broadcast along the RREQ way.
                    Nodes, RQp = forward_Pkt(Nodes, Adj, RQp, Pkt, Nbrs, i)
                if Pkt.Flag == 3:
                    Nodes, RQp = forward_RRER(Nodes, Adj, RQp, Pkt, Nbrs, i) 

    return Nodes, RQp

def forward_RREQ(Nodes, Adj, RQp, Pkt, Nbrs, i):            
    #if Pkt.Prb >= random():
    Flag = 0
    for Nbr in Nbrs:
        if Pkt.Dst == Nbr and len(RQp[i])<=Nodes[i].Buffer:
            Flag = 1
            Nodes[Nbr].update_RREQ_State(Pkt)
            NPkt = RREQ(Pkt.Id, Pkt.Src, Pkt.Dst, Pkt.TTL - 1, 1)
            NPkt.add_HopCount()
            Nodes[Nbr].update_RREQ_RoutingTable(NPkt.Src, NPkt.Id, i)
            RQp[Nbr].insert(0,NPkt)
            # if there is the destination node in the neighbor node.

    if Flag == 0:
        for Nbr in Nbrs:
            if [Pkt.Id, Pkt.Src] not in Nodes[Nbr].RREQ_State and len(RQp[i]) <= Nodes[i].Buffer :

            # if the node hasn't recieve the same RREQ.
                Nodes[Nbr].update_RREQ_State(Pkt)
                NPkt = RREQ(Pkt.Id, Pkt.Src, Pkt.Dst, Pkt.TTL - 1, 1)
                NPkt.add_HopCount()
                Nodes[Nbr].update_RREQ_RoutingTable(NPkt.Src, NPkt.Id, i)
                RQp[Nbr].insert(0,NPkt)

    return Nodes, RQp

def forward_RREP(Nodes, Adj, RQp, Pkt, Nbrs, i):
    
    for j in xrange(0, len(Nodes[i].RREQ_RoutingTable)):
        if Pkt.Dst == Nodes[i].RREQ_RoutingTable[j][0]:
        # if there is destination information in Routingtable
            k = Nodes[i].RREQ_RoutingTable[j][2]
            if k in Nbrs:
                if len(RQp[k]) <= Nodes[k].Buffer:
                    NPkt = RREP(Pkt.Id, Pkt.Src, Pkt.Dst, Pkt.TTL - 1)
                    Nodes[k].update_RREP_RoutingTable(NPkt.Src, NPkt.Id, i)
                    RQp[k].insert(0, NPkt)
            else:
                for Nbr in Nbrs:
                    if len(RQp[Nbr]) <= Nodes[Nbr].Buffer:
                        NPkt = RREP(Pkt.Id, Pkt.Src, Pkt.Dst, Pkt.TTL - 1)
                        Nodes[k].update_RREP_RoutingTable(NPkt.Src, NPkt.Id, i)
                        NPkt.add_HopCount()
                        RQp[Nbr].insert(0, NPkt)

    return Nodes, RQp

def forward_Pkt(Nodes, Adj, RQp, Pkt, Nbrs, i):
    '''The packet was forward base on the RREP routing Table maintain in the local node.
    '''
    for j in xrange(0, len(Nodes[i].RREP_RoutingTable)):
        if Pkt.Dst == Nodes[i].RREP_RoutingTable[j][0]:
            k = Nodes[i].RREP_RoutingTable[j][2]
            if k in Nbrs and len(RQp[k]) <= Nodes[k].Buffer:

                NPkt =  Packet(Pkt.Id, Pkt.Src, Pkt.Dst, 1, Pkt.TTL - 1, Pkt.Num)
                Nodes[k].update_Packet_State(Pkt.Id, Pkt.Num, Pkt.Src)
                RQp[k].insert(0,NPkt)
        else:
            if Pkt.RouteRecovery > 0:
                RQp[i].insert(0,Pkt)
                # if the Packet couldn't find the routing information stored in the routing table, firt push the 
                # Packet in to the node's buffer,
                Pkt.RouteRecovery -= 1
                if Pkt.RouteRecovery == 4:
                    for Nbr in Nbrs:
                        if len(RQp[Nbr])<=Nodes[Nbr].Buffer:
                            NPkt =  RRER(Pkt.Id, Pkt.Src, Pkt.Dst, 3)
                            RQp[Nbr].insert(0,NPkt)
            else:
                NPkt = RREQ(Pkt.Id+1, Pkt.Src, Pkt.Dst, 30, 1)
                RQp[Pkt.Src].insert(0, NPkt)
                if (Pkt.Src, Pkt.Id) in Nodes[Pkt.Src].PktRemain:
                    Nodes[Pkt.Src].PktRemain[Pkt.Src, Pkt.Id] = 0
                # the Nodes won't generate the Packet.      
                    
                        
    return Nodes, RQp
     
def forward_RRER(Nodes, Adj, RQp, Pkt, Nbrs, i):
    Flag = 0
    for j in xrange(0, len(Nodes[i].RREP_RoutingTable)):
        if Pkt.Dst == Nodes[i].RREP_RoutingTable[j][0]:
            NPkt = RREP(Pkt.Id, Pkt.Src, Pkt.Dst, 5)
            RQp[i].insert(0, NPkt)     
            Flag = 1

    if Flag == 0:
        for Nbr in Nbrs:
            if len(RQp[Nbr]) <= Nodes[Nbr].Buffer:
                NPkt = RRER(Pkt.Id, Pkt.Src, Pkt.Dst, Pkt.TTL - 1)
                RQp[Nbr].insert(0, NPkt)

    return Nodes, RQp



def sim(Length=1000, Width=1000, Number=100, MaxIter = 100, TxRng = 50, RQNum = 20, InitTTL = 30):
    Nodes, RQ = init(Length, Width, Number)
    Nodes, RQ = run(Nodes, RQ, RQNum, MaxIter, TxRng, InitTTL)
    check_PktRecieve(Nodes)

def check_PktRecieve(Nodes):
    Num = 0
    for i in xrange(0, len(Nodes)):
        temp = []
        print Nodes[i].Recieve_Pkt
        for j in xrange(0, len(Nodes[i].Recieve_Pkt)):
            if [Nodes[i].Recieve_Pkt[j][0], Nodes[i].Recieve_Pkt[j][2]] not in temp:
                temp.append([Nodes[i].Recieve_Pkt[j][0], Nodes[i].Recieve_Pkt[j][2]])
        
        Num += len(temp)

    print Num                 
