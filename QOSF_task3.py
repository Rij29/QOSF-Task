#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__doc__=""" QOSF_task3.py: It takes an input quantum circuit and changes the gates
                            to only the set(Rx,Ry,CZ)

    USAGE: 

        Run QOSF_task3.py
        
"""
"""
Created on Fri Sep 25 19:48:07 2020

@author: Rijul Sachdeva
"""

import numpy as np
from scipy.optimize import minimize
import math
import json,pickle


########### Defining all quantum gates ###############
def I():
    """Identity gate"""
    return (np.matrix([[1,0],[0,1]]), 'I')
   
def H():
    """Hadamard gate"""
    return (np.matrix([[1,1],[1,-1]])/np.sqrt(2),'H')

def X():
    """Pauli X-gate"""
    return np.matrix([[0,1],[1,0]]),'X'

def Y():
    """Pauli Y-gate"""
    return np.matrix([[0,-1j],[1j,0]]),'Y'

def Z():
    """Pauli Z-gate"""
    return np.matrix([[1,0],[0,-1]]),'Z'

def Rx(x):
    """Rotational X-gate with angle x"""
    return np.matrix([[math.cos(x/2),-1j*math.sin(x/2)],[-1j*math.sin(x/2),math.cos(x/2)]])

def Ry(x):
    """Rotational Y-gate with angle x"""
    return np.matrix([[math.cos(x/2),-math.sin(x/2)],[math.sin(x/2),math.cos(x/2)]])

def Rz(x):
    """Rotational Z-gate with angle x"""
    return np.matrix([[np.exp(-1j*x/2),0],[0,np.exp(1j*x/2)]])

def CNOT():
    """The CNOT gate"""
    return np.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

def CZ():
    """The C-phase gate"""
    return np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])



########## This is the circuit template that contains only Rx, Rz and CZ gates ##########

def output_template1(x,j):
    
    """
    It returns a circuit template with a combination of Rz and Rx gates to test against
        the circuit
    
    param:
        :x: This is the list of angles for the Rz and Rx gates
        
    """
    #p0,p1 = np.pi*x[0],np.pi*x[1]
    p0,p1 = x[0],x[1]
    
    if(j==1):
        U = np.exp(1j*p0)*Rz(p1)
    elif(j==2):
        #p2 = np.pi*x[2]
        p2 = x[2]
        
        U = np.exp(1j*p0)*(Rx(p2)@Rz(p1))
    elif (j==3):
        #p2, p3 = np.pi*x[2], np.pi*x[3]
        p2, p3 = x[2], x[3]
        
        U = np.exp(1j*p0)*(Rz(p3)@Rx(p2)@Rz(p1))
    return U

def kron(M):
    i = I()
    return np.kron(i[0],M)

def output_template2(x,j):
    """
      It returns a circuit template with a combination of Rz, Rx and CZ gates to test against
        the circuit
    
    param:
        :x: This is the list of angles for the Rz and Rx gates
        
    """
   # p0,p1,p2 = np.pi*x[0],np.pi*x[1],np.pi*x[2]
    p0,p1,p2 = x[0],x[1],x[2]
    
    if(j==1):
        a = kron(Rz(p2)) @ CZ() @ kron(Rz(p1))
        U = np.exp(1j*p0)*a
    
    if(j==2):
        #p3,p4 = np.pi*x[3], np.pi*x[4]
        p3,p4 = x[3], x[4]
        
        a = kron(Rx(p4)) @ kron(Rz(p3))
        b = kron(Rx(p2)) @ kron(Rz(p1))
        c = a @ CZ() @ b
        U = np.exp(1j*p0)*c
        
    if(j==3):
        #p3, p4, p5, p6 = np.pi*x[3], np.pi*x[4], np.pi*x[5], np.pi*x[6]
        p3, p4, p5, p6 = x[3], x[4], x[5], x[6]
        
        a = kron(Rz(p6)) @ kron(Rx(p5)) @ kron(Rz(p4))
        b = kron(Rz(p3)) @ kron(Rx(p2)) @ kron(Rz(p1))
        c = a @ CZ() @ b
        U = np.exp(1j*p0)*c
        
    return U
    
     
    
def mprint(M, dig=2, w=None):
    """ 
    Print matrix M with significant digit 2 
    
    param:
    :M: matrix to be printed
    :dig: significant digit
    :w: width of the whole number
    """
    # printing parameters 
    row = len(M)
    col = int(np.size(M)/row)

    if w is None : 
        a = M[0,0] #sample
        if isinstance(a, complex):
            w = 8+dig*2
        else : 
            w = 4+dig
    # print the matrix 
    fmt = "{0:"+str(w)+"."+str(dig)+"f}"
    for i in range(row):
        for j in range(col):
            print(fmt.format(M[i,j]), end='')
        print("\n",end='')
        
        
        
def objective_function(x,M,j,qubitgate):
    """The function which needs to be minimized
    
        param:
        :x: The optimization parameter (angles of U)
        :M: The required matrix
        :qubitgate: The number of qubits the gate acts on(Can be 1 or 2)
        """
    objective_val = 0.0
     
    if(qubitgate==1):
        U = output_template1(x,j)
    elif(qubitgate==2):
        U = output_template2(x,j)
        
    dim = len(M)
    
    for i in range(dim):
        for j in range(dim):
             if M[i,j] != 0.0 :
                 objective_val += np.abs(M[i,j] - U[i,j])**2
                
    return objective_val


def optimization_func(M,qubitgate,count=0): 
    """
    This function optimizes the circuit
    
    param:
        :M: np.array, the desired unitary
        :qubitgate: The number of qubits the gate acts on(Can be 1 or 2) 
    """
    
    j=1
    angles = []
    while j<4:
        succ_run = False 
        while not succ_run :
            if(qubitgate == 1):
                b=2*np.pi
                x = np.random.uniform(low=0.0,high=b,size=j+1)
                bounds = [(0,b)] *(j+1)
                result = minimize(objective_function, x, args=(M,j,qubitgate), method='L-BFGS-B', bounds=bounds, options={'gtol': 1e-06})
                succ_run = result.success
                #print('j=',j,', success:',succ_run)
            elif(qubitgate==2):
                b=2*np.pi
                x = np.random.uniform(low=0.0,high=b,size=(2*j)+1)
                #x = [1]+[0.5]*(2*j)
                bounds = [(0,b)] *((2*j)+1)
                result = minimize(objective_function, x, args=(M,j,qubitgate), method='TNC', bounds=bounds, options={'gtol': 1e-06})
                succ_run = result.success
                #print('j=',j,', success:',succ_run)
        if(math.floor(result.fun*100)==0):
                #print('fun:',result.fun)
                angles = result.x.tolist()
                angles = [round(x,2) for x in angles]
                #print('angles', angles)
                if(qubitgate==1):
                    a = output_template1(angles,j)
                if(qubitgate==2):
                    a = output_template2(angles,j)
                #mprint(a)
                break
        j+=1
    if(len(angles)!=0):
        return angles
    elif(count<10):
        return optimization_func(M,qubitgate,count+1)
    


def new_configuration(gate,qubitgate):
    """
     This returns the new configuration of the input gate in terms of Rz, Rx and CZ gates
     
     param:
         :gate: This is the gate to be optimized
         :qubitgate: The number of qubits the gate acts on(Can be 1 or 2)
    """
   
    angles = optimization_func(gate,qubitgate)
    
    if(qubitgate==1):
        if (len(angles)==2):
            return [('Rz',angles[1])]
        elif(len(angles)==3):
            return [('Rz',angles[1]),('Rx',angles[2])]
        elif(len(angles)==4):
            return [('Rz',angles[1]),('Rx',angles[2]),('Rz',angles[3])]

    elif(qubitgate==2):
        if(len(angles)==3):
            return [('Rz',angles[1]),('CZ',),('Rz',angles[2])]
        elif(len(angles)==5):
            return [('Rz',angles[1]),('Rx',angles[2]),('CZ',),('Rz',angles[3]),('Rx',angles[4])]
        elif(len(angles)==7):
            return [('Rz',angles[1]),('Rx',angles[2]),('Rz',angles[3]),('CZ',),('Rz',angles[4]),('Rx',angles[5]),('Rz',angles[5])]


def reducing_overhead(new_circuit):
    """
    This reduces the number of gates in the final circuit
    
    param:
        :new_circuit: This is the list of final gates without being reduced
    """
    
    condn = True
    while condn:
        start = len(new_circuit)
       
        for i in new_circuit:
            if ((i[0]=='Rx' or i[0]=='Rz') and (i[1]==0 or 6.25<=i[1]<=round(2*np.pi,2))):
                new_circuit.remove(i)
                
        #print('\n',new_circuit, '\nLength:',len(new_circuit))
        
        for i in range(len(new_circuit)-1):
            if(new_circuit[i][0]=='Rz' and new_circuit[i+1][0]=='CZ'):
                if(new_circuit[i][2]==new_circuit[i+1][1] or new_circuit[i][2]==new_circuit[i+1][2]):
                    new_circuit[i],new_circuit[i+1] = new_circuit[i+1],new_circuit[i]
        
        #print('\n',new_circuit, '\nLength:',len(new_circuit))
        
        i=0
        while i<len(new_circuit)-1:   
            if(new_circuit[i][0]=='Rz' and new_circuit[i+1][0]=='Rz'):
                if(new_circuit[i][2]==new_circuit[i+1][2]):
                    l = list(new_circuit[i])
                    summ = new_circuit[i][1] + new_circuit[i+1][1]
                    if(summ > 6.28):
                        q = summ//(2*np.pi)
                        a = summ - (q*2*np.pi)
                        l[1] = round(a,2)
                    else:
                        l[1] = round(summ,2)
                    new_circuit[i] = tuple(l)
                    del new_circuit[i+1]
                    
            if(new_circuit[i][0]=='Rx' and new_circuit[i+1][0]=='Rx'):
                if(new_circuit[i][2]==new_circuit[i+1][2]):
                    l = list(new_circuit[i])
                    summ = new_circuit[i][1] + new_circuit[i+1][1]
                    if(summ > 6.28):
                        q = summ//(2*np.pi)
                        a = summ - (q*2*np.pi)
                        l[1] = round(a,2)
                    else:
                        l[1] = round(summ,2)
                    new_circuit[i] = tuple(l)
                    del new_circuit[i+1]
            
            if(new_circuit[i][0]=='CZ' and new_circuit[i+1][0]=='CZ'):
                if(new_circuit[i][1]==new_circuit[i+1][1] and new_circuit[i][2]==new_circuit[i+1][2]):
                    del new_circuit[i+1]
                    
            i+=1
        end = len(new_circuit)
        if(end>=start):
            condn = False    
    return new_circuit
    


def output_circuit(new_gates,input_circuit):
    """
    This gives the new circuit containing only Rz, Rx and CZ gates
    It further calls the reducing_overhead function to get the final circuit after removal of 
        unnecesary gates
    
    param:
        :new_gates: dictionary that contains the new configuration for input gates
        :input_circuit: This is the input circuit

    """
    
    keys = list(new_gates.keys())
    new_circuit = []
    
    for i in input_circuit:
        if((i[0] in keys) and (i[0]!='CNOT')):
            new = [elem + (i[1],) for elem in new_gates[i[0]]]
            new_circuit.extend(new)
        elif((i[0] in keys) and (i[0]=='CNOT')):
            new = []
            for j in new_gates[i[0]]:
                if(j[0]!='CZ'):
                    new.append(j+(i[2],))
                elif(j[0]=='CZ'):
                    new.append(j + (i[1],)+(i[2],))     
            new_circuit.extend(new)
        elif(i[0]=='Rz' or i[0]=='Rx' or i[0]=='CZ'):
            new_circuit.append(i)
        elif(i[0]=='Ry'):
            prenew = new_configuration(Ry(i[1]),1)
            new = [elem + (i[2],) for elem in prenew]
            new_circuit.extend(new)
    
    print('\nOutput Circuit without reduction:\n ', new_circuit, '\nLength:',len(new_circuit))
    return reducing_overhead(new_circuit)




def main():   
    
    # gates_one_qubit = [I(),H(),X(),Y(),Z()]
    # gates_two_qubit = [CNOT()]
    # new_gates = dict()
    
    # for i in range(len(gates_one_qubit)):
    #     g = gates_one_qubit[i]
    #     new_gates[g[1]] = new_configuration(g[0],1)
    
    # #new_gates['CNOT'] = new_configuration(CNOT(), 2)
    # new_gates['CNOT'] =   [('Rz',1.57),('Rx',1.57),('Rz',1.57),('CZ',),('Rz',1.57),('Rx',1.57),('Rz',1.57)]
            
    # print('dict:',new_gates)
    
    
    # with open('new_gates_dict.pickle','wb') as outfile:
    #     pickle.dump(new_gates, outfile)
        


    with open('new_gates_dict.pickle', 'rb') as f:
          new_gates = pickle.load(f)
    
    #print(new_gates)
    
    input_circuit = [('I',0),('H',0),('X',0),('Y',1),('Z',1),('Rx',round(np.pi,2),0),('Rz',round(np.pi/4,2),0),('CZ',0,1),('CNOT',1,0),('Ry',round(np.pi/2,2),0)]
    print('Input Ciruit:\n ', input_circuit,'\nLength:',len(input_circuit))
    
    new_circuit = output_circuit(new_gates, input_circuit)  
    print('\nFinal output Circuit:\n',new_circuit,'\nLength:',len(new_circuit))    
    
    

if __name__ == "__main__" :
    main()
    