# QOSF-Task3

This program converts one quantum circuit into another, using a restricted set of gates (Rx, Rz, CZ).
The basic gates that can be present in the input circuit are: (I, H, X, Y, Z, Rx, Ry, Rz, CNOT, CZ).

## Working of the program

### 1. Finding the equivalent circuit for all input gates

To find the equivalent circuit for an input gate, I use the [DiVincenzo-Smolin scheme](https://arxiv.org/abs/cond-mat/9409111). This is an exhaustive approach that performs decomposition of a unitary matrix in SU(8) into a sequence of two-qubit gates (SU(4) matrices). Even though it focuses on three-qubit gates, it is extendable to unitary in arbitrary dimensions. The non-linear function that needs to be minimized is:

<img src="https://render.githubusercontent.com/render/math?math=f = \displaystyle\sum_{i}\displaystyle\sum_{j}|M_{ij} - U_{ij}|^{2}">

where M is the required matrix i.e. the matrix of the input gate and U is the matrix of the circuit containing Rx,Rz and CZ gates that we test against the required matrix. 

Our job is to minimize the function f numerically using any non-linear minimization tool. I have used Broyden-Fletcher-Goldfarb-Shanno algorithm (BFGS) optimization function. The successful optimization of the function, i.e. <img src="https://render.githubusercontent.com/render/math?math=|f|\leq\epsilon"> returns the required topology and the optimized measurement angles. Thus, the total operation from network gate U approximates M within error <img src="https://render.githubusercontent.com/render/math?math=\epsilon">.

#### Creation of U-gate
Now, I needed a U-gate composed of Rz, Rx and CZ gates only. I use the fact that all one qubit gates can be defined using the universal SU(2) unitary matrix upto a phase. I also make use of the fact that the SU(2) unitary gate can be broken down in terms of three rotation gates, i.e. <img src="https://render.githubusercontent.com/render/math?math=Rz(\alpha)Rx(\beta)Rz(\gamma)">. Thus, for all one qubit gates, my U can be <img src="https://render.githubusercontent.com/render/math?math=Rz(\alpha)"> or <img src="https://render.githubusercontent.com/render/math?math=Rz(\alpha)Rx(\beta)"> or <img src="https://render.githubusercontent.com/render/math?math=Rz(\alpha)Rx(\beta)Rz(\gamma)">. I find out the matrices for these combination of gates and test it against the required matrix using the optimization function. 

Similarly, for two-qubit gates, I just add a U-gate on either side of a CZ gate i.e. <img src="https://render.githubusercontent.com/render/math?math=(I\otimes U)-CZ-(I\otimes U)"> and test the matrix of the same against the matrix of the required two-qubit gate. 

###### Example
The Hadamard gate comes out to be Rz(pi/2)Rx(pi/2)Rz(pi/2)

I then save the results for the same in a dictionary which I add to a file (new_gates_dict.pickle). The result looks like:
```
{'I': [('Rz', 0.0)], 'H': [('Rz', 1.57), ('Rx', 1.57), ('Rz', 1.57)], 'X': [('Rz', 0.0), ('Rx', 3.13)], 'Y': [('Rz', 0.0), ('Rx', 3.12), ('Rz', 3.14)], 'Z': [('Rz', 3.14)], 'CNOT': [('Rz', 1.57), ('Rx', 1.57), ('Rz', 1.57), ('CZ',), ('Rz', 1.57), ('Rx', 1.57), ('Rz', 1.57)]}
```
### 2. Taking Input Circuit and Replacing Gates
The next step is to take in the input circuit with all the allowed basic gates of the form:
('I',q),('H',q),('X',q),('Y',q),('Z',q),('Rx',x,q),('Ry',x,q),('Rz',x,q),('CNOT',q1,q2),('CZ',q1,q2)
Here, q represents the qubit the gate acts on, x represents the angle for the rotational gates and q1,q2 represent control and target qubit respectively.

I then replace each input gate with the new configuration from the dictionary, adding in the qubit as well. I let the Rx, Rz and CZ gate configuartions stay the same as given in the input circuit and to replace the Ry gate, I run the optimization depending on the angle of Ry gate and then replace it. 

### 3. Reducing the overhead of the new circuit
I consider the depth of the circuit here. Since I am using a restricted gate set, the number of gates in the new circuit would be much more than the input circuit. To counter this issue, I optimize my new circuit as much as possible to reduce the number of gates. To do this I use the following reductions:

- Remove all Rx and Rz gates with angles 0 or 2pi since these result in an Identity gate which can be ignored
- Remove the occurences of two consecutive CZ gates acting on the same qubit pair
- Combine all Rx and Rz gates acting on the same qubits such that <img src="https://render.githubusercontent.com/render/math?math=Rn(\alpha)Rn(\beta)=Rn(\alpha%2B\beta)">

(I also move any Rz gates accross CZ gates acting on the same qubit as the Rz gates (as they are commutative) to have more opportunities of combining gates. For eg., if my configuration was ('Rz',x1,1),('CZ',0,1),('Rz',x2,1) then moving the Rz gate accross CZ would result in ('CZ',0,1),('Rz',x1,1),('Rz',x2,1). Then both the Rz gates can be combined.)

I do these reductions multiple times till I notice there is no more possibility of reductions. 


## An Example run on the Program
Input Ciruit:
```
  [('I', 0), ('H', 0), ('X', 0), ('Y', 1), ('Z', 1), ('Rx', 3.14, 0), ('Rz', 0.79, 0), ('CZ', 0, 1), ('CNOT', 1, 0), ('Ry', 1.57, 0)] 
Length: 10
```
Restricted Circuit without reduction:
```
  [('Rz', 0.0, 0), ('Rz', 1.57, 0), ('Rx', 1.57, 0), ('Rz', 1.57, 0), ('Rz', 0.0, 0), ('Rx', 3.13, 0), ('Rz', 0.0, 1), ('Rx', 3.12, 1), ('Rz', 3.14, 1), ('Rz', 3.14, 1), ('Rx', 3.14, 0), ('Rz', 0.79, 0), ('CZ', 0, 1), ('Rz', 1.57, 0), ('Rx', 1.57, 0), ('Rz', 1.57, 0), ('CZ', 1, 0), ('Rz', 1.57, 0), ('Rx', 1.57, 0), ('Rz', 1.57, 0), ('Rz', 4.71, 0), ('Rx', 1.57, 0), ('Rz', 1.57, 0)] 
Length: 23
```

Final restricted Circuit:
```
 [('Rz', 1.57, 0), ('Rx', 1.57, 0), ('Rz', 1.57, 0), ('Rx', 3.13, 0), ('Rx', 3.12, 1), ('Rx', 3.14, 0), ('CZ', 0, 1), ('Rz', 2.36, 0), ('Rx', 1.57, 0), ('CZ', 1, 0), ('Rz', 3.14, 0), ('Rx', 3.14, 0), ('Rz', 1.57, 0)] 
Length: 13
```

As shown here, the input circuit with 10 gates is changed to a circuit with restricted gate set. Without any reductions, the circuit is very involved with 23 gates but on reductions, we are left with only 13 gates which is only 3 more than out input circuit. 



















