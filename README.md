# QOSF-Task3

This program converts one quantum circuit into another, using a restricted set of gates (Rx, Rz, CZ).
The basic gates that can be present in the input circuit are: (I, H, X, Y, Z, Rx, Ry, Rz, CNOT, CZ).

## Working of the program

### 1. Finding the equivalent circuit for all input gates

To find the equivalent circuit for an input gate, I use the [DiVincenzo-Smolin scheme](https://arxiv.org/abs/cond-mat/9409111). This is an exhaustive approach that performs decomposition of a unitary matrix in SU(8) into a sequence of two-qubit gates (SU(4) matrices). Even though it focuses on three-qubit gates, it is extendable to unitary in arbitrary dimensions. The non-linear function that needs to be minimized is:

<img src="https://render.githubusercontent.com/render/math?math=f = \sum_{i}\sum_{j}|M_{ij} - S_{ij}|^{2}">














The DiVincenzo-Smolin scheme in Algorithm 6.1 performs a decomposition of an SU(2n) unitary matrix M, sought within combinations of N two-qubit gates, and
ones; this can be done by a particular search scheme, e.g., Algorithm 6.2. One may choose a suitable non-linear minimization scheme. For instance, DiVincenzo-Smolin scheme employs the Broyden-Fletcher-Goldfarb-Shanno algorithm (BFGS) [165] with the objective function defined as
f :=􏰜􏰜|Mij −Sij|2, (6.2) ij
⃗⃗ where S is the total resulting matrix from gate network N(Φ), where N(Φ) is
constructed from T according to gates S. The successful optimization, i.e., |f| ≤ ε, ⃗
operation from network gate N(Φ) approximates M within error ε.
The DiVincenzo-Smolin scheme allows one to obtain arbitrary precision of de- composition. The scheme works particularly well for a small depth with a small number of qubits, as the algorithm grows exponentially in the depth. Moreover, a larger size topology exhausts the optimization, as the free parameter grows linearly in the depth, which quadratically increases the number of elements of the Hessian matrix.
3 The value of ε is determined by the outcome of the chosen objective function, i.e., the Euclidean distance matrix like f. For this case, the best ε will be the machine precision, e.g., the smallest floating point representable in Python is 2.225 × 10−308 with 53 digits precision.
All of possible topologies are captured in variable T, where |T| ≤
within precision ε.
􏰄n􏰅N , namely less than all possible topologies before eliminating the equivalent






















