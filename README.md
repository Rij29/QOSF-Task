# QOSF-Task3

This program converts one quantum circuit into another, using a restricted set of gates (Rx, Rz, CZ).
The basic gates that can be present in the input circuit are: (I, H, X, Y, Z, Rx, Ry, Rz, CNOT, CZ).

## Working of the program

### 1. Finding the equivalent circuit for all input gates











The DiVincenzo-Smolin scheme [55] is a systematic, exhaustive approach, that performs a decomposition of a unitary matrix in SU(8) into a sequence of two-qubit gates (SU(4) matrices) that is called a two-bit gate network. The scheme focuses on the decomposition of three-qubit gates, but it is extendable to unitary in arbitrary dimensions, viz., SU(2n) for n ∈ N. It had been believed that a three-qubit gate is necessary for universality, implied by the result of Deutsch [164], until some years later, DiVincenzo proved that two-qubit gates are also universal [16]. Moreover, in [164], the efficiency (in the gate-depth sense) of implementing a SU(2n) matrix is not addressed. The DiVincenzo-Smolin scheme addressed Deutsch’s issues by providing an explicit implementation to obtain a two-qubit network of a SU(8) matrix. Here, we provide an extension of the DiVincenzo-Smolin scheme that works for an arbitrary dimensional unitary.
The gist of DiVincenzo-Smolin’s scheme lies in the non-linear minimization pro- cedure of a gate network. A gate network is characterized by its two-bit gate topology that is defined in Definition 6.1. Now, let us discuss gate topologies in more detail.







The DiVincenzo-Smolin scheme in Algorithm 6.1 performs a decomposition of an SU(2n) unitary matrix M, sought within combinations of N two-qubit gates, and
3
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






















