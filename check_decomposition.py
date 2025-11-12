import numpy as np
import sys

if len(sys.argv) <= 1:
    print("Checks whether the projective decomposition of the noisy version of a POVM is valid.")
    print("")
    print("Usage: " + sys.argv[0] + " <filename>")
    print("where <filename> points to a .npy file containing the decomposition.")
    exit()

#Helper functions:

def dm(state):
    """ given a ket, returns its associated density matrix (outer product)
        args:      ndarray, vector/ket, e.g. ket('01')
        returns:    ndarray, density matrix
    """
    return np.outer(state, state.conj())


def povm_fiducial(psi):
    """ Constructs POVM from fiducial vector by applying Heisenberg Weyl operators
        psi: the (normalized) d-dimensional fiducial vector from which the POVM should be constructed
        returns: The POVM as list of d**2 effects M_i
    """
    d = len(psi)
    
    povm = []
    for i in range(d):
        for j in range(d):
            psi_ij = HW_X(d, i) @ HW_Z(d, j) @ psi
            povm.append(1/d * dm(psi_ij))
    
    return povm


def off_diagonal(d, u):
    """ d-dim matrix with offdiagonal at the u'th offdiagonal """
    return np.diagflat(np.ones(d-abs(u)), u)


def HW_X(D,k=1):
    """ k-th power of Heisenberg Weyl X (shift) matrix 
        Eq 15 in Scott, arXiv:quant-ph/0310137
    """
    return off_diagonal(D, -k) + off_diagonal(D, D-k)

def HW_Z(D, k=1):
    """ k-th power of Heisenberg Weyl Z (clock) matrix 
        Eq 15 in Scott, arXiv:quant-ph/0310137
    """
    return np.diag(np.exp(2j*np.pi*k*np.arange(D)/D) )


# In order to compare the simulation with the proper SIC POVM, define their fiducial vectors...
psi4 = np.array([0.48571221409126403909152153176812197109, 0.60043369656069688700611847041568366744-0.44989636690811813902417022753091663501j,-.20118858648686589293456281596678826706j,-0.39924511007383099407155565444889540038-0.035815847183145900067351304237205336087j])
psi5 = np.array([0.39104489402214774638257588694092612854,-0.28486558319586666154004262263615905414-0.64712933282796239400249892482493465752j,-0.23188384736899577826443055554623498413-0.19820390755555243362222174127453043647j,+0.13193857997561206409072011157833671545-0.10939599651964327439560846264470517498j,+0.43096743921096438598841830212227390108+0.19747405581954685663309013848038222203j])

# ... and build the corresponding SIC POVMs from it:
sic4 = povm_fiducial(psi4)
sic5 = povm_fiducial(psi5)

# Load npy file. It should contain the projective measurement decomposition of a POVM (M_1,...,M_d^2) in form of a dictionary. The keys are tuples indicating the row indices where the d projectors of the corresponding entry should contribute to. The values are lists of d matrices which should be projections, scaled by the probability to perform this measurement.

Ps = np.load(sys.argv[1], allow_pickle = True).item()
d = next(iter(Ps.values()))[0].shape[0]

num_measurements = len(Ps)

print("Loaded "+sys.argv[1]+", found a decomposition of a "+str(d)+"-dimensional POVM in terms of " +  str(num_measurements)+ " measurements.")

# To check if the decomposition really yields a noisy version of a SIC POVM, we calculate the row sums. Initialize them with zeros:
rowsum = {}
for i in range(d**2):
	rowsum[i] = np.zeros( (d,d), dtype = complex)


print("Checking if it really is a decomposition in terms of rank-1 projective measurements...")
# Iterate through all measurements..
for key in Ps:
    #key is a tuple indicating on which rows the projectors of this particular measurement should act.
    
    # The entries are scaled by a factor of num_measurements to avoid numerical precission problems. Undo this scaling
    Pk = [Ps[key][i] / num_measurements for i in range(len(Ps[key]))]
    
    # Calculate the probability of performing this measurement
    prob = np.trace(Pk[0]).real
    
    #Check if the measurement constitutes a proper PVM:
    
    for P in Pk:
        # 1. Check if matrices are rank-1
        assert( abs(np.linalg.eigvalsh(P)[-2]) < 1e-6)
        
        # 2. Check if matrices /after scaling by probability) are projections
        assert( np.linalg.norm( (P/prob) @ (P/prob) - P/prob)**2 < 1e-6 )
        
	# 3. Elements must sum to identity. To check that, form the sum:
    colsum = sum(Pk) / prob
    assert( np.linalg.norm(colsum - np.eye(d))**2 < 1e-6)
    
    # Add projectors to corresponding row sum to check similarity with noisy SIC POVM in the end.
    for i in range(len(Pk)):
        rowsum[key[i]] += Pk[i]
        
print("... it is!")
print("Checking if it sums to a noisy version of the SIC POVM...")
# Finally, check if projective simulation constitutes a noisy version of a SIC POVM. If that really is the case, then we can calculate the visibility from the smaller eigenvalue lambda_min via v = 1-d**2 * lambda_min

v = 1 - d**2 * np.linalg.eigvalsh(rowsum[0])[0]
print("... visibility would be v =",v)

sic = sic4 if d == 4 else sic5 if d == 5 else None
for i in range(d*d):
    Mnoisy = v * sic[i] + (1-v)/d**2 * np.eye(d)
    assert( np.linalg.norm(Mnoisy - rowsum[i]) < 1e-6)

print("... it is!")
print("Summary: The decomposition in "+sys.argv[1]+" constitutes a projective decomposition of a noisy version of a "+str(d)+"-dimensional SIC POVM with visibility parameter v =",v)
