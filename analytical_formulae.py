import numpy as np
import sympy as sym

psi = sym.Symbol('psi')	# telescope angle
phi = sym.Symbol('phi')	# HWP angle
deltaA = sym.Symbol('deltaA')	# ...
deltaB = sym.Symbol('deltaB')	# ...

def rot(theta):
    mat = sym.Matrix([[1,0,0],[0,sym.cos(2*theta),sym.sin(2*theta)],[0,-sym.sin(2*theta),sym.cos(2*theta)]])
    return mat
    
mat_det = sym.Matrix([[1,1,0],[1,1,0],[0,0,0]])/2	# ideal polarimeter

def A_row(psi,phi,delta,xi):
    rotat = rot(xi+delta+psi)
    det_rotat = mat_det*rotat
    return [det_rotat[0,0], det_rotat[0,1], det_rotat[0,2]]
    
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('one detector')
    
A_mat = sym.Matrix([A_row(0         ,phi,deltaA,0),
                    A_row(sym.pi/2  ,phi,deltaA,0),
                    A_row(sym.pi/4  ,phi,deltaA,0),
                    A_row(3*sym.pi/4,phi,deltaA,0)])

A_hat = sym.Matrix([A_row(0         ,phi,0,0),
                    A_row(sym.pi/2  ,phi,0,0),
                    A_row(sym.pi/4  ,phi,0,0),
                    A_row(3*sym.pi/4,phi,0,0)])

ATA = A_hat.transpose()*A_hat
print(ATA.inv()*A_hat.transpose())

ATAinvATA = ATA.inv()*A_hat.transpose()*A_mat

print('II: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,0]))))
print('IQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,1]))))
print('IU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,2]))))
print()              
print('QI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,0]))))
print('QQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,1]))))
print('QU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,2]))))
print()
print('UI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,0]))))
print('UQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,1]))))
print('UU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,2]))))

print('')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('two identical detectors')

A_mat = sym.Matrix([A_row(0         ,phi,deltaA,0),
                    A_row(sym.pi/2  ,phi,deltaA,0),
                    A_row(sym.pi/4  ,phi,deltaA,0),
                    A_row(3*sym.pi/4,phi,deltaA,0),
                    A_row(0         ,phi,-deltaA,0),
                    A_row(sym.pi/2  ,phi,-deltaA,0),
                    A_row(sym.pi/4  ,phi,-deltaA,0),
                    A_row(3*sym.pi/4,phi,-deltaA,0)])

A_hat = sym.Matrix([A_row(0         ,phi,0,0),
                    A_row(sym.pi/2  ,phi,0,0),
                    A_row(sym.pi/4  ,phi,0,0),
                    A_row(3*sym.pi/4,phi,0,0),
                    A_row(0         ,phi,0,0),
                    A_row(sym.pi/2  ,phi,0,0),
                    A_row(sym.pi/4  ,phi,0,0),
                    A_row(3*sym.pi/4,phi,0,0)])

ATA = A_hat.transpose()*A_hat
print(ATA.inv()*A_hat.transpose())

ATAinvATA = ATA.inv()*A_hat.transpose()*A_mat

print('II: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,0]))))
print('IQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,1]))))
print('IU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,2]))))
print()              
print('QI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,0]))))
print('QQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,1]))))
print('QU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,2]))))
print()
print('UI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,0]))))
print('UQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,1]))))
print('UU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,2]))))

print('')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('two different detectors')

A_mat = sym.Matrix([A_row(0         ,phi,deltaA,0),
                    A_row(sym.pi/2  ,phi,deltaA,0),
                    A_row(sym.pi/4  ,phi,deltaA,0),
                    A_row(3*sym.pi/4,phi,deltaA,0),
                    A_row(0         ,phi,-deltaB,0),
                    A_row(sym.pi/2  ,phi,-deltaB,0),
                    A_row(sym.pi/4  ,phi,-deltaB,0),
                    A_row(3*sym.pi/4,phi,-deltaB,0)])

A_hat = sym.Matrix([A_row(0         ,phi,0,0),
                    A_row(sym.pi/2  ,phi,0,0),
                    A_row(sym.pi/4  ,phi,0,0),
                    A_row(3*sym.pi/4,phi,0,0),
                    A_row(0         ,phi,0,0),
                    A_row(sym.pi/2  ,phi,0,0),
                    A_row(sym.pi/4  ,phi,0,0),
                    A_row(3*sym.pi/4,phi,0,0)])

ATA = A_hat.transpose()*A_hat
print(ATA.inv()*A_hat.transpose())

ATAinvATA = ATA.inv()*A_hat.transpose()*A_mat

print('II: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,0]))))
print('IQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,1]))))
print('IU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[0,2]))))
print()              
print('QI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,0]))))
print('QQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,1]))))
print('QU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[1,2]))))
print()
print('UI: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,0]))))
print('UQ: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,1]))))
print('UU: '+ str(sym.simplify(sym.trigsimp(ATAinvATA[2,2]))))
