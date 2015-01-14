# This script generates c++-code for the three particle force in MD silica code.

from sympy.utilities.codegen import codegen
from sympy import *

xij, yij, zij, xik, yik, zik = symbols('xij, yij, zij, xik, yik, zik')
B, costheta0, r0a, r0b, ksia, ksib = symbols('B, cosTheta0, r0a, r0b, ksia, ksib')
rij, rik, rijDotRik, cosTheta = symbols('rij, rik, rijDotRik, cosTheta')

V = B*exp(ksia / (sqrt(xij**2 + yij**2 + zij**2) - r0a) + ksib/( sqrt(xik**2 + yik**2 + zik**2) - r0b))*( (xij*xik + yij*yik + zij*zik) / ( sqrt(xij**2 + yij**2 + zij**2)*sqrt(xik**2 + yik**2 + zik**2)) - costheta0)**2

dVdXij = diff(V, xij)
dVdXij = dVdXij.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdXij = dVdXij.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdXij = dVdXij.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdXij = dVdXij.subs(rijDotRik/(rij*rik), cosTheta)

dVdYij = diff(V, yij)
dVdYij = dVdYij.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdYij = dVdYij.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdYij = dVdYij.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdYij = dVdYij.subs(rijDotRik/(rij*rik), cosTheta)

dVdZij = diff(V, zij)
dVdZij = dVdZij.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdZij = dVdZij.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdZij = dVdZij.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdZij = dVdZij.subs(rijDotRik/(rij*rik), cosTheta)

dVdXik = diff(V, xik)
dVdXik = dVdXik.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdXik = dVdXik.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdXik = dVdXik.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdXik = dVdXik.subs(rijDotRik/(rij*rik), cosTheta)

dVdYik = diff(V, yik)
dVdYik = dVdYik.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdYik = dVdYik.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdYik = dVdYik.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdYik = dVdYik.subs(rijDotRik/(rij*rik), cosTheta)

dVdZik = diff(V, zik)
dVdZik = dVdZik.subs(sqrt(xij**2 + yij**2 + zij**2), rij)
dVdZik = dVdZik.subs(sqrt(xik**2 + yik**2 + zik**2), rik)
dVdZik = dVdZik.subs(xij*xik + yij*yik + zij*zik, rijDotRik)
dVdZik = dVdZik.subs(rijDotRik/(rij*rik), cosTheta)

print codegen(("dVdXij", simplify(dVdXij)), "C", "file")[0][1]
print codegen(("dVdYij", simplify(dVdYij)), "C", "file")[0][1]
print codegen(("dVdZij", simplify(dVdZij)), "C", "file")[0][1]

print codegen(("dVdXik", simplify(dVdXik)), "C", "file")[0][1]
print codegen(("dVdYik", simplify(dVdYik)), "C", "file")[0][1]
print codegen(("dVdZik", simplify(dVdZik)), "C", "file")[0][1]