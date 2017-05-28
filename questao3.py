#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 13:15:32 2017

@author: thamybessa & brunohry
"""
from __future__ import print_function
import numpy as np
import numpy.linalg as npl
import scipy.linalg as spl
import scipy as sp

def Separador ():
    print(" ")
    print(42*"-")
    print(" ")

x_dif = np.array([2,3,5,8,13,21]) #estilo fibonacci 
x_igual = np.array([2,4,6,8,10,12]) #p.a. de razao 2

def acha_yn(xzinho):
	yn = []
	for numero in xzinho:
		y = np.sin(numero)
		yn.append(y)
	return yn

y_dif = acha_yn(x_dif) #calculando os yn para x com expaçamento diferente
y_igual = acha_yn(x_igual) #calculando os yn para x com expaçamento igual

#trazendo os pares
def pares_funcao(x=[], y=[]):
    pares = []
    for itemx, itemy in zip(x, y):
        pares.append([itemx, itemy])
    print("Os pares da função são: ")
    print(pares)
    return pares
def lusempivot(A, ptol, y): #decomposicao LU sem pivot

    m,n = np.shape(A)
    for i in np.arange(0,n):
        pivot = A[i,i]
        if abs(pivot) < ptol:
           print("nenhum pivot encontrado")
           break
        for k in np.arange(i+1,n):
            A[k,i] = A[k,i]/pivot
            A[k,i+1:n] = A[k,i+1:n] - A[k,i]*A[i,i+1:n]
    Lnp = np.eye(n)+np.tril(A,-1)
    Unp = np.triu(A)
    yzinho_np = npl.solve(Lnp, y)
    solucao_np = npl.solve(Unp, yzinho_np)
    # print(Lnp)
    # print(Unp)
    print(solucao_np)
    return Lnp,Unp, yzinho_np, solucao_np

def calcula_polinomio(solucao, arrayx):
	expoente = 5
	solu = solucao.tolist()
	x_array = arrayx.tolist()
	p1 = 0
	p2 = 0
	p3 = 0
	p4 = 0
	p5 = 0
	p6 = 0
	polinomio = []
	for i in solu:
		p1 += i*(x_array[0]**(expoente))
		p2 += i*(x_array[1]**(expoente))
		p3 += i*(x_array[2]**(expoente))
		p4 += i*(x_array[3]**(expoente))
		p5 += i*(x_array[4]**(expoente))
		p6 += i*(x_array[5]**(expoente))
		expoente -= 1
	polinomio = [p1, p2, p3, p4, p5, p6]
	print(polinomio)


Vander_v = np.vander(x_dif) #matriz de Vandermonde para x com espaçamento diferente
Vander_u = np.vander(x_igual) #matriz de Vandermonde para x com espaçamento diferente
cond_v = npl.cond(Vander_v) #numero de condicionamento da matriz de Vandermonde para x com espaçamento diferente
cond_u = npl.cond(Vander_u)	#numero de condicionamento da matriz de Vandermonde para x com espaçamento igual
epsilon = np.finfo(dtype=np.float64).eps
digitos_precisos_v = -np.round(np.log10(npl.cond(Vander_v)*epsilon),1) #digitos precisos da matriz v
digitos_precisos_u = -np.round(np.log10(npl.cond(Vander_u)*epsilon),1) #digitos precisos da matriz u

coeficientes_de_V = np.linalg.solve(Vander_v, y_dif) #resolução do sistema linear para matriz V
coeficientes_de_U = x = np.linalg.solve(Vander_u, y_igual) #resolução do sistema linear para matriz U

#Solucao por decomposição LU com Pivoteamento
Pv, Lv, Uv = spl.lu(Vander_v) #decompondo a matriz V
Pu, Lu, Uu = spl.lu(Vander_u) #decompondo a matriz U
yzinho_v = np.linalg.solve(Lv, y_dif) #computando a primeira parte do resultado de V
yzinho_u = np.linalg.solve(Lu, y_igual) #computando a primeira parte do resultado de U
Solucao_LU_V = np.linalg.solve(Uv, yzinho_v) #solução para V
Solucao_LU_U = np.linalg.solve(Uu, yzinho_u) #solução para U

#Solução por Decomposição QR com pivoteamento
a = spl.qr(Vander_v, pivoting=True) #decompondo a matriz V
Q_v = a[0]
R_v = a[1]
y_v = sp.dot(Q_v.T, y_dif)
xQR_v = np.linalg.solve(R_v,y_v)
b = spl.qr(Vander_u, pivoting=True) #decompondo a matriz U
Q_u = b[0]
R_u = b[1]
y_u = sp.dot(Q_u.T, y_igual)
xQR_u = np.linalg.solve(R_u,y_u)

#Solução por Decomposição QR sem pivoteamento
Qnp_v, Rnp_v = spl.qr(Vander_v) #decompondo a matriz V
ynp_v = sp.dot(Qnp_v.T, y_dif)
xQRnp_v = np.linalg.solve(Rnp_v, ynp_v)
Qnp_u, Rnp_u = spl.qr(Vander_u) #decompondo a matriz U
ynp_u = sp.dot(Qnp_u.T, y_igual)
xQRnp_u = np.linalg.solve(Rnp_u, ynp_u)

#Solução por Decomposição Polar
u_polar_v,p_polar_v = spl.polar(Vander_v) #decompondo a matriz V
Y_polar_v = np.linalg.solve(u_polar_v, y_dif)
sol_polar_v = np.linalg.solve(p_polar_v, Y_polar_v) #Solucao para V
u_polar_u,p_polar_u = spl.polar(Vander_u) #decompondo a matriz U
Y_polar_u = np.linalg.solve(u_polar_u, y_igual)
sol_polar_u = np.linalg.solve(p_polar_u, Y_polar_u) #solucao para U


tolerancia = 50*epsilon

print("Os pares da funcao de v:")
pares_funcao(x_dif, y_dif)
print("A matriz vandermonde de v:")
print(Vander_v)
print("E seu número de condicionamento:")
print(cond_v)
print("E quantos dígitos de precisão:")
print(digitos_precisos_v)
print("Coeficientes da interpoladora:")
print(coeficientes_de_V)
print("Conferindo...:")
print("y da função original:")
print(y_dif)
print("y a partir da solução:")
calcula_polinomio(coeficientes_de_V, x_dif)
Separador()
print("Solução da matriz V com decomposição LU com pivoteamento:")
print(Solucao_LU_V)
print("y da função original:")
print(y_dif)
print("y a partir da solução:")
calcula_polinomio(Solucao_LU_V, x_dif)

Separador()
print("Solução da matriz V com decomposição LU sem pivoteamento:")
sol_np_V = lusempivot(Vander_v, tolerancia, y_dif)
Separador()
print("Solução da matriz V com decomposição QR sem pivoteamento:")
print(xQRnp_v)
print("Conferindo...:")
print("y da função original:")
print(y_dif)
print("y a partir da solução:")
calcula_polinomio(xQRnp_v, x_dif)
Separador()
print("Solução da matriz V com decomposição QR com pivoteamento:")
print(xQR_v)
print("Conferindo...:")
print("y da função original:")
print(y_dif)
print("y a partir da solução:")
calcula_polinomio(xQR_v, x_dif)
print("Solução da matriz V com decomposição Polar:")
print(sol_polar_v)
print("Conferindo...:")
print("y da função original:")
print(y_dif)
print("y a partir da solução:")
calcula_polinomio(sol_polar_v, x_dif)
Separador()
Separador()



print("Os pares da funcao de u:")
pares_funcao(x_igual, y_igual)
print("A matriz vandermonde de u:")
print(Vander_u)
print("E seu número de condicionamento:")
print(cond_u)
print("E quantos dígitos de precisão:")
print(digitos_precisos_u)
print("Coeficientes da interpoladora:")
print(coeficientes_de_U)
print("Conferindo...:")
print("y da função original:")
print(y_igual)
print("y a partir da solução:")
calcula_polinomio(coeficientes_de_U, x_igual)
Separador()
print("Solução da matriz U com decomposição LU com pivoteamento:")
print(Solucao_LU_U)
print("Conferindo...:")
print("y da função original:")
print(y_igual)
print("y a partir da solução:")
calcula_polinomio(Solucao_LU_U, x_igual)
Separador()
print("Solução da matriz U com decomposição LU sem pivoteamento:")
sol_np_U = lusempivot(Vander_u, tolerancia, y_igual)
Separador()
print("Solução da matriz U com decomposição QR sem pivoteamento:")
print(xQRnp_u)
print("Conferindo...:")
print("y da função original:")
print(y_igual)
print("y a partir da solução:")
calcula_polinomio(xQRnp_u, x_igual)
Separador()
print("Solução da matriz U com decomposição QR com pivoteamento:")
print(xQR_u)
print("Conferindo...:")
print("y da função original:")
print(y_igual)
print("y a partir da solução:")
calcula_polinomio(xQR_u, x_igual)
Separador()
print("Solução da matriz U com decomposição Polar:")
print(sol_polar_u)
print("Conferindo...:")
print("y da função original:")
print(y_igual)
print("y a partir da solução:")
calcula_polinomio(sol_polar_u, x_igual)
Separador()



