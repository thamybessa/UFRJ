#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 2017

@author: thamybessa
"""
from __future__ import print_function
from __future__ import division
import math

def confere(n, m): #função de simples conferência do resultado... formulinha de matcomb
	soma = n + m
	resultado = math.factorial(soma)/(math.factorial(n)*math.factorial(m))
	print(resultado)

	n = raw_input("digite um numero para a linha da matriz: ")
	n = int(n)
	m = raw_input("digite um numero para a coluna da matriz: ")
	m = int(m)
#confere(n, m)

def computa_caminhos(n): 
	matriz = [[0 for i in range(n+1)] for j in range(n+1)]
	matriz[1][1] = 2
	for i in range(2, n+1): #prepara caso base
		matriz[1][i] = (matriz[1][i-1]) + 1
		matriz[i][1] = (matriz[i-1][1]) + 1

	for i in range(2, n+1): #percorro cada linha de cada coluna ##colunas
		for j in range(2, n+1): ##linhas
			matriz[i][j] = matriz[i][j-1] + matriz[i-1][j]
	print("O número de caminhos possíveis é: ")		
	print(matriz[n][n])

n = raw_input("digite um numero n para a matriz quadrada: ")
n = int(n)
computa_caminhos(n)