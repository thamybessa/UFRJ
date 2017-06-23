#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 2017

@author: thamybessa
"""
from __future__ import print_function

def N_bits(n):
	array = [None]*(n+1)
	array[1] = 2
	array[2] = 3
	for i in range(3, n+1):
		k = i-1
		j = i-2
		array[i] = (array[k] + array[j])

	print("O número de números é:")
	print(array[n])

numero = raw_input("Digite um número inteiro para obter o resultado: ")
numero = int(numero)
N_bits(numero)