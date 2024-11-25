from sympy import sympify, symbols, lambdify, diff, sin, cos, exp
from sympy.utilities.lambdify import lambdify
import numpy as np

def funcion(funcionTexto, x, derivada=False):
    simbolo = symbols('x')  # Declara la variable simbólica
    # Convertir el texto en una expresión simbólica con funciones reconocidas
    expresion = sympify(funcionTexto, locals={"sin": sin, "cos": cos, "exp": exp})  
    
    if derivada:
        expresion = diff(expresion, simbolo)  # Calcula la derivada
    
    # Convertir la expresión simbólica en una función evaluable
    func = lambdify(simbolo, expresion, "numpy")
    
    return func(x)  # Evalúa la función con el valor de x

f = "-(exp(-x))-1"
xi=0

numerador = funcion(f,xi)
denominador = funcion(f,xi,derivada=True)
x = xi - (numerador/denominador)

print(numerador)
print(denominador)

print(exp(-0))