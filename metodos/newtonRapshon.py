from sympy import sympify, symbols, lambdify, diff, sin, cos, exp
from sympy.utilities.lambdify import lambdify
import numpy as np
import matplotlib.pyplot as plt

def funcion(funcionTexto, x, derivada=False):
    simbolo = symbols('x')  # Declara la variable simbólica
    # Convertir el texto en una expresión simbólica con funciones reconocidas
    expresion = sympify(funcionTexto, locals={"sin": sin, "cos": cos, "exp": exp})  
    
    if derivada:
        expresion = diff(expresion, simbolo)  # Calcula la derivada
    
    # Convertir la expresión simbólica en una función evaluable
    func = lambdify(simbolo, expresion, "numpy")
    
    return func(x)  # Evalúa la función con el valor de x

def Ea(xNuevo,xAnterior):
  return abs(((xNuevo-xAnterior)/xNuevo))*100

def newtonRaphson(funcionTexto, xi, errorPorcentual):
  i = 0
  while(True):
    x = xi - (funcion(funcionTexto,xi)/funcion(funcionTexto,xi,derivada=True))
    E = Ea(x,xi)
    if i > 1 and E < (errorPorcentual/100):
      break
    xi = x
    i += 1
  return xi

def grafica(funcionTexto, xi):
  x = np.linspace(-2, 4, 100)
  y = funcion(funcionTexto, x)

  plt.plot(x, y, label= funcionTexto)
  plt.axhline(0, color='black', linewidth=0.5)
  plt.axvline(0, color='black', linewidth=0.5)

  yPunto = funcion(funcionTexto,xi)
  plt.scatter(xi, yPunto, color='black', zorder=5, label=f'{xi}')

  plt.xlim(-2, 4)
  plt.ylim(-10, 10)
  plt.title(f'Grafica de {funcionTexto}')
  plt.xlabel('x')
  plt.ylabel('f(x)')
  plt.legend()
  plt.grid(True)
  plt.show()



f = "(exp(-x))-x"
grafica(f,newtonRaphson(f,0,1))
