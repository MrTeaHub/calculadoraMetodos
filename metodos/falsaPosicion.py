import numpy as np
import matplotlib.pyplot as plt
from sympy import sympify, symbols, lambdify, diff, sin, cos, exp, sqrt

# Jesús David Peña Chivatá 160004828

# f(x)= x**10 - 1
# Error 0.1%
# [0 - 1.3]

def funcion(funcionTexto, x, derivada=False):
    simbolo = symbols('x')  
    # Convertir el texto en una expresión simbólica con funciones reconocidas
    expresion = sympify(funcionTexto, locals={"sin": sin, "cos": cos, "exp": exp, "sqrt":sqrt})  
    
    if derivada:
        expresion = diff(expresion, simbolo) 
    func = lambdify(simbolo, expresion, "numpy")
    return func(x)  # Evalúa la función con el valor de x

def Ea(xNuevo, xAnterior):
  return ((xNuevo-xAnterior)/xNuevo)*100

def falsaPosicion(funcionTexto, xl, xu, errorPorcentual):
    xrAnterior = 0
    rta = []
    while True:
        xr = xu-((funcion(funcionTexto, xu)*(xl-xu))/(funcion(funcionTexto, xl) -funcion(funcionTexto, xu)))
        x = round(funcion(funcionTexto, xl)*funcion(funcionTexto, xr),8)
        e = Ea(xr,xrAnterior)
        if x < 0:
            xu = xr
        elif x > 0:
            xl = xr
        elif x == 0 and e < (errorPorcentual/100):
            rta.append(e)
            rta.append(xr)
            break
        xrAnterior=xr
    return rta

def grafica(f,xPuntos):
  x = np.linspace(-2, 4, 100)  
  y = funcion(f,x)

  plt.plot(x, y, label=f)
  plt.axhline(0, color='black',linewidth=0.5)  
  plt.axvline(0, color='black',linewidth=0.5) 
  
  yPunto = funcion(f,x)
  plt.scatter(xPuntos, yPunto, color='red', zorder=5)

  plt.xlim(-2, 4)
  plt.ylim(-10, 10)
  plt.title(f)
  plt.xlabel('x')
  plt.ylabel('f(x)')
  plt.legend()
  plt.grid(True)
  plt.show()

f="(-0.5*x**2)+2.5*x+4.5"
raiz= falsaPosicion(f,6,7,1)
print(raiz[1])