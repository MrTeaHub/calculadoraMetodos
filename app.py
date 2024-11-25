from flask import Flask, render_template, request, Response, redirect, url_for
from sympy import sympify, symbols, lambdify, diff, sin, cos, tan, exp, sqrt, integrate
import numpy as np
import matplotlib.pyplot as plt
import io
import base64

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

# Funciones de los métodos
def Ea(xNuevo,xAnterior):
    return abs(((xNuevo-xAnterior)/xNuevo))*100

def integrarFuncionDefinida(funcionTexto, a, b):
    x = symbols('x')
    func = sympify(funcionTexto)
    return integrate(func, (x, a, b))
    
def funcion(funcionTexto, x, derivada=False):
    simbolo = symbols('x')  
    expresion = sympify(funcionTexto, locals={"sin": sin, "cos": cos, "tan": tan, "exp": exp, "sqrt":sqrt})  
    
    if derivada:
        expresion = diff(expresion, simbolo) 
    func = lambdify(simbolo, expresion, "numpy")
    return func(x)

def newtonRaphson(funcionTexto, xi, errorPorcentual):
    i = 0
    rta = []
    raices = [xi]
    errores = []
    while(True):
        x = xi - (funcion(funcionTexto,xi)/funcion(funcionTexto,xi,derivada=True))
        e = Ea(x,xi)
        raices.append(x)
        errores.append(e)
        if i > 1 and e < (errorPorcentual):
            rta.append(raices)
            rta.append(errores)
            break
        xi = x
        i += 1
    return rta

def falsaPosicion(funcionTexto, xl, xu, errorPorcentual):
    xrAnterior = 0
    rta = []
    errores = []
    Xl = []
    Xu = []
    Xr = []
    while True:
        xr = xu-((funcion(funcionTexto, xu)*(xl-xu))/(funcion(funcionTexto, xl) -funcion(funcionTexto, xu)))
        x = round(funcion(funcionTexto, xl)*funcion(funcionTexto, xr),8)
        e = Ea(xr,xrAnterior)
        errores.append(e)
        Xr.append(xr)
        Xu.append(xu)
        Xl.append(xl)
        if x < 0:
            xu = xr
        elif x > 0:
            xl = xr
            
        if x == 0 or e < errorPorcentual:
            rta.append(Xr)
            rta.append(Xl)
            rta.append(Xu)
            rta.append(errores)
            break
        xrAnterior=xr
    return rta

def reglaTrapecio(funcionTexto, a, b, n, h):
    integralAprox = (funcion(funcionTexto,a) + funcion(funcionTexto,b)) / 2
    
    for i in range(1, n):
        integralAprox += funcion(funcionTexto, a + i * h)
    
    integralAprox *= h
    return integralAprox

def diferenciacionNumerica(funcionTexto, xi, h, orden):
    simbolo = symbols('x')
    expresion = sympify(funcionTexto)
    func = lambdify(simbolo, expresion, "numpy")
    
    primerDerivada = diff(expresion, simbolo)
    segundaDerivada = diff(primerDerivada, simbolo)
    
    funcionPrimerDerivada = lambdify(simbolo, primerDerivada, "numpy")
    funcionSegundaDerivada = lambdify(simbolo, segundaDerivada, "numpy")
    
    xiMenos3 = xi-(3*h)
    xiMenos2 = xi-(2*h)
    xiMenos1 = xi-h
    xiMas1 = xi+h
    xiMas2 = xi+(2*h)
    xiMas3 = xi+(3*h)
    
    rta = []
    derivadas = []
    errores = []
    valorVerdadero = 0
    
    if orden == "primerOrden":
        adelante = (-func(xiMas2)+(4*func(xiMas1))-(3*func(xi)))/(2*h)
        atras = ((3*func(xi))-(4*func(xiMenos1))+func(xiMenos2))/(2*h)
        centrada = (-func(xiMas2)+(8*func(xiMas1))-(8*func(xiMenos1))+func(xiMenos2))/(12*h)
        
        errores.append(Ea( funcionPrimerDerivada(xi),adelante))
        errores.append(Ea( funcionPrimerDerivada(xi),atras))
        errores.append(Ea( funcionPrimerDerivada(xi),centrada))
        
        derivadas.append(adelante)
        derivadas.append(atras)
        derivadas.append(centrada)
        
        valorVerdadero = funcionPrimerDerivada(xi)
    else:
        adelante = (-func(xiMas3)+(4*func(xiMas2))-(5*func(xiMas1))+(2*func(xi)))/(h**2)
        atras = ((2*func(xi))-(5*func(xiMenos1))+(4*func(xiMenos2))-func(xiMenos3))/(h**2)
        centrada = (-func(xiMas2)+(16*func(xiMas1))-(30*func(xi))+(16*func(xiMenos1))-func(xiMenos2))/(12*h**2)
        
        errores.append(Ea( funcionSegundaDerivada(xi),adelante))
        errores.append(Ea( funcionSegundaDerivada(xi),atras))
        errores.append(Ea( funcionSegundaDerivada(xi),centrada))
        
        derivadas.append(adelante)
        derivadas.append(atras)
        derivadas.append(centrada)
        
        valorVerdadero = funcionSegundaDerivada(xi)
    
    rta.append(derivadas)
    rta.append(errores)
    rta.append(valorVerdadero)
    return rta

def grafica(funcionTexto,xi):
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

    # Guardar el gráfico en un objeto en memoria (BytesIO)
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plt.close()

    # Codificar la imagen en base64
    img_base64 = base64.b64encode(img.getvalue()).decode('utf-8')
    return img_base64
    
@app.route('/metodos/newton_raphson', methods=["GET","POST"])
def newton_raphson():
    if request.method == "POST":
        try:
            funcion = request.form["funcion"]
            xi = float(request.form["xi"])
            errorPorcentual = float(request.form["errorPorcentual"])

            respuesta = newtonRaphson(funcion, xi, errorPorcentual)
            raices = respuesta[0]
            errores = respuesta[1]
            iteraciones = list(range(len(raices))) 
            tablaDatos = list(zip(iteraciones,raices,errores))
            img_base64 = grafica(funcion, raices[-2])
            print(raices[-2])
            return render_template('metodos/newton_raphson.html', funcion=funcion, tablaDatos=tablaDatos, raiz=raices[-2],grafica=img_base64, error_mensaje=None)

        except ValueError:
            return render_template('metodos/newton_raphson.html',
            error_mensaje="Error: Asegurarse de ingresar valores numéricos válidos para 'xi' y 'error porcentual'.")

        except Exception as e:
            return render_template('metodos/newton_raphson.html', error_mensaje=f"Error inesperado: {str(e)}")
        
    return render_template('metodos/newton_raphson.html')

@app.route('/metodos/falsa_posicion', methods=["GET","POST"])
def falsa_posicion():
    if request.method == "POST":
        try:
            funcion = request.form["funcion"]
            xl = float(request.form["xl"])
            xu = float(request.form["xu"])
            errorPorcentual = float(request.form["errorPorcentual"])
            respuesta = falsaPosicion(funcion, xl, xu, errorPorcentual)
            Xr = respuesta[0]
            Xl = respuesta[1]
            Xu = respuesta[2]
            errores = respuesta[3]
            
            iteraciones = list(range(len(Xr))) 
            tablaDatos = list(zip(iteraciones,Xl, Xu, Xr, errores))
            img_base64 = grafica(funcion, Xr[-1])
            
            print(Xr)
            print(errores)
            print(len(Xr))
            print(len(errores))
            
            return render_template('metodos/falsa_posicion.html', funcion=funcion, tablaDatos=tablaDatos, raiz=Xr[-1], grafica=img_base64,error_mensaje=None)

        except ValueError:
            return render_template('metodos/falsa_posicion.html',
            error_mensaje="Error: Asegurarse de ingresar valores numéricos válidos para 'xl, xu' y 'error porcentual'.")

        except Exception as e:
            return render_template('metodos/falsa_posicion.html', error_mensaje=f"Error inesperado: {str(e)}")
        
    return render_template('metodos/falsa_posicion.html')

@app.route('/metodos/regla_trapecio', methods=["GET","POST"])
def regla_trapecio():
    if request.method == "POST":
        try:
            funcion = request.form["funcion"]
            a = float(request.form["a"])
            b = float(request.form["b"])
            n = int(request.form["subIntervalos"])
            h = (b - a) / n
            
            respuesta = reglaTrapecio(funcion, a, b, n, h)
            valorVerdadero = integrarFuncionDefinida(funcion, a, b)
            error = Ea(valorVerdadero, respuesta)
            return render_template('metodos/regla_trapecio.html', funcion=funcion,integral=respuesta, a=a, b=b, n=n, h=h, error=error, valorVerdadero=valorVerdadero)

        except ValueError:
            return render_template('metodos/regla_trapecio.html',
            error_mensaje="Error: Asegurarse de ingresar valores numéricos válidos para los limites y el número de subintervalos.")

        except Exception as e:
            return render_template('metodos/regla_trapecio.html', error_mensaje=f"Error inesperado: {str(e)}")
        
    return render_template('metodos/regla_trapecio.html')

@app.route('/metodos/diferenciacion_numerica', methods=["GET","POST"])
def diferenciacion_numerica():
    if request.method == "POST":
        try:
            funcion = request.form["funcion"]
            x = float(request.form["x"])
            h = float(request.form["h"])
            orden = request.form["ordenDerivada"]
            
            respuesta = diferenciacionNumerica(funcion, x, h, orden)
            derivadas = respuesta[0]
            errores = respuesta[1]
            valorVerdadero = respuesta[2]

            tablaDatos = list(zip(derivadas, errores))
            return render_template('metodos/diferenciacion_numerica.html', funcion=funcion, tablaDatos=tablaDatos, valorVerdadero=valorVerdadero, orden=orden)
            
        except ValueError:
            return render_template('metodos/diferenciacion_numerica.html',
            error_mensaje="Error: Asegurarse de ingresar valores numéricos válidos para los limites y el número de subintervalos.")

        except Exception as e:
            return render_template('metodos/diferenciacion_numerica.html', error_mensaje=f"Error inesperado: {str(e)}")
        
    return render_template('metodos/diferenciacion_numerica.html')
