import numpy as np
from scipy.integrate import simps, quad, dblquad
import matplotlib.pyplot as plt


wavelength = 1e-6
aperture_width = 1e-5
k = 2 * np.pi / wavelength

# All of the functions are defined below

def integral_kernel(x_co, x_co_ap, k, z):     # A function that defines the inside of the frensel integral which represents the complex exponetial part
     return np.exp(1j * k / (2 * z) * (x_co - x_co_ap) ** 2)

def fresnel_integral(x_co, k, z, x_co_ap_vals): # A function that integrates the kernel 
     return simps(integral_kernel(x_co, x_co_ap_vals, k, z), x_co_ap_vals, dx = delta_x)

def real_kernel(x_co, x_co_ap, k, z): # A functino that represents the real part of the complex exponential of the frensel integral in 1 dimension
    return np.cos((k/(2*z)) * ((x_co - x_co_ap)**2))

def imag_kernel(x_co, x_co_ap, k, z): # A functino that represents the imaginary part of the complex exponential of the fresnel integral in 1 dimension
    return np.sin((k/(2*z)) * ((x_co - x_co_ap)**2))

def Fresnel2dreal(yp, xp, y, x, k, z): # A function that represents the kernel of the real part of the 2-D frensel equation
    return np.cos(k*((x-xp)**2 + (y-yp)**2)/(2*z))

def Fresnel2dimag(yp, xp, y, x, k, z) : # A function that represents the kernel of the imaginary part of the 2-D frensel equation
    return np.sin(k*((x-xp)**2 + (y-yp)**2)/(2*z))

def lower_limit(yp): # A function that represents the lower limit for the circular aperture and the 1e-5 is the radius of the circle
    return -np.sqrt(1e-5**2 - (yp)**2)

def upper_limit(yp): # A function that represents the upper limit for the circular aperture
    return np.sqrt(1e-5**2 - (yp)**2)


MyInput = 0 # This is a menu
while MyInput != 'q':
    print('Choice "a" is for 1D diffraction using simpson rule for integration')
    print('Choice "b" is for 1D diffraction using quadrature integration')
    print('Choice "c" is for 2D diffraciton using a square aperture')
    print('Choice "d" is for 2D diffraction using a circular aperture')
    print('') 
    MyInput = input ( 'Enter a choice, "a", "b", "c", "d" or "q" to quit: ')
    print('')
    print('You entered the choice: ', MyInput)
    print('') 
    if MyInput == 'a' :
        print ('You have chosen part (a)' )
        print('')
        N = 101 # Number used for the 
        delta_x = (aperture_width)/(N-1) 
        x_co_ap_vals = np.linspace(-aperture_width / 2, aperture_width / 2, N) # The creation of an array for every position of the aperture
        x_co_vals = np.linspace(-0.005, 0.005, N) # The creation of an array for every position of the screen 
        fresnel_integral_vals = np.zeros(x_co_vals.shape, dtype=complex) # An array for each integral 
        z = float(input('Enter a value for z (screen distance) in meters: '))

        for i, x_co in enumerate(x_co_vals): # It iterates through each value of the screen position  
            fresnel_integral_vals[i] = (fresnel_integral(x_co, k, z, x_co_ap_vals))/z # The integral function above is used and the value is stored in an array


        mod_e = np.abs(fresnel_integral_vals) # The modulus of the integral is found 
        intensity = np.square(mod_e) # The square of the modulus is found
        plt.plot(x_co_vals, intensity) # The graph is plotted
        plt.title('Relative intensity against screen position for 1D diffraction', y = 1.1) # The title was shifted up as it was interfering with the plot
        plt.xlabel('Screen position / m')
        plt.ylabel('Relative intensity / (W/m^2)')
        plt.show()
    elif MyInput == 'b':
        print ('You have chosen part (b)') 
        print('')
        N = 101
        x_co_ap_vals = np.linspace(-aperture_width / 2, aperture_width / 2, N)
        x_co_vals = np.linspace(-0.005, 0.005, N)
        fresnel_integral_vals = np.zeros(len(x_co_vals), dtype=complex)
        z = float(input('Enter a value for z (screen distance) in meters: '))

        for i, x_co in enumerate(x_co_vals): # iterating through each value of the screen array and integrating both real and imaginary part
            real_val, real_error = quad(real_kernel, x_co_ap_vals[0], x_co_ap_vals[-1], args=(x_co, k, z)) # -1 is the last value in the array 
            imag_val, imag_error = quad(imag_kernel, x_co_ap_vals[0], x_co_ap_vals[-1], args=(x_co, k, z))
            fresnel_integral_vals[i] = (real_val + 1j * imag_val)/z # adding the imaginary values and the real parts together, dividing by z 

        mod_e = abs(fresnel_integral_vals)
        intensity = np.square(mod_e)

        plt.plot(x_co_vals, intensity)
        plt.title('Relative intensity against screen position for 1D diffraction', y = 1.1) # The title was shifted up as it was interfering with the plot
        plt.xlabel('Screen position / m')
        plt.ylabel('Relative intensity / (W/m^2)')
        plt.show()
        
    elif MyInput == 'c':
        print ('You have chosen part (c)' )
        print ('')
        
        N = 101
        x = np.linspace(-5e-3, 5e-3, N) # Array for the x coordinates for the screen
        y = np.linspace(-5e-3, 5e-3, N) # Array for the y coordinates for the screen
        xp = np.linspace(-1e-5, 1e-5, N) # Array for the x coordinates of the aperture
        yp = np.linspace(-1e-5, 1e-5, N) # Array for the y coordinates of the aperture
        
        z = float(input('Enter a value for z (screen distance) in meters: '))

        mod = np.zeros((x.size, y.size))

        # Iterating through both x and y values of the screen for each integral
        for i, x_val in enumerate(x): 
            for j, y_val in enumerate(y):
                realpart, real_err = dblquad(Fresnel2dreal ,yp[0],yp[-1],xp[0],xp[-1],args=(y_val, x_val, k, z)) # integral calculated for the real part and stored into an variable for specific screen coordinate
                imagpart, imag_err = dblquad(Fresnel2dimag,yp[0],yp[-1],xp[0],xp[-1],args=(y_val, x_val, k, z)) # integral calculated for the imaginary part and stored into an variable for specific screen coordinate
                mod[i, j] = (np.abs(realpart + 1j * imagpart))/z

        acc_intensity = np.square(mod)

        plt.imshow(acc_intensity, cmap='plasma', extent=[x[0], x[-1], y[0], y[-1]])
        plt.title('Relative intensity for 2D diffraction with a square aperture', y = 1.1) # The title was shifted up as it was interfering with the plot
        plt.xlabel('x screen coordinate / m')
        plt.ylabel('y screen coordinate / m')
        cbar = plt.colorbar()
        cbar.set_label('Relative intensity / (W/m^2)')
        cbar.ax.set_position([0.85, 0.1, 0.05, 0.80]) # Adjusting the position of the colour bar so it doesnt interfere with the plot
        plt.show()
        
    elif MyInput == 'd':
        print('You have chosen part (d)')
        print('')
        
        N = 101
        screen_width = 1e-2
        x = np.linspace(-screen_width/2, screen_width/2, N)
        y = np.linspace(-screen_width/2, screen_width/2, N)
        xp = np.linspace(-1e-5, 1e-5, N)
        yp = np.linspace(-1e-5, 1e-5, N)
        k = 2 * np.pi / 1e-6
        z = float(input('Enter a value for z (screen distance) in meters: '))
        x0 = 0
        y0 = 0
        R = 1e-5


        intensity = np.zeros((x.size, y.size))
        
        for i, x_val in enumerate(x):
            for j, y_val in enumerate(y):
                realpart, real_err = dblquad(Fresnel2dreal , yp[0], yp[-1], lower_limit, upper_limit, args=(y_val, x_val, k, z)) # For real part of the integral the x limits are replaced by the function describing the lower and upper limit for a circle
                imagpart, imag_err = dblquad(Fresnel2dimag, yp[0], yp[-1], lower_limit, upper_limit, args=(y_val, x_val, k, z)) # For imaginary part of the integral the x limits are replaced by the function describing the lower and upperlimit for a circle
                intensity[i, j] = (np.abs(realpart + 1j * imagpart))/z
        acc_intensity = np.square(intensity)

        plt.imshow(acc_intensity, cmap='plasma', extent=[x[0], x[-1], y[0], y[-1]])
        plt.title('Relative intensity for 2D diffraction with a circular aperture', y = 1.1) # The title was shifted up as it was interfering with the plot
        plt.xlabel('x screen coordinate / m')
        plt.ylabel('y screen coordinate / m')
        cbar = plt.colorbar()
        cbar.set_label('Relative intensity / (W/m^2)')
        cbar.ax.set_position([0.85, 0.1, 0.05, 0.80])
        plt.show()
    elif MyInput != 'q':
        print('This is not a valid choice')
        print('You have chosen to finish - goodbye.')
