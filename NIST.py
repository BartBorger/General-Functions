# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:15:50 2024

@author: jaborger
source for Al Cu Fe W:https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir84-3007.pdf 
source for others data range 4-300K: https://trc.nist.gov/cryogenics/materials/materialproperties.htm 
source for phosfor bronze based on lakeshore wire https://ntrs.nasa.gov/api/citations/20090032058/downloads/20090032058.pdf
Feel free to add interpolations for other metals not yet included in this function
"""


def NIST_Lambda(T,RRR,Material): 
    """function returning thermal conductivity of Material as a function of temperature and RRR (if applicable). Curves taken from the  NIST database
    
        :T:     Temperature [K]
        :RRR:   Residual Resistivity Ratio [-]
        :Material: String to identify material can be[Al; Cu; Fe; W; 304; 316; G-10N; G-10W; Kapton; Brass;Mylar;P-Bronze]
    """
    import numpy as np
    import warnings
    

    # source:https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir84-3007.pdf
    L_0 = 2.443e-8            # lorentz ratio  
    match Material.lower():
        case "al":
            resistivity_273 = 24.8e-9 
            resistivity_0 = (resistivity_273)/(RRR-1)
            
            #constants
            beta = resistivity_0/L_0
            W_0 = beta/T
            
            P1 = 4.716e-8
            P2 = 2.446
            P3 = 623.6
            P4 = -0.16
            P5 = 130.9
            P6 = 2.5
            P7 = 0.8168
            #resistivity_c is a term for mathmatical residuals
            W_c=-0.0005*np.log(T/330)*np.exp(-(np.log(T/380)/0.6)**2)-0.0013*np.log(T/110)*np.exp(-(np.log(T/94)/0.5)**2)
            
            W_i = P1*T**P2/(1+P1*P3*T**(P2+P4)*np.exp(-(P5/T)**P6))+W_c
            W_io = P7*W_i*W_0/(W_0+W_i)
            lambda_Al = 1/(W_0+W_i+W_io)
            return lambda_Al
        case "cu":
            resistivity_273 = 15.4e-9 # 15.4 nOhm/m 
            resistivity_0 = (resistivity_273)/(RRR-1)
            
            #constants
            beta = resistivity_0/L_0
            W_0 = beta/T
            beta_r = beta/0.0003
            P1 = 1.754e-8
            P2 = 2.763
            P3 = 1102
            P4 = -0.165
            P5 = 70
            P6 = 1.765
            P7 = 0.838/beta_r**0.1661
            #resistivity_c is a term for mathmatical residuals
            W_c=-0.00012*np.log(T/420)*np.exp(-(np.log(T/470)/0.7)**2)-0.00016*np.log(T/73)*np.exp(-(np.log(T/87)/0.45)**2)-0.00002*np.log(T/18)*np.exp(-(np.log(T/21)/0.5)**2)
            
            W_i = P1*T**P2/(1+P1*P3*T**(P2+P4)*np.exp(-(P5/T)**P6))+W_c
            W_io = P7*W_i*W_0/(W_0+W_i)
            lambda_Cu = 1/(W_0+W_i+W_io)
            return lambda_Cu
        case "fe":
            resistivity_273 = 87e-9 
            resistivity_0 = (resistivity_273)/(RRR-1)
            
            #constants
            beta = resistivity_0/L_0
            W_0 = beta/T
            P1 = 166.9e-8
            P2 = 1.868
            P3 = 1.503e5
            P4 = -1.22
            P5 = 238.6
            P6 = 1.392
            P7 = 0
            #resistivity_c is a term for mathmatical residuals
            W_c=-0.004*np.log(T/440)*np.exp(-(np.log(T/650)/0.8)**2)-0.002*np.log(T/90)*np.exp(-(np.log(T/90)/0.45)**2)
            
            W_i = P1*T**P2/(1+P1*P3*T**(P2+P4)*np.exp(-(P5/T)**P6))+W_c
            W_io = P7*W_i*W_0/(W_0+W_i)
            lambda_Fe = 1/(W_0+W_i+W_io)
            return lambda_Fe
        case "w":
            resistivity_273 = 48.4e-9 
            resistivity_0 = (resistivity_273)/(RRR-1)
            
            #constants
            beta = resistivity_0/L_0
            W_0 = beta/T
            P1 = 31.70e-8
            P2 = 2.29
            P3 = 541.3
            P4 = -0.22
            P5 = 69.94
            P6 = 3.557
            P7 = 0.0
            #resistivity_c is a term for mathmatical residuals
            W_c=-0.00085*np.log(T/130)*np.exp(-(np.log(T/230)/0.7)**2)+0.00015*np.exp(-(np.log(T/3500)/0.8)**2)+0.0006*np.log(T/90)*np.exp(-(np.log(T/80)/0.4)**2)+0.0003*np.log(T/24)*np.exp(-(np.log(T/33)/0.5)**2)
            
            W_i = P1*T**P2/(1+P1*P3*T**(P2+P4)*np.exp(-(P5/T)**P6))+W_c
            W_io = P7*W_i*W_0/(W_0+W_i)
            lambda_W = 1/(W_0+W_i+W_io)
            return lambda_W
        case "304":
            #note 304L has the same coefficients
            if any(T>300) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<300")
            a =-1.4087
            b=1.3982
            c=0.2543
            d=-0.6260
            e=0.2334
            f=0.4256
            g=-0.4658 	
            h=0.1650 	
            i=-0.0199
            lambda_304 =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_304
        case "316":
            if any(T>300) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<300")
            a =-1.4087
            b=1.3982 	
            c=0.2543
            d=-0.6260
            e=0.2334
            f=0.4256
            g=-0.4658 	
            h=0.1650 	
            i=-0.0199
            lambda_316 =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_316
        case "g-10n":
            if any(T>300) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<300")
            a =-4.1236
            b=13.788 	
            c=-26.068
            d=26.272
            e=-14.663
            f=4.4954
            g=-0.6905 	
            h=0.0397
            i=0
            lambda_G10N =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_G10N
        case "g-10w":
            if any(T>300) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<300")
            a =-2.64827
            b=8.80228 	
            c=-24.8998
            d=41.1625 	
            e=-39.8754 	
            f=23.1778 	
            g=-7.95635
            h=1.48806 	
            i=-0.11701 	
            lambda_G10W =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_G10W
        case "kapton":
            if any(T>300) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<300")
            a =5.73101 	
            b=-39.5199
            c=79.9313
            d=-83.8572 	
            e=50.9157 	
            f=-17.9835 	
            g=3.42413
            h=-0.27133
            i=0 	
            lambda_Kapton =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_Kapton
        case "brass":
            #UNS C2600
            if any(T>116) or any(T<5):
                warnings.warn("Temperature outside range of validity: 5<T<116")
            a =0.021035
            b=-1.01835
            c=4.54083
            d=-5.03374
            e=3.20536
            f=-1.12933
            g=0.174057
            h=-0.27133
            i=-0.0038151	
            lambda_Brass =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_Brass
        case "mylar":
            if any(T>83) or any(T<1):
                warnings.warn("Temperature outside range of validity: 1<T<83")
            a=-1.37737
            b=-3.40668
            c=20.5842
            d=-53.1244
            e=73.2476
            f=-57.6546
            g=26.1192
            h=-6.34790
            i=0.640331
            lambda_mylar =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
            return lambda_mylar
        case "p-bronze":
            warnings.warn("No NIST data available, fit from wires in JWST")
            if any(T>295) or any(T<4):
                warnings.warn("Temperature outside range of validity: 4<T<295")
            a=-10.9482 
            b=28.4752
            c=-32.3378
            d=20.9036 
            e=-8.05399 
            f=1.90329 
            g=-0.271774
            h=0.0215998 
            i=-7.35095e-4
            lambda_p_bronze =np.exp(a+b*(np.log(T)) + c*(np.log(T))**2 + d*(np.log(T))**3 + e*(np.log(T))**4 + f*(np.log(T))**5 + g*(np.log(T))**6 + h*(np.log(T))**7 + i*(np.log(T))**8 )
            return lambda_p_bronze
        
        case _:
            raise ValueError("Unknown material please check spelling or add it to this function.")

def NIST_Specific_Heat(T,RRR,Material): 
    """function returning specific heat of Material as a function of temperature. Curves taken from the  NIST database
        Data to be added 
    
        :T:     Temperature [K]
        :Material: String to identify material can be[Al; Cu; Fe; W; 304; 316; G-10N; G-10W; Kapton; Brass;Mylar;P-Bronze]
    """
    import numpy as np
    import warnings
    
    match Material.lower():
        case "304":
                #note 304L has the same coefficients
                if any(T>300) or any(T<4):
                    warnings.warn("Temperature outside range of validity: 4<T<300")
                a= 22.0061
                b= -127.5528
                c= 303.647
                d= -381.0098
                e= 274.0328
                f= -112.9212
                g= 24.7593
                h= -2.239153
                i= 0 
                Cv_304 =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
                return Cv_304
            case "316":
                if any(T>300) or any(T<4):
                    warnings.warn("Temperature outside range of validity: 4<T<300")
                if T>50: 
                    a=-1879.464
                    b= 3643.198
                    c= 76.70125
                    d= -6176.028
                    e= 7437.6247
                    f= -4305.7217
                    g= 1382.4627
                    h= -237.22704
                    i= 17.05262
                else:
                    a= 12.2486
                    b= -80.6422
                    c= 218.743
                    d= -308.854
                    e= 239.5296
                    f= -89.9982
                    g= 3.15315 
                    h= 8.44996
                    i= -1.91368 
                Cv_316 =10**(a+b*(np.log10(T)) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8 )
                return Cv_316
