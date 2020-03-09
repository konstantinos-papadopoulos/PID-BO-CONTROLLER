
# Copyright (c) 2016 Konstantinos G. Papdopoulos. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License v1.0 which accompanies this distribution,
# and is available at http://www.eclipse.org/legal/epl-v10.html

"""
Controlled process
"""

from collections import defaultdict
from scipy.interpolate import pade
import control
import logging
import matplotlib
import matplotlib.pyplot as plt
import random

def create_plant(args):
    """Creates the controlled process determined by the user input
       Plant creation supports poles (max:5), zeros (max:5), time delay

    :param poles
    :param zeros
    :param time_delay
    :param user_defined_plant
    :
    :type
    :type
    :type
    :type
    """
    kp = random.random();
    polesOrder = args.poles
    zerosOrder = args.zeros
    time_delay = args.time_delay
    poles = []; polesSorted = [];
    zeros = []; zerosSorted = [];
    # poles generation
    # --------------------------------------------------------------------------
    for i in range(int(polesOrder)):
        if polesOrder == 0:
            poles = [0]
        else:
            poles.append(random.random())
    # zeros generation
    # --------------------------------------------------------------------------
    for j in range(int(zerosOrder)):
        if zerosOrder == 0:
            zeros = [0]
        else:
            zeros.append(random.random())
    # delay generation
    # --------------------------------------------------------------------------
    timeDelay =  random.random()
    
    # Sort poles and zeros to identify dominant pole and dominant zeros
    # --------------------------------------------------------------------------
    polesSorted = sorted(poles, reverse=True)
    zerosSorted = sorted(zeros, reverse=True)

    # plant pole/zeros coefficients using list comprehensions to create tzi,tpj
    # --------------------------------------------------------------------------
    tzi = {}; tpj = {}; tzi_names = []; tpj_names = [];
    if int(zerosOrder) > 0:
        tzi_names = ['tz'+str(i+1) for i in range(int(zerosOrder))]
        for i in range(len(tzi_names)):
            tzi[tzi_names[i]] = zerosSorted[i]
            
    if int(polesOrder) > 0:    
        tpj_names = ['tp'+str(j+1) for j in range(int(polesOrder))]  
        for j in range(len(tpj_names)):
            tpj[tpj_names[j]] = polesSorted[j]        
    print(polesOrder)
        
    q0 = 1 ; q1 = sum(tzi.values());
      
    print('tzi-dictionary:',tzi)
    print('-------------------------------------------------------------------')
    q0 = 1
    print('-------------------------------------------------------------------')
    print('-------------------------------------------------------------------')
    q1 = 0 
    q1 = sum(tzi.values());
    print('-------------------------------------------------------------------')
    q2 = 0 
    for i in range(int(zerosOrder)):
        for j in range(int(zerosOrder)):
            if  (i != j) and (i < j): 
                print('i:',i+1,'j:',j+1); 
                q2 = tzi['tz'+str(i+1)]*tzi['tz'+str(j+1)] + q2
                print('q2->',q2)
    print('-------------------------------------------------------------------')
    q3 = 0 
    for i in range(int(zerosOrder)):
        for j in range(int(zerosOrder)):
            for k in range(int(zerosOrder)):
               if  (i != j != k != i) and (i < j < k): 
                   print('i:',i+1,'j:',j+1,'k:',k+1);
                   q3 = tzi['tz'+str(i+1)]*tzi['tz'+str(j+1)]*tzi['tz'+str(k+1)] + q3
                   print('q3->',q3)
    print('-------------------------------------------------------------------')
    q4 = 0
    for i in range(int(zerosOrder)):
        for j in range(int(zerosOrder)):
            for k in range(int(zerosOrder)):
                for l in range(int(zerosOrder)):
                    if  (i != j != k != l != i) and (i < j < k < l): 
                        print('i:',i+1,'j:',j+1,'k:',k+1,'l:',l+1);
                        q4 = tzi['tz'+str(i+1)]*tzi['tz'+str(j+1)]*tzi['tz'+str(k+1)]*tzi['tz'+str(l+1)] + q4
                        print('q4->',q4)
                        
    print('tpj-dictionary:',tpj)
    print('-------------------------------------------------------------------')
    pGp0 = 1
    print('-------------------------------------------------------------------')
    pGp1 = 0 
    pGp1 = sum(tpj.values());
    print('-------------------------------------------------------------------')
    pGp2 = 0 
    for i in range(int(polesOrder)):
        for j in range(int(polesOrder)):
            if  (i != j) and (i < j): 
                print('i:',i+1,'j:',j+1); 
                pGp2 = tpj['tp'+str(i+1)]*tpj['tp'+str(j+1)] + pGp2
                print('pGp2->',pGp2)
    print('-------------------------------------------------------------------')
    pGp3 = 0
    for i in range(int(polesOrder)):
        for j in range(int(polesOrder)):
            for k in range(int(polesOrder)):
               if  (i != j != k != i) and (i < j < k): 
                   print('i:',i+1,'j:',j+1,'k:',k+1);
                   pGp3 = tpj['tp'+str(i+1)]*tpj['tp'+str(j+1)]*tpj['tp'+str(k+1)] + pGp3
                   print('pGp3->',pGp3)

    print('-------------------------------------------------------------------')
    pGp4 = 0
    for i in range(int(polesOrder)):
        for j in range(int(polesOrder)):
            for k in range(int(polesOrder)):
                for l in range(int(polesOrder)):
                    if  (i != j != k != l != i) and (i < j < k < l): 
                        print('i:',i+1,'j:',j+1,'k:',k+1,'l:',l+1);
                        pGp4 = tpj['tp'+str(i+1)]*tpj['tp'+str(j+1)]*tpj['tp'+str(k+1)]*tpj['tp'+str(l+1)] + pGp4
                        print('pGp4->',pGp4)
    print('-------------------------------------------------------------------')
    pGp5 = 0
    for i in range(int(polesOrder)):
        for j in range(int(polesOrder)):
            for k in range(int(polesOrder)):
                for l in range(int(polesOrder)):
                    for m in range(int(polesOrder)):
                        if  (i != j != k != l != m != i) and (i < j < k < l < m): 
                            print('i:',i+1,'j:',j+1,'k:',k+1,'l:',l+1,'m:',m+1);
                            pGp5 = tpj['tp'+str(i+1)]*tpj['tp'+str(j+1)]*tpj['tp'+str(k+1)]*tpj['tp'+str(l+1)]*tpj['tp'+str(m+1)] + pGp5
                            print('pGp5->',pGp5)
    print('-------------------------------------------------------------------')                            
    print('pGp0->',pGp0)
    print('pGp1->',pGp1)
    print('pGp2->',pGp2)
    print('pGp3->',pGp3)
    print('pGp4->',pGp4)
    print('pGp5->',pGp5)
    
    # Pade approximation
    # --------------------------------------------------------------------------
    n1 = 0
    n2 = 0
    n3 = pow(timeDelay,3) / 6
    
    if timeDelay != 0:
        n1 = timeDelay
        n2 = pow(timeDelay,2) / 2.0
        n3 = pow(timeDelay,3) / 6.0
    e_exp = [1.0, 1.0, (1.0/2.0) * pow(timeDelay,2), \
             (1.0/6.0) * pow(timeDelay,3), (1.0/24.0) * pow(timeDelay,4), \
             (1.0/120.0) * pow(timeDelay,5)]                     
    p, q = pade(e_exp, 5)
    print(p)
    print(q)
    
    # plant poles coefficients dynamics
    # TODO: create a dynamic dictionary to tp_i based on the order of zeros 
    # coming from the user    
    # --------------------------------------------------------------------------
 #==============================================================================
 #    pGp0 = 1; pGp1 = sum(tpj.values());
 #     
 #    pGp2 = tp1*tp2 + tp1*tp3 + tp1*tp4 + tp1*tp5 + tp2*tp3 + \
 #           tp2*tp4 + tp2*tp5 + tp3*tp4 + tp3*tp5 + tp4*tp5;
 # 
 #    pGp3 = tp1*tp2*tp3 + tp1*tp2*tp4 + tp1*tp2*tp5 + tp1*tp3*tp4 + \
 #           tp1*tp3*tp5 + tp1*tp4*tp5 + tp2*tp3*tp4 + tp2*tp3*tp5 + tp2*tp4*tp5 \
 #           + tp3*tp4*tp5;
 # 
 #    pGp4 = tp1*tp2*tp3*tp4 + tp1*tp2*tp3*tp5 + tp1*tp2*tp4*tp5 + \
 #           tp1*tp3*tp4*tp5 + tp2*tp3*tp4*tp5;
 #    pGp5 = tp1*tp2*tp3*tp4*tp5
 # 
    numGp = [q4, q3, q2, q1, q0];
    kpnumGp = [i * kp for i in numGp]
    denGp  = [pGp5, pGp4, pGp3, pGp2, pGp1, pGp0];
    plantGp = control.tf(kpnumGp,denGp);
 #==============================================================================
    
    (T , yout) = control.step_response(plantGp)
    plotStepResponse(T, yout, kp)
    
    return plantGp

def plantRandom(args):
    alpha = 0.1

    pass

def plantDetermined(args):

    pass

def plantAlpha(args):
    
    pass

def plotStepResponse(T, yout, kp):
    """ The methdo takes the vector of a step response in T, yout form plus the
    dc gain kp at steady state and plot on Figure with the step response plus a
    straight horizontal line based in kp  
    """
    kpLine = [1] * len(T)
    fig, ax = plt.subplots()
    ax.plot(T, yout)
    kpLine = [i * kp for i in kpLine]
    ax.plot(T, kpLine)
    ax.set(xlabel='time (s)', ylabel='$y_{out}$',
           title='Plant $G_p$ Open Loop step response')
    ax.legend([ 'step response','dc-gain: kp:' + str(format(kp, '.2f'))])
    plt.grid(color='k', linestyle='-.', linewidth=0.1)
    plt.show()


def plant_estimation(poles,zeros,time_delay,user_defined_plant):
    """The method takes the real plant as an input and returns the settling
       time of the plant's step response, the estimated DC gain (plant's steady
       state gain, an estimation of the plant's unmodelled dynamics and the
       plant transfer function [first order mode] based on the step response
       measurements). All above measurements will be used as an input for
       initializing the automatic tuning algorith,

       Plant creation supports poles (max:5), zeros (max:5), time delay

    :param poles
    :param zeros
    :param time_delay
    :param user_defined_plant
    :
    :type
    :type
    :type
    :type
    """

    pass

def plant_normalize():
    """The method takes the real plant as an input and returns the settling
       time of the plant's step response, the estimated DC gain (plant's steady
       state gain, an estimation of the plant's unmodelled dynamics and the
       plant transfer function [first order mode] based on the step response
       measurements). All above measurements will be used as an input for
       initializing the automatic tuning algorith,

       Plant creation supports poles (max:5), zeros (max:5), time delay
       Plant creation supports poles (max:5), zeros (max:5), time delay

    :param poles
    :param zeros
    :param time_delay
    :param user_defined_plant
    :
    :type
    :type
    :type
    :type
    """       
       
       
    pass
