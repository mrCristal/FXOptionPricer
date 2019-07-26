import scipy.stats as si
from math import log,exp
import datetime
import numpy as np
import matplotlib as mtl
import matplotlib.pyplot as plt
from pandas import Series
from PricerClass import MCpricePayment
from dateutil import parser
from PricerClass import BarrierAnalyticalPrice,BSpriceS

class ClassError(Exception):
    pass

class BarrierContingentVanilla():
    '''Pays out a vanilla payoff conditional on a barrier being hit (or not)'''
    def __init__(self,spot,strike,expiry,drate,frate,vol,barrier,Type,flvr='call',cpair='Not set',notional=100000):
        self.spot = spot
        self.strike = strike
        self.expiry = expiry
        self.drate = drate
        self.frate = frate
        self.vol = vol
        self.flvr = self.__getFlavour__(flvr)
        self.cpair = cpair
        self.notional = notional
        self.barrier = barrier
        self.tau = self.to_maturity()
        self.mu = self.__get_mu__()
        self.lbda = (self.mu**2 + 2*self.drate/self.vol**2)**0.5
        self.nu = self.__get_nu__(S=spot,H=barrier)
        self.Type = Type
    
    def __get_mu__(self,drate=None,frate=None,vol=None):
        if drate is None: drate = self.drate
        if frate is None: frate = self.frate
        if vol is None: vol = self.vol
        return ((drate-frate) - 0.5*vol**2)/vol**2
    
    def __get_nu__(self,S,H):
        if S>H:
            return 1
        return -1

    def __getFlavour__(self,flavour):
        if type(flavour) is not str:
            raise ClassError('Please specify call or put')
        if 'c' in flavour.lower():
            return 1
        else:
            return -1
    def __formatDate__(self,date):
        if type(date) is datetime.date:
            return date
        elif type(date) is str:
            return parser.parse(date).date()
        else:
            raise ClassError('Invalid expiry format')
    
    def to_maturity(self, date=None):
        if date is None: 
            return (self.__formatDate__(self.expiry) - datetime.date.today()).days/365
        return (self.__formatDate__(date)-datetime.date.today()).days/365
    
    def format_price(self,price, premium='domestic'):
        price = self.notional*price
        if premium.lower() == 'domestic':
            return round(price,2)
        elif premium.lower() == 'foreign':
            return round(price/self.spot,2)
        elif premium.lower() == 'domestic/foreign':
            return round(price/self.notional,6)
        elif premium.lower() == 'foreign/domestic':
            return round(price/(self.notional*self.strike*self.spot),6)
        elif premium.lower() == '%domestic':
            return round(price/(self.notional*self.strike),6)
        elif premium.lower() == '%foreign':
            return round(price/(self.notional*self.spot),6)

    def price(self,flvr=None,S=None,fi=None,di=None,tau=None,K=None,vol=None,mu=None,H=None,nu=None,barrierType=None,premium='domestic'):
        '''Returns the analytical price'''
        if flvr is None: flvr = self.flvr
        if   S  is None: S   =  self.spot
        if   H  is None: H  =   self.barrier
        if  fi  is None: fi  =  self.frate
        if  di  is None: di  =  self.drate
        if  tau is None: tau =  self.tau
        if   K  is None: K  =   self.strike
        if  vol is None: vol =  self.vol
        if  mu  is None: mu  =  self.mu
        if  nu  is None: nu  =  self.nu
        if barrierType is None: barrierType = self.Type

        domPrice = BarrierAnalyticalPrice(spot=S,strike=K,barrier=H,_type=barrierType,drate=di,frate=fi,vol=vol,nu=nu,flvr=flvr,mu=mu,tau=tau).price()
        return self.format_price(price=domPrice,premium=premium)

    def plot_dvds(self, premium='domestic',with_corrVanilla=False):
        dvds = [[],[]]
        for x in np.arange(0.75*self.barrier,1.25*self.barrier,0.01):
            dvds[0].append(round(self.price(premium=premium,S=x),6))
            dvds[1].append(round(x,6))
        plt.plot(Series(dvds[0],dvds[1]))
        if with_corrVanilla:
            dvds2 = [[],[]]
            for x in np.arange(0.9*self.spot,1.1*self.spot,0.01):
                dvds2[0].append(round(self.format_price(BSpriceS(spot=x,strike=self.strike,dt=self.tau,dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).price(),premium=premium),6))
                dvds2[1].append(round(x,6))
                Series(dvds2[0],dvds2[1]).plot(linestyle='--',secondary_y=True)
        plt.ylabel('Option Value in {}'.format(premium.capitalize()))
        plt.xlabel('Spot Value')
        plt.show()
        plt.clf()
  
class BarrierContingentPayment():
    '''Pays out a fixed amount depending on whether a barrier(s) has been hit or not'''
    def __init__(self,S,fr,dr,vol,expiry,Type,notional=100000,**kwargs):
        barrierLoc = ['upper','lower','down','up']
        validKeys = [loc+'barrier' for loc in barrierLoc]+[loc+'_barrier' for loc in barrierLoc]+['barrier','h']
        self.S = S
        self.fr = fr
        self.dr = dr
        self.vol = vol
        self.expiry = expiry
        self.notional = notional
        self.Type=Type
        self.dt = self.to_maturity()
        for k in kwargs.keys():
            if k.lower() in validKeys:
                continue
            else:
                raise ClassError('Please specify barrier as either upperbarrier or lowerbarrier or both')
        if len(kwargs) == 2:
            self.upper = kwargs[max(kwargs, key=kwargs.get)]
            self.lower = kwargs[min(kwargs, key=kwargs.get)]
        elif len(kwargs)==1:
            self.barrier = list(kwargs.values())[0]
        else:
            raise ClassError('Please specify a max of 2 barriers')

    def __formatDate__(self,date):
        if type(date) is datetime.date:
            return date
        elif type(date) is str:
            return parser.parse(date).date()
        else:
            raise ClassError('Invalid expiry format')
    
    def to_maturity(self, date=None):
        if date is None: 
            return (self.__formatDate__(self.expiry) - datetime.date.today()).days/365
        return (self.__formatDate__(date)-datetime.date.today()).days/365
    
    def price(self,sims=10000,premium='domestic',steps=None):
        '''Returns the price generated by the Monte-Carlo numerical method'''
        if steps is None:
            steps = int(self.dt*365)
        if hasattr(self,'upper'):
            domPrice = self.notional * MCpricePayment(S=self.S,upper=self.upper,lower=self.lower,dt=self.to_maturity(),dr=self.dr,fr=self.fr,vol=self.vol,Type=self.Type,sims=sims).price(steps=steps)
        else:
            domPrice = self.notional * MCpricePayment(S=self.S,barrier=self.barrier,dt=self.to_maturity(),dr=self.dr,fr=self.fr,vol=self.vol,Type=self.Type,sims=sims).price(steps=steps)

        return self.price_formatter(DomPrice=domPrice,premium=premium)
    
    def priceWithBB(self,sims=10000,premium='domestic',steps=None):
        '''Computes the Monte Carlo price of the option with the use of Brownian Bridges'''
        if steps is None:
            steps = int(self.dt*365)
        if hasattr(self,'upper'):
            domPrice = self.notional * MCpricePayment(S=self.S,upper=self.upper,lower=self.lower,dt=self.to_maturity(),dr=self.dr,fr=self.fr,vol=self.vol,Type=self.Type,sims=sims).priceWithBB(steps=steps)
        else:
            domPrice = self.notional * MCpricePayment(S=self.S,barrier=self.barrier,dt=self.to_maturity(),dr=self.dr,fr=self.fr,vol=self.vol,Type=self.Type,sims=sims).priceWithBB(steps=steps)
        return self.price_formatter(DomPrice=domPrice,premium=premium)

    def price_formatter(self,premium, DomPrice):
        if premium.lower() == 'domestic':
            return round(DomPrice,2)
        elif premium.lower() == 'foreign':
            return round(DomPrice/self.S,2)