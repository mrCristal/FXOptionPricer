from ContractClass import Contract
import datetime
import scipy.stats as si
from math import log,exp
import numpy as np
import matplotlib.pyplot as plt
from pandas import Series
from PricerClass import BSpriceF,BSpriceS,MCprice
from dateutil import parser

class ClassError(Exception):
    pass

class Vanilla(Contract):
    '''A derivative contract which derives its price from an underlying asset, in this case the spot exchange rate. This contract can be either a call or put'''
    def __init__(self,spot,strike,expiry,drate,frate,vol,flvr='call',forward=None,notional=100000):
        super().__init__(spot,strike,expiry,drate,notional)
        self.frate = frate
        self.forward = forward
        self.vol = vol
        self.flvr = self.__getFlavour__(flvr)

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
        else:
            return (self.__formatDate__(date)-datetime.date.today()).days/365
    
    def __price_formatter__(self,premium,DomPrice):
        if premium.lower() == 'domestic':
            return round(DomPrice,2)
        elif premium.lower() == 'foreign':
            return round(DomPrice/self.spot,2)
        elif premium.lower() == 'domestic/foreign':
            return round(DomPrice/self.notional,6)
        elif premium.lower() == 'foreign/domestic':
            return round(DomPrice/(self.notional*self.strike*self.spot),6)
        elif premium.lower() == '%domestic':
            return round(DomPrice/(self.notional*self.strike),6)
        elif premium.lower() == '%foreign':
            return round(DomPrice/(self.notional*self.spot),6)
    
    def spotPrice(self,spot=None,strike=None,vol=None,expiry=None,premium='domestic'):
        '''Returns the analytical price with a spot hedge'''
        if spot is None: spot = self.spot
        if strike is None: strike = self.strike
        if vol is None: vol = self.vol
        to_maturity = self.to_maturity(date=expiry)
        domPrice = self.notional * BSpriceS(spot=spot,strike=strike,dt=to_maturity,dr=self.drate,fr=self.frate,vol=vol,flvr=self.flvr).price()
        return self.__price_formatter__(DomPrice=domPrice,premium=premium)
        

    def forwardPrice(self,forward=None,strike=None,expiry=None,vol=None, premium='domestic'):
        '''Returns the analytical price with a forward hedge if it is given'''
        if self.forward is not None:
            if forward is None: forward = self.forward
            if strike is None: strike = self.strike
            if vol is None: vol = self.vol
            to_maturity = self.to_maturity(date=expiry)
            domPrice = self.notional * BSpriceF(forward=forward,strike=strike,dt=to_maturity,dr=self.drate,vol=vol,flvr=self.flvr).price()
            return self.__price_formatter__(DomPrice=domPrice,premium=premium)
        raise ClassError('No forward argument was specified')
    
    def MCprice(self,sims=10000,premium='domestic',withLog=False):
        '''Returns the numerically derived price with a spot hedge only. This Monte-Carlo method employs anti-thetic variables and quasi-random numbers'''
        domPrice = self.notional * MCprice(S=self.spot,K=self.strike,dt=self.to_maturity(),fr=self.frate,dr=self.drate,vol=self.vol,flvr=self.flvr,sims=sims).price(withLog=withLog)
        return self.__price_formatter__(DomPrice=domPrice,premium=premium)
    
    def standard(self,premium='domestic'):
        '''Yields both prices hedged with spot and forward'''
        return ('Priced with Spot: {} ; Priced with Forward: {}'.format(self.spotPrice(premium=premium),self.forwardPrice(premium=premium)))
    
    def delta(self,spot=None,forward=None,hedge='spot'):
        if spot is None:
            spot = self.spot
        if forward is None:
            forward = self.forward
        if hedge.lower() == 'spot':
            d1 = BSpriceS(spot=spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d1()
            if self.flvr == 1:
                return exp(-self.frate*self.to_maturity()) * si.norm.cdf(d1,0,1)
            else:
                return exp(-self.frate*self.to_maturity()) * (si.norm.cdf(d1,0,1)-1)
        elif hedge.lower() == 'forward':
            d1 = BSpriceF(forward=forward,strike=self.strike,dt=self.to_maturity(),dr=self.drate,vol=self.vol,flvr=self.flvr).d1()
            if self.flvr == 1:
                return si.norm.cdf(d1,0,1)
            else:
                return si.norm.cdf(d1,0,1)-1
        else:
            raise ClassError('Please input either spot or forward as the hedge')
    
    def gamma(self):
        d1 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d1()
        return (si.norm.pdf(d1, 0,1)*exp(-self.frate)*self.to_maturity())/self.spot * self.vol * self.to_maturity()**0.5
    
    def vega(self):
        d1 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d1()
        return self.spot*exp(-self.frate*self.to_maturity()) * si.norm.pdf(d1,0,1) * self.to_maturity()**0.5

    def theta(self, flavour=None):
        if flavour is None:
            flavour = self.flvr
        d1 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d1()
        d2 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d2()
        if str(flavour) == '1':
            return -(self.spot*exp(-self.frate * self.to_maturity()) *si.norm.pdf(d1,0,1)*self.vol)/(2*self.to_maturity()**0.5) -  -self.frate*self.spot*exp(-self.frate*self.to_maturity())*si.norm.cdf(d1,0,1) - self.drate*self.strike*exp(-self.drate*self.to_maturity())*si.norm.cdf(d2,0,1)
        else:
            return -(self.spot*exp(-self.frate * self.to_maturity()) *si.norm.pdf(d1,0,1)*self.vol)/(2*self.to_maturity()**0.5) +  -self.frate*self.spot*exp(-self.frate*self.to_maturity())*si.norm.cdf(-d1,0,1) + self.drate*self.strike*exp(-self.drate*self.to_maturity())*si.norm.cdf(-d2,0,1)
    
    def rho(self,flavour=None):
        if flavour is None:
            flavour = self.flvr
        d2 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d2()
        if str(flavour) == '1':
            return self.to_maturity()*self.strike*exp(-self.drate*self.to_maturity()) * si.norm.cdf(d2,0,1)
        else:
            return -self.to_maturity()*self.strike*exp(-self.drate*self.to_maturity()) * si.norm.cdf(-d2,0,1)
    
    def zeta(self, flavour=None):
        if flavour is None:
            flavour = self.flvr
        d2 = BSpriceS(spot=self.spot,strike=self.strike,dt=self.to_maturity(),dr=self.drate,fr=self.frate,vol=self.vol,flvr=self.flvr).d2()
        if str(flavour) == '1':
            return si.norm.cdf(d2,0,1)
        else:
            return si.norm.cdf(-d2,0,1)
    
    def getGreeks(self,flavour=None):
        if flavour is None:
            flavour = self.flvr
        print('Delta : {}'.format(self.delta()) + '\nGamma: {}'.format(self.gamma()) + '\nVega: {}'.format(self.vega()) + '\nTheta: {}'.format(self.theta(flavour=flavour)) + '\nRho: {}'.format(self.rho(flavour=flavour)) + '\nZeta: {}'.format(self.zeta(flavour=flavour)))
    
    def plot_dvds(self, premium='domestic',withPayoff=False):
        dvds = [[],[]]
        for x in np.arange(0.75*self.strike,1.25*self.strike,0.01):
            dvds[0].append(self.spotPrice(premium=premium,spot=x))
            dvds[1].append(round(x,2))
        plt.plot(Series(dvds[0],dvds[1]))
        if withPayoff and premium=='domestic':
            plt.plot(dvds[1],[max(self.flvr*(x-self.strike),0)*self.notional for x in dvds[1]])
        plt.ylabel('Option Value in {}'.format(premium.capitalize()))
        plt.xlabel('Spot Value')
        plt.show()
        plt.clf()
    
    def plot_dvdt(self, premium='%domestic'):
        dvdt = [[],[]]
        for x in range(1,357):
            dvdt[0].append(self.spotPrice(expiry=datetime.date.today()+datetime.timedelta(days=x),premium=premium))
            dvdt[1].append(x)
        plt.plot(Series(dvdt[0],dvdt[1]))
        plt.ylabel('Option Value in {}'.format(premium.capitalize()))
        plt.xlabel('Days to maturity')
        plt.show()
        plt.clf()
    
    def plot_dvdvol(self,premium='%domestic'):
        dvdvol = [[],[]]
        for x in np.arange(0.001,0.3,0.01):
            dvdvol[0].append(self.spotPrice(vol=x,premium=premium))
            dvdvol[1].append(x)
        plt.plot(Series(dvdvol[0],dvdvol[1]))
        plt.ylabel('Option Value in {}'.format(premium.capitalize()))
        plt.xlabel('Volatility Value')
        plt.show()
        plt.clf()
    
    def plot_delta(self):
        delta = [[],[]]
        for x in np.arange(0.75*self.strike,1.25*self.strike,0.01):
            delta[0].append(self.delta(spot=x))
            delta[1].append(x)
        plt.plot(Series(delta[0],delta[1]))
        plt.xlabel('Delta')
        plt.show()
        plt.clf()