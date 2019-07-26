import scipy.stats as si
from math import log,exp,cos,pi
import numpy as np


class ClassError(Exception):
    pass


class BSpriceS():
    def __init__(self,spot,strike,dt,dr,fr,vol,flvr):
        self.S = spot
        self.K = strike
        self.dr = dr
        self.fr = fr
        self.vol = vol
        self.dt = dt
        self.flvr = flvr

    def d1(self,S=None,K=None,dr=None,fr=None,dt=None,vol=None):
        if S is None: S = self.S
        if K is None: K = self.K
        if dr is None: dr = self.dr
        if fr is None: fr = self.fr
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        return (log(S/K)+(dr-fr+0.5*vol**2)*dt)/(vol*(dt**0.5))
    
    def d2(self,S=None,K=None,dr=None,fr=None,dt=None,vol=None):
        if S is None: S = self.S
        if K is None: K = self.K
        if dr is None: dr = self.dr
        if fr is None: fr = self.fr
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        return (log(S/K)+(dr-fr-0.5*vol**2)*dt)/(vol*(dt**0.5))
    
    def price(self,S=None,K=None,dr=None,fr=None,dt=None,vol=None,flvr=None):
        if S is None: S = self.S
        if K is None: K = self.K
        if dr is None: dr = self.dr
        if fr is None: fr = self.fr
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        if flvr is None: flvr = self.flvr
        Nd1 = si.norm.cdf(flvr*self.d1(),0,1)
        Nd2 = si.norm.cdf(flvr*self.d2(),0,1)
        return flvr*(S*exp(-1*fr*dt)*Nd1 - K*exp(-1*dr*dt)*Nd2)

class BSpriceF():
    def __init__(self,forward,strike,dt,dr,vol,flvr):
        self.F = forward
        self.K = strike
        self.dr = dr
        self.vol = vol
        self.dt = dt
        self.flvr = flvr

    def d1(self,F=None,K=None,dt=None,vol=None):
        if F is None: F = self.F
        if K is None: K = self.K
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        return (log(F/K)+(0.5*vol**2)*dt)/(vol*(dt**0.5))
    
    def d2(self,F=None,K=None,dt=None,vol=None):
        if F is None: F = self.F
        if K is None: K = self.K
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        return (log(F/K)-(0.5*vol**2)*dt)/(vol*(dt**0.5))
    
    def price(self,F=None,K=None,dr=None,dt=None,vol=None,flvr=None):
        if F is None: F = self.F
        if K is None: K = self.K
        if dt is None: dt = self.dt
        if dr is None: dr = self.dr
        if vol is None: vol = self.vol
        if flvr is None: flvr = self.flvr
        Nd1 = si.norm.cdf(flvr*self.d1(F=F,K=K,vol=vol,dt=dt),0,1)
        Nd2 = si.norm.cdf(flvr*self.d2(F=F,K=K,vol=vol,dt=dt),0,1)
        return flvr*exp(-1*dr*dt)*(F*Nd1 - K*Nd2)

class MCprice():
    def __init__(self,S,fr,dr,vol,dt,K,flvr,sims=10000):
        self.S    = S
        self.fr   = fr
        self.dr   = dr
        self.vol  = vol
        self.dt   = dt
        self.K    = K
        self.sims = sims
        self.flvr = flvr

    def HaltonNrs(self,n,b):
        n0, H, f = n, 0, 1/b
        while n0 > 0:
            n1 = int(n0/b)
            r = n0 - n1*b
            H = H+f*r
            f = f/b
            n0 = n1
        return H
    
    def price(self, S=None, dr=None, fr=None, vol=None, dt=None, K=None, flvr=None, sims=10000, withLog=False):
        if S is None: S = self.S
        if K is None: K = self.K
        if dr is None: dr = self.dr
        if fr is None: fr = self.fr
        if dt is None: dt = self.dt
        if vol is None: vol = self.vol
        if flvr is None: flvr = self.flvr
        sumpayOff = 0
        drift = (dr-fr - 0.5*vol**2)*dt
        vSqrdt = vol*dt**0.5
        if withLog:
            X=log(S)
            for sim in range(1,sims+1):
                eps=self.BoxMuller(self.HaltonNrs(sim,3),self.HaltonNrs(sim,5))
                Xt1 = X + drift + vSqrdt*eps
                Xt2 = X + drift + vSqrdt*eps*-1
                sumpayOff += (max(flvr*(exp(Xt1)-K),0) + max(flvr*(exp(Xt2)-K),0))/2
        else:
            for sim in range(1,sims+1):
                eps=self.BoxMuller(self.HaltonNrs(sim,3),self.HaltonNrs(sim,5))
                St1 = S*exp(drift+vSqrdt*eps)
                St2 = S*exp(drift+vSqrdt*eps*-1)
                sumpayOff+=(max(flvr*(St1-K),0)+max(flvr*(St2-K),0))/2
        return exp(-1*dr*dt)*sumpayOff/sims
    
    def BoxMuller(self,z,y):
        return (-2*log(z))**0.5 * cos(2*pi*y)


class MCpricePayment():
    def __init__(self,S,fr,dr,vol,dt,Type,sims,**kwargs):
        self.S=S
        self.fr=fr
        self.dr=dr
        self.vol=vol
        self.dt=dt
        self.sims=sims
        self.payoff=self.__formatType__(Type)
        self.Type = Type
        for k,v in kwargs.items():
            setattr(self,k,v)
    
    def __formatType__(self,paymentType):
        if paymentType.lower() not in ['touch','notouch']:
            raise ClassError('Please specify touch or notouch')
        if paymentType.lower() == 'touch':
            return 1
        elif paymentType.lower() == 'notouch':
            return -1

    def price(self,steps):
        
        step = self.dt/steps

        if self.payoff == 1:
            sumpayOff=0
        elif self.payoff==-1:
            sumpayOff=2*self.sims # *2 because I am generating 2 paths for the antithetic variables
        else:
            raise ClassError('Something went wrong')

        drift = (self.dr-self.fr-0.5*self.vol**2)*step
        vSqrdt = self.vol*step**0.5
        logSpot = log(self.S)
        
        if hasattr(self,'barrier'):
            above = self.barrier > self.S
            log_barrier=log(self.barrier)
            for _sim in range(self.sims):
                X = logSpot
                X1 = logSpot
                Xpath = []
                X1path = []
                for _i in range(steps):
                    eps = np.random.normal(0,1)
                    dX = drift+vSqrdt*eps
                    dX1 = drift + vSqrdt*eps*-1
                    X += dX
                    X1 += dX1
                    Xpath.append(X)
                    X1path.append(X1)
                if above:
                    if max(Xpath) >= log_barrier:
                            sumpayOff += self.payoff
                    
                    if max(X1path) >= log_barrier:
                            sumpayOff += self.payoff
                else:
                    if min(Xpath) <= log_barrier:
                            sumpayOff += self.payoff
                    
                    if min(X1path) <= log_barrier:
                            sumpayOff += self.payoff
        else:
            upperLog = log(self.upper)
            lowerLog = log(self.lower)
            for _sim in range(self.sims):
                X = logSpot
                X1 = logSpot
                Xpath = []
                X1path = []
                for _i in range(steps):
                    eps = np.random.normal(0,1)
                    dX = drift+vSqrdt*eps
                    dX1 = drift + vSqrdt*eps*-1
                    X += dX
                    X1 += dX1
                    Xpath.append(X)
                    X1path.append(X1)
                if max(Xpath) >= upperLog or min(Xpath) <=lowerLog:
                            sumpayOff += self.payoff
                    
                if max(X1path) >= upperLog or min(X1path) <=lowerLog:
                            sumpayOff += self.payoff


        return exp(-1*self.dr*self.dt)*sumpayOff/(self.sims*2)

    def priceWithBB(self,steps):
        
        step = self.dt/steps #step is effectively 1/365 , can be seen if you check the self.dt computation as well as for steps in BarrierClass.Payment.price
        if self.payoff == 1:
            sumpayOff=0
        elif self.payoff == -1:
            sumpayOff=self.sims
        else:
            raise ClassError('Something went wrong')
        
        drift = (self.dr-self.fr-0.5*self.vol**2)*step
        vSqrdt = self.vol*step**0.5

        boundProb = lambda x : x if 0<=x<=1 else (0 if x<0 else 1) 

        def probTrigger(x,y,z,):
            prob = round(exp(-2*log(x/z)*log(y/z)/(step*self.vol**2)),5)
            return boundProb(prob)

        if hasattr(self,'barrier'):
            above = self.barrier > self.S
            for _sim in range(self.sims):
                X = self.S
                Xpath = [X]
                for _i in range(steps):
                    eps = np.random.normal(0,1)
                    dX = drift+vSqrdt*eps
                    X *= exp(dX)
                    Xpath.append(X)
                    if above:
                        if X>=self.barrier:
                            sumpayOff+=self.payoff
                            break
                        else:
                            if round(np.random.random(),5) <= probTrigger(Xpath[-1],Xpath[-2],self.barrier):
                                sumpayOff+=self.payoff
                                break
               
                    else:
                        if X<=self.barrier:
                            sumpayOff+=self.payoff
                            break
                        else:
                            if round(np.random.random(),5) <= probTrigger(Xpath[-1],Xpath[-2],self.barrier):
                                sumpayOff+=self.payoff
                                break
        else:
            for _sim in range(self.sims):
                X = self.S
                Xpath = [X]
                for _i in range(steps):
                    eps = np.random.normal(0,1)
                    dX = drift+vSqrdt*eps
                    X *= exp(dX)
                    Xpath.append(X)
                    if X>=self.upper or X<=self.lower:
                        sumpayOff+= self.payoff
                        break
                    elif round(np.random.random(),5) <= probTrigger(Xpath[-1],Xpath[-2],self.upper) or round(np.random.random(),5) <= probTrigger(Xpath[-1],Xpath[-2],self.lower):
                        sumpayOff+=self.payoff
                        break
        return exp(-1*self.dr*self.dt)*sumpayOff/(self.sims)

class BarrierAnalyticalPrice():
    def __init__(self,spot,strike,drate,frate,vol,barrier,_type,nu,flvr,mu,tau):
        self.spot = spot
        self.strike = strike
        self.drate = drate
        self.frate = frate
        self.vol = vol
        self.flvr = flvr
        self.barrier = barrier
        self.tau = tau
        self.mu = mu
        self.lbda = (self.mu**2 + 2*self.drate/self.vol**2)**0.5
        self.nu = nu
        self._type_ = _type
    
    def __x1__(self,S,K,vol,tau,mu):
        return log(S/K)/(vol*tau**0.5) + (1+mu)*vol*tau**0.5
    
    def __x2__(self,S,H,vol,tau,mu):
        return log(S/H)/(vol*tau**0.5) + (1+mu)*vol*tau**0.5

    def __y1__(self,H,S,K,vol,tau,mu):
        return log(H**2/(S*K))/(vol*tau**0.5) + (1+mu)*vol*tau**0.5

    def __y2__(self,H,S,vol,tau,mu):
        return log(H/S)/(vol*tau**0.5) + (1+mu)*vol*tau**0.5

    def __A__(self,flvr,S,fr,dr,tau,K,vol,mu):
        return flvr*S*exp(-1*fr*tau) * si.norm.cdf(flvr*self.__x1__(S,K,vol,tau,mu),0,1)-flvr*K*exp(-1*dr*tau)*si.norm.cdf(flvr*self.__x1__(S,K,vol,tau,mu)-flvr*vol*tau**0.5,0,1)
    
    def __B__(self,flvr,S,H,fr,dr,tau,K,vol,mu):
        return flvr*S*exp(-1*fr*tau) * si.norm.cdf(flvr*self.__x2__(S,H,vol,tau,mu),0,1)-flvr*K*exp(-1*dr*tau)*si.norm.cdf(flvr*self.__x2__(S,H,vol,tau,mu)-flvr*vol*tau**0.5,0,1)

    def __C__(self,flvr,S,fr,dr,tau,K,vol,mu,H,nu):
        return flvr*S*exp(-1*fr*tau)*(H/S)**(2*(mu+1))*si.norm.cdf(nu*self.__y1__(H,S,K,vol,tau,mu),0,1)-flvr*K*exp(-1*dr*tau)*(H/S)**(2*mu)*si.norm.cdf(nu*self.__y1__(H,S,K,vol,tau,mu)-nu*vol*tau**0.5,0,1)
    
    def __D__(self,flvr,S,fr,dr,tau,K,vol,mu,H,nu):
        return flvr*S*exp(-1*fr*tau)*(H/S)**(2*(mu+1)) * si.norm.cdf(nu*self.__y2__(H,S,vol,tau,mu),0,1)-flvr*K*exp(-1*dr*tau)*(H/S)**(2*mu)*si.norm.cdf(nu*self.__y2__(H,S,vol,tau,mu)-nu*vol*tau**0.5,0,1)
    
    def price(self,flvr=None,S=None,fr=None,dr=None,tau=None,K=None,vol=None,mu=None,H=None,nu=None,barrierType=None,premium='domestic'):
        if flvr is None: flvr = self.flvr
        if   S  is None: S   =  self.spot
        if   H  is None: H  =   self.barrier
        if  fr  is None: fr  =  self.frate
        if  dr  is None: dr  =  self.drate
        if  tau is None: tau =  self.tau
        if   K  is None: K  =   self.strike
        if  vol is None: vol =  self.vol
        if  mu  is None: mu  =  self.mu
        if  nu  is None: nu  =  self.nu
        if barrierType is None: barrierType = self._type_
        
        if barrierType.lower() == 'cin' or barrierType.lower() == 'down and in call':
            if S<=H:
                return BSpriceS(spot=S,strike=K,dt=tau,dr=dr,fr=fr,vol=vol,flvr=flvr).price()
            if K>H: 
                return self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) + self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
        
        elif barrierType.lower() == 'cni' or barrierType.lower() == 'up and in call':
            if S>=H:
                return BSpriceS(spot=S,strike=K,dt=tau,dr=dr,fr=fr,vol=vol,flvr=flvr).price()
            if K>H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu)
            elif K<H: 
                return self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) - self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu) + self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
        
        elif barrierType.lower() == 'pin' or barrierType.lower() == 'down and in put':
            if S<=H:
                return BSpriceS(spot=S,strike=K,dt=tau,dr=dr,fr=fr,vol=vol,flvr=flvr).price()
            if K>H: 
                return self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) - self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu) + self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu)
        
        elif barrierType.lower() == 'pni' or barrierType.lower() == 'up and in put':
            if S>=H:
                return BSpriceS(spot=S,strike=K,dt=tau,dr=dr,fr=fr,vol=vol,flvr=flvr).price()
            if K>H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) + self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H: 
                return self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
        
        elif barrierType.lower() == 'con' or barrierType.lower() == 'down and out call':
            if S<=H:
                return 0
            if K>H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H: 
                return self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) - self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
        
        elif barrierType.lower() == 'cno' or barrierType.lower() == 'up and out call':
            if S>=H:
                return 0
            if K>H: 
                return 0
            elif K<H:
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) + self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu) - self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
        
        elif barrierType.lower() == 'pon' or barrierType.lower() == 'down and out put':
            if S<=H:
                return 0
            if K>H: 
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) + self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu) - self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H:
                return 0
        
        elif barrierType.lower() == 'pno' or barrierType.lower() == 'up and out put':
            if S>=H:
                return 0
            if K>H: 
                return self.__B__(flvr,S,H,fr,dr,tau,K,vol,mu) - self.__D__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)
            elif K<H:
                return self.__A__(flvr,S,fr,dr,tau,K,vol,mu) - self.__C__(flvr,S,fr,dr,tau,K,vol,mu,H,nu)