from VanillaClass import Vanilla,ClassError
from BarrierClass import BarrierContingentVanilla,BarrierContingentPayment

class Option():

    '''
    Returns an instance of a derivative contract based on the specified arguments
    '''
    def __new__(self,**kwargs):
        self.cpair='Not set'
        self.notional=100000
        self.forward = None
        barrierLoc = ['upper','lower','down','up']
        validKwargs = ['spot', 'forward', 'strike', 'expiry', 'vol',
                       '_type_', 'flvr', 'frate', 'drate', 'notional', 'cpair', 'barrier']
        validKwargs += [loc+'barrier' for loc in barrierLoc]+[loc+'_barrier' for loc in barrierLoc]
        for k in kwargs.keys():
            if k.lower() in validKwargs:
                continue
            else:
                raise ClassError('{} is not a valid parameter, please specify one of the following: \n{}'.format(k,validKwargs))
        for k,v in kwargs.items():
            if k.lower() in ['lowerbarrier', 'lower', 'down', 'lower_barrier']:
                self.lower = kwargs[k]
            elif k.lower() in ['upperbarrier', 'upper', 'up', 'upper_barrier']:
                self.upper=kwargs[k]
            else:
                setattr(self,k.lower(),v)

        if hasattr(self, 'barrier') and hasattr(self, 'strike'):
            return BarrierContingentVanilla(spot=self.spot,
                                            strike=self.strike,
                                            expiry=self.expiry,
                                            drate=self.drate,
                                            frate=self.frate,
                                            vol=self.vol,
                                            barrier=self.barrier,
                                            Type=self._type_,
                                            flvr=self.flvr,
                                            cpair=self.cpair,
                                            notional=self.notional)

        elif hasattr(self,'barrier') and not hasattr(self, 'strike'):
            return BarrierContingentPayment(S=self.spot,
                                            fr=self.frate,
                                            dr=self.drate,
                                            vol=self.vol,
                                            expiry=self.expiry,
                                            Type=self._type_,
                                            notional=self.notional,
                                            barrier=self.barrier)

        elif (hasattr(self, 'upper') and hasattr(self, 'lower')) and not hasattr(self, 'strike'):
            return BarrierContingentPayment(S=self.spot,
                                            fr=self.frate,
                                            dr=self.drate,
                                            vol=self.vol,
                                            expiry=self.expiry,
                                            Type=self._type_,
                                            notional=self.notional,
                                            lowerbarrier=self.lower,
                                            upperbarrier=self.upper)

        elif hasattr(self, 'strike') and \
                hasattr(self, 'flvr') and not \
                (hasattr(self, 'barrier') or hasattr(self, 'upper') or hasattr(self, 'lower')):

            return Vanilla(spot=self.spot,
                           strike=self.strike,
                           expiry=self.expiry,
                           drate=self.drate,
                           frate=self.frate,
                           vol=self.vol,
                           flvr=self.flvr,
                           notional=self.notional,
                           forward=self.forward)
        else:
            raise ClassError('Specified arguments do not match any type of an option instance')