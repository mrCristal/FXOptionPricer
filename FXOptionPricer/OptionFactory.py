from VanillaClass import Vanilla,ClassError
from BarrierClass import BarrierContingentVanilla,BarrierContingentPayment

def Option(**kwargs):

    '''
    Returns an instance of a derivative contract based on the specified arguments
    '''
    if 'notional' not in kwargs.keys():
        kwargs['notional'] = 100000
    if 'forward' not in kwargs.keys():
        kwargs['forward'] = None
    barrierLoc = ['upper','lower','down','up']
    validKwargs = ['spot', 'forward', 'strike', 'expiry', 'vol',
                    '_type_', 'flvr', 'frate', 'drate', 'notional', 'cpair', 'barrier']
    validKwargs += [loc+'barrier' for loc in barrierLoc]+[loc+'_barrier' for loc in barrierLoc]
    for k in kwargs.keys():
        if k.lower() in validKwargs:
            continue
        else:
            raise ClassError('{} is not a valid parameter, please specify one of the following: \n{}'.format(k,validKwargs))

    if 'barrier' in kwargs.keys() and 'strike' in kwargs.keys():
        return BarrierContingentVanilla(spot=kwargs['spot'],
                                        strike=kwargs['strike'],
                                        expiry=kwargs['expiry'],
                                        drate=kwargs['drate'],
                                        frate=kwargs['frate'],
                                        vol=kwargs['vol'],
                                        barrier=kwargs['barrier'],
                                        Type=kwargs['_type_'],
                                        flvr=kwargs['flvr'],
                                        cpair=kwargs['cpair'],
                                        notional=kwargs['notional'])

    elif 'barrier' in kwargs.keys() and 'strike' not in kwargs.keys():
        return BarrierContingentPayment(S=kwargs['spot'],
                                        fr=kwargs['frate'],
                                        dr=kwargs['drate'],
                                        vol=kwargs['vol'],
                                        expiry=kwargs['expiry'],
                                        Type=kwargs['_type_'],
                                        notional=kwargs['notional'],
                                        barrier=kwargs['barrier'])

    elif 'upper' in kwargs.keys() or 'lower' in kwargs.keys() and 'strike' not in kwargs.keys():
        return BarrierContingentPayment(S=kwargs['spot'],
                                        fr=kwargs['frate'],
                                        dr=kwargs['drate'],
                                        vol=kwargs['vol'],
                                        expiry=kwargs['expiry'],
                                        Type=kwargs['_type_'],
                                        notional=kwargs['notional'],
                                        lowerbarrier=kwargs['lower'],
                                        upperbarrier=kwargs['upper'])

    elif 'strike' in kwargs.keys() and 'flvr' in kwargs.keys()  and ('barrier' not in kwargs.keys() or 'upper' not in kwargs.keys() or 'lower' in kwargs.keys()):

        return Vanilla(spot=kwargs['spot'],
                        strike=kwargs['strike'],
                        expiry=kwargs['expiry'],
                        drate=kwargs['drate'],
                        frate=kwargs['frate'],
                        vol=kwargs['vol'],
                        flvr=kwargs['flvr'],
                        notional=kwargs['notional'],
                        forward=kwargs['forward'])
    else:
        raise ClassError('Specified arguments do not match any type of an option instance')
