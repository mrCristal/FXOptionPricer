import datetime

class Contract:
    def __init__(self, spot, strike, expiry, drate, notional=100000,cpair='Not Set'):
        self.spot = spot
        self.strike = strike
        self.expiry = expiry
        self.drate = drate
        self.notional = notional
        self.cpair = cpair
    
    def to_maturity(self, date=None):
        if date is None:
            return (self.expiry - datetime.date.today()).days/365
        return (date - datetime.date.today()).days/365
    