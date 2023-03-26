# Try to practice the OOP

class Reactor:
    tipo = "BATCH"

    def __init__(self,T,N,P) -> None:
        self.T = T
        self.N = N
        self.P = P
    
    # Instance method
    def time(self):
        N = self.N
        T = self.T
        P = self.P
        t = N/T*P
        return t


R1 = Reactor(365,2,9)
R2 = Reactor(285,6,4)

R1.time()
R2.time()