from abc import ABC, abstractmethod


class MonteCarloSim(ABC):
    
    @abstractmethod
    def set_par(self):
        pass
    
    @abstractmethod
    def f(self):
        pass
    
    @abstractmethod
    def t_span(self):
        pass
    
    @abstractmethod
    def set_IC(self):
        pass
    #Otrzymuje parametry, równania, tspan, IC
    #Jak potrzeba losuje par i IC <--- to ręcznie ustawiane?
    #Wywołuje obie SameImportedClass przesyłając par, f, tspan, IC 