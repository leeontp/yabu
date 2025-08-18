import time
import os

fox_banner = """
               /\\_/\\      
              ( o.o )     Y A B U
               > ^ <      Your All-purpose Bio-structural Utility
   ==================================================================
                           W E L C O M E
   ==================================================================
"""

fox_confused = """
               /\\_/\\   
              ( o.O )?        YABU is confused...
               > ^ <  
"""

# Animación del zorrito comiendo
fox_eating_frames = [
    r"""
               /\_/\      
              ( o.o )      
              ( >🍪 )      YABU is eating...
    """,
    r"""
               /\_/\      
              ( -.- )      
              ( 🍪< )      *munch munch*
    """,
    r"""
               /\_/\      
              ( o.o )      
              ( >🍪 )      YABU is eating...
    """,
    r"""
               /\_/\      
              ( -.- )      
              (  <  )      *gulp!* 
    """
]

def animate_fox_eating(cycles: int = 2, delay: float = 0.5):
    """
    Animación del zorrito comiendo.
    
    Parameters
    ----------
    cycles : int
        Número de veces que repite la animación.
    delay : float
        Tiempo entre frames.
    """
    for _ in range(cycles):
        for frame in fox_eating_frames:
            os.system("clear")
            print(frame)
            time.sleep(delay)

