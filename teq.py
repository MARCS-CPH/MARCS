import numpy as np

r_sun_AU = 0.00465047

r_star = 1.01 #in solar radii
t_star = 5770.0
d = 0.06 #AU

teq = t_star*((r_star*r_sun_AU) /(2.*d))**0.5

print(teq)