# Frecuencias alélicas
m = (2*144 + 201)/(2*459)
n = (2*114 + 201)/(2*459)

# Consideramos M = p y N = q

# Frecuencias genotípicas
p_square = m^2
p_q_2 = (2*m*n)
q_square = n^2
frecuencias_genoripicas = p_square + p_q_2 + q_square

# Cálculo de x cuadrada (gi cuadrada)

gi_square_MM = (144 - 130)^2/130 
gi_square_MN = (201 - 228)^2/228
gi_square_NN = (114 - 100)^2/100

GI_SQUARE = gi_square_MM + gi_square_MN + gi_square_NN



# EJERCICIO II

# frecuencias alelicas

f = (2*320 + 80 + 125)/(2*1000)
M = (2*120 + 80 + 120)/(2*1000)
S = (2*235 + 125 + 120)/(2*1000)

# Freciencia Genotípica
# Consideramos f = p, M = r y S = q

frecuencias_genoripica2 = f^2 + S^2 + M^2 + (2*f*S) + (2*S*M) + (2*f*M)


# Indice de fijación 

He = 1 - (f^2 + M^2 + S^2)
# Frecuencias genótipicas observadas
Ho = (80 + 125 + 120)/1000

index_fijation = (He - Ho)/He
