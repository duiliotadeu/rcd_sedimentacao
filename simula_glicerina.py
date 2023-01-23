# arquivo para execucao das simulacoes iniciais com glicerina


import subprocess

print('Rodando python script para multiplas simulações:')
print('================================================')

casos = ['glicerina2_3particula_15p.in', 'glicerina2_3particula_10p.in', \
    'glicerina3_1particula_25p.in', 'glicerina3_1particula_20p.in', 'glicerina3_1particula_15p.in', 'glicerina3_1particula_10p.in', 
    'glicerina3_2particula_25p.in', 'glicerina3_2particula_20p.in',
    'glicerina3_2particula_15p.in', 'glicerina3_2particula_10p.in', 
    'glicerina3_3particula_25p.in', 'glicerina3_3particula_20p.in',
    'glicerina3_3particula_15p.in', 'glicerina3_3particula_10p.in',]

for c in casos:
    print('Executando Caso ' + c)
    comando = './sedimentacao -i '+ c
    subprocess.run([comando], shell=True)
    print('*****  Caso : '+c +' executado')