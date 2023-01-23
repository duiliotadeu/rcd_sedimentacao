
import subprocess
print("Python script of multiple simulations")

simulation_name = ["-i caso_a_3particulas.in", "-i caso_a_8particulas.in", "-i caso_a_10particulas.in", \
    "-i caso_aLongo_3particulas.in", "-i caso_aLongo_10particulas.in"]

for s in simulation_name:
    cmd = './sedimentacao '+s
    print('running case  ', s)
    subprocess.run([cmd], shell=True)


