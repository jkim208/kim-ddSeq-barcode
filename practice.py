import pyximport  # install Cython to access pyximport
pyximport.install(build_in_temp=False)
from editDistance import edit_distance
import distance


stringBase = 'JOHNHANCOCK'
stringTest = 'JOHNHANCOCK'
print('For Cython')
print(edit_distance(stringBase, stringTest))

print('For distance pack')
print(distance.hamming(stringBase, stringTest))





