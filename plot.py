import time
import pickle



##### name for saved dictionary
name = 'jezebel_7_0'


time_1 = time.perf_counter()


with open(name,'rb') as file:
    main = pickle.load(file) 


main.plot_matrix()


print('done in ' + str(round(time.perf_counter()-time_1,3)) + ' seconds')