import time
import pickle



##### name for saved dictionary
name = 'jezebel_8'


time_1 = time.perf_counter()


with open('temp','rb') as file:
    main = pickle.load(file) 
main.evaluate_variance()

with open(name,'wb') as file:
    pickle.dump(main,file,pickle.HIGHEST_PROTOCOL)   
    
    
  
with open(name,'rb') as file:
    main = pickle.load(file) 

main.print_total_variance()

main.plot_sensitivity('dU')
### options:
    #'dU' - per unit lethargy 
    #'dE' - per MeV

main.plot_matrix()


main.plot_variance('std')
### options:
    #'var' - variance 
    #'std' - standard deviation

print('done in ' + str(round(time.perf_counter()-time_1,3)) + ' seconds')
