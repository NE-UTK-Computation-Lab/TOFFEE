import time
import pickle



##### name for saved dictionary
name = 'jezebel_8'


time_1 = time.perf_counter()


with open('temp','rb') as file:
    main = pickle.load(file) 


##### These are the specifications for Tally identification. This does nothing for KEFF cases.
main.tally_identifier = '1tally        4' 
main.location_identifier = 'cell  100' 
main.tally_multiplier = False
main.multiplier_identifier = 'multiplier bin:   1.00000E-24'
main.tally_bin = False
main.erg_bin_identifier = 'total'


main.evaluate_variance()

#main.check_spd()
#main.super_matrix()

with open(name+'.pk','wb') as file:
    pickle.dump(main,file,pickle.HIGHEST_PROTOCOL)   
    
    
  
with open(name+'.pk','rb') as file:
    main = pickle.load(file) 

main.print_total_variance()

#main.plot_sensitivity('dE')
### options:
    #'dU' - per unit lethargy 
    #'dE' - per MeV

main.plot_uncertainties()


#main.plot_variance('std')
### options:
    #'var' - variance 
    #'std' - standard deviation

print('done in ' + str(round(time.perf_counter()-time_1,3)) + ' seconds')