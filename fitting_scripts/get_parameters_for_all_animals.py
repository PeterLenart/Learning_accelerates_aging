import pandas as pd
import numpy as np
from scipy.optimize import minimize
import os
import cma
import csv

matplotlib.rcParams.update({'font.size': 17})
matplotlib.rcParams['figure.figsize'] = [15, 15]

# Change to the folder where your animal mortality data are saved
folders_path = "./mortality_data/animals_and_plants/animals_for_parameters/"

# Change to the file where you want to save the parameters
params_csv = "./parameters_animals.csv"

folders = sorted([os.path.join(folders_path, x) for x in os.listdir(folders_path)])
# Get the list of all studied species
animals = []
for folder in folders:
    files = sorted([os.path.join(folder, x) for x in os.listdir(folder)])
    for file in files :
        data_path = file
        animal = os.path.splitext(os.path.basename(data_path))[0]
        latin_name = ' '.join(animal.split(' ')[0:2])
        animals.append(animal)


# function computing the contribution of aging to the hazard rate
@np.vectorize
def aging_gompertz_makaham(x, a, b, c):
    return c + a*np.exp((b*x))

# function computing the contribution of learning to the hazard rate
@np.vectorize
def learning(x, Lmax, k_learning, n):
    return (Lmax/(1 + np.exp(n*(x-k_learning)))) - Lmax

# function computing the contribution of growth to the hazard rate
@np.vectorize
def growth(x, Gmax, growth_speed, age_stop):
    results = []
    if (x < age_stop):
        results.append(Gmax/(1 + x**growth_speed) - Gmax ) 
    else : 
        results.append(Gmax/(1 + age_stop**growth_speed) - Gmax ) 
    return np.array(results)

# function summing all the contributions and returning the log of this sum. 
@np.vectorize
def log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed, age_stop):
    res = aging_gompertz_makaham(x, a, b, c) + learning(x, Lmax, k_learning, n) + growth(x, Gmax, growth_speed, age_stop)
    if res > 0 : 
        return np.log10(res)
    else :
        return np.log10(1.)

# cost function used for fitting
def errlearning_gompertz_makaham(params, x, y, weights, lm=0, ord=0):
    from numpy.linalg import norm
    a, b, c, Lmax, k_learning, n, Gmax, growth_speed, age_stop = params

    # penalization for negative values
    if any(e<0 for e in params) == True:
        return 100

    err = norm(weights*(log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed, age_stop) - y))
    return err + lm * norm(params, ord)

def errlearning_gompertz_makaham_cma(params, x, y, weights, bounds, lm=0, ord=0):
    from numpy.linalg import norm
    a, b, c, Lmax, k_learning, n, Gmax, growth_speed, age_stop = params

    penality = 0
    #checks if params are in bound
    for i in range(len(params)):
        if params[i] > bounds[i][1] : 
            penality += (1+ params[i] - bounds[i][1])**2
        if params[i] < bounds[i][0]:
            penality += (2 + bounds[i][0] - params[i])**2


    err = norm(weights*(log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed, age_stop) - y))
    return err + lm * norm(params, ord) + penality

# function used for computing the s measure (average distance between the data and the fitted model in %)
def s_measure(pred, data):
    distance = np.absolute(pred - data)
    distance_percentage = (distance/np.absolute(data))*100
    return np.mean(distance_percentage)
    # return np.average(distance_percentage, weights=weights)


# fitting function using the minimize method from scipy
def fit(x, mortality_array, guess, bounds):

    # weight computation, can sometimes improve the fit but may also impact it negatively
    lx = [10000]
    for i in range(1, mortality_array.shape[0]):
        lx.append(lx[i-1]-lx[i-1]*mortality_array[i])
    lx = np.array(lx)
    surviving_fraction = lx/lx[0]
    surviving_fraction = np.square(surviving_fraction)

    # comment this line to use the computed weights
    surviving_fraction = np.ones_like(surviving_fraction)

    params = minimize(errlearning_gompertz_makaham, guess, method="SLSQP", tol = 1e-15, args=(x, np.log10(mortality_array), surviving_fraction, 1e-4, 2), bounds = bounds)
    return params.x

# fitting function using the CMA method, more accurate but slower than the one using minimize from scipy
def fit_cma(x, mortality_array, sigma, bounds):
    lx = [10000]
    for i in range(1, mortality_array.shape[0]):
        lx.append(lx[i-1]-lx[i-1]*mortality_array[i])
    lx = np.array(lx)
    surviving_fraction = lx/lx[0]
    surviving_fraction = np.square(surviving_fraction)
    surviving_fraction = np.ones_like(surviving_fraction)
    global f 
    def f(params):
        return errlearning_gompertz_makaham_cma(params, x, np.log10(mortality_array), surviving_fraction, bounds)

    es = cma.CMAEvolutionStrategy([4.76614499e-03, 2.29715390e-02, 1.87466206e-02, 5.58711082e-03, 3.63829644e+01, 1.38484469e-01, 4.34710796e-02, 0.1, 30], sigma)
    es.optimize(f, iterations=5000, n_jobs=-1)
    return es

def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)
        
if __name__ == "__main__":
    
    # Create the csv file where the parameters for this country will be saved and write the header
    file = open(params_csv, 'a+', newline='')
    file.close()
    header = ['group', 'name', 'a', 'b', 'c', 'Lmax', 'k_learning', 'n', 'Gmax', 'growth_speed', 'age_stop']
    append_list_as_row(params_csv, header)
    
    # Loop over the folders contained in the main data folder (one folder for one group of species)
    for folder in folders:
        files = sorted([os.path.join(folder, x) for x in os.listdir(folder)])
        group = os.path.basename(folder)
        # Loop over all the different species
        for file in files :
            
            # Initialize the output row
            row = []

            data_path = file
            animal = os.path.splitext(os.path.basename(data_path))[0]
            name = ' '.join(animal.split(' ')[0:-1])

            # Append the group and the specie's name to the output row
            row.append(group)
            row.append(name)

            # read text file into pandas DataFrame
            base_dataframe = pd.read_csv(data_path, sep=",")
            base_dataframe = base_dataframe[['x', 'qx']]
            df1 = base_dataframe.astype(float)

            years = df1[['x']].to_numpy().squeeze()
            mortality_array = df1[['qx']].to_numpy().squeeze()

            beg = np.min(years) #starting x 
            fin = np.max(years) #final x
            steps = years.size
            x = np.linspace(beg, fin, steps)
            
            # Initial guess and bounds for the parameters
            guess = [4.17399393e-03, 4.29746495e-02, 2.86076474e-02, 3.91364426e-02, 4.61448543e+01, 8.85350563e-02, 5.92916574e-02, 0.1, 30]
            bounds = ((0,2),(0,2),(0,2),(0,1),(0,fin),(0,3),(0,1),(0, 1),(0, fin))

            # Fit with both fitting methods (CMA and scipy.minimize)
            es = fit_cma(x, mortality_array, 0.05, bounds)
            params = es.result.xfavorite
            params2 = fit(x, mortality_array, guess, bounds)

            sm1 = s_measure(log_mortality_gompertz_makaham(x, params[0], params[1], params[2],  params[3], params[4], params[5], params[6], params[7], params[8]), np.log10(mortality_array))
            sm2 = s_measure(log_mortality_gompertz_makaham(x, params2[0], params2[1], params2[2],  params2[3], params2[4], params2[5], params2[6], params2[7], params2[8]), np.log10(mortality_array))

            # Keep the parameters of the best fit and save them
            if sm1 < sm2 : 
                for idx in range(params.shape[0]):
                    row.append(params[idx])
            else:
                for idx in range(params2.shape[0]):
                    row.append(params2[idx])

            append_list_as_row(params_csv, row)
    
    

