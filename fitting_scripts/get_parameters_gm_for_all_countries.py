import pandas as pd
import numpy as np
from scipy.optimize import minimize
import os
import scipy.integrate as integrate
import csv

# Male, Female or Total, studied gender of the mortality data
gender = "Total"

# cohort or population, type of study
data_type = "population"

# Change to the folder when the mortality data is saved
data_folder = "./mortality_data/new_population_all_countries/"

# List all the files in this folder
data_files = [os.path.join(data_folder, x) for x in os.listdir(data_folder)]

# Folder when the parameters will be saved
folder = "./gm_country_parameters_"+data_type+"/"+gender+"/"

if not os.path.exists(folder):
    os.makedirs(folder)
    
# studied years
year_list = range(1948, 2016)
# year_list = range(1948, 1953)

# age at which the data used to fit the model stops
age_stop = 65

# age at which the data used to evaluate the fit stops
age_stop_test = 65

# function computing the contribution of aging to the hazard rate
@np.vectorize
def aging_gompertz_makaham(x, a, b, c):
    return c + a*np.exp((b*x))

# function computing the contribution of learning to the hazard rate
@np.vectorize
def learning(x, Lmax, k_learning, n):
    return (Lmax/(1 + np.exp(n*(x-k_learning)))) - Lmax

# function computing the "amount of learning" during the studied time period
def learning_amount(age_stop, Lmax, k_learning, n):
    x = np.linspace(0, age_stop, 1000)
    return integrate.trapezoid(learning(x, Lmax, k_learning, n), x)

# function computing the contribution of growth to the hazard rate
@np.vectorize
def growth(x, Gmax, growth_speed):
    results = []
    if (x < 30):
        results.append(Gmax/(1 + x**growth_speed) - Gmax ) 
    else : 
        results.append(Gmax/(1 + 30**growth_speed) - Gmax ) 
    return np.array(results)

# function summing all the contributions and returning the log of this sum. 
@np.vectorize
def log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed):
    res = aging_gompertz_makaham(x, a, b, c) + learning(x, Lmax, k_learning, n) + growth(x, Gmax, growth_speed)
    if res > 0 : 
        return np.log10(res)
    else :
        return np.log10(1.)

# cost function used for fitting
def errlearning_gompertz_makaham(params, x, y, weights, lm=0, ord=0):
    from numpy.linalg import norm
    a, b, c, Lmax, k_learning, n, Gmax, growth_speed = params

    # penalization for negative values
    if any(e<0 for e in params) == True:
        return 100

    err = norm(weights*(log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed) - y))
    return err + lm * norm(params, ord)

# cost function used for fitting Gompertz-Makeham model
def err_gompertz_makaham(params, x, y, weights, lm=0, ord=0):
    from numpy.linalg import norm
    a, b, c = params

    # penalization for negative values
    if any(e<0 for e in params) == True:
        return 100

    err = norm(weights*(np.log10(aging_gompertz_makaham(x, a, b, c)) - y))
    return err + lm * norm(params, ord)

# function used for computing the s measure (average distance between the data and the fitted model in %)
def s_measure(pred, data):
    distance = np.absolute(pred - data)
    distance_percentage = (distance/np.absolute(data))*100
    return np.mean(distance_percentage)
    # return np.average(distance_percentage, weights=weights)

def fit_gm(x, mortality_array, guess, bounds):

    # weight computation, can sometimes improve the fit but may also impact it negatively
    lx = [10000]
    for i in range(1, mortality_array.shape[0]):
        lx.append(lx[i-1]-lx[i-1]*mortality_array[i])
    lx = np.array(lx)
    surviving_fraction = lx/lx[0]
    surviving_fraction = np.square(surviving_fraction)

    # comment this line to use the computed weights
    surviving_fraction = np.ones_like(surviving_fraction)

    params = minimize(err_gompertz_makaham, guess, method="SLSQP", tol = 1e-10, args=(x, np.log10(mortality_array), surviving_fraction, 1e-4, 2), bounds = bounds)
    return params.x

def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)
        
if __name__ == "__main__":
    # Loop over all the data files
    for data_file in data_files:
        # Get country name from file name
        country = os.path.basename(data_file).split("_population")[0]
        print(country)

        # Create the csv file where the parameters for this country will be saved and write the header
        country_csv = folder+country+"_parameters.csv"

        with open(country_csv, 'a+', newline='') as file:
            if os.stat(country_csv).st_size == 0:
                header = ['year', 'a', 'b', 'c']
                append_list_as_row(country_csv, header)

        for j in range(len(year_list)):
            
            year = year_list[j]
            print(year)
            row = [year]
            
            data_path = data_file
            # read text file into pandas DataFrame
            base_dataframe = pd.read_csv(data_path, sep=" ")

            base_dataframe = base_dataframe[['Year', 'Age', gender]]
            base_dataframe = base_dataframe[base_dataframe['Year'] >= year]
            base_dataframe = base_dataframe[base_dataframe['Year'] <= year]
            base_dataframe = base_dataframe[~(base_dataframe['Age'] == "110+")]
            base_dataframe = base_dataframe.astype({"Year" : int, "Age" : int})

            df1 = base_dataframe[base_dataframe['Age'] <= age_stop]
            df1 = df1.astype(float)

            df1_test = base_dataframe[base_dataframe['Age'] <= age_stop_test]
            df1_test = df1_test.astype(float)

            beg = 0 #starting x 
            fin = age_stop #final x
            steps = age_stop+1
            x = np.linspace(beg, fin, steps)
            x_test = np.linspace(beg, age_stop_test, age_stop_test+1)

            dataframe = df1[df1['Year'] == year]
            dataframe = dataframe[gender]
            mortality_array = dataframe.to_numpy().squeeze()

            dataframe_test = df1_test[df1_test['Year'] == year]
            dataframe_test = dataframe_test[gender]
            mortality_array_test = dataframe_test.to_numpy().squeeze()

            bounds = ((0,1),(0,1),(0,1))
            guess = [1.1558369632812522e-05, 0.011827608531113518, 0.000009738166836723708]
            params = fit_gm(x, mortality_array, guess, bounds)
            for idx in range(params.shape[0]):
                row.append(params[idx])

            append_list_as_row(country_csv, row)
