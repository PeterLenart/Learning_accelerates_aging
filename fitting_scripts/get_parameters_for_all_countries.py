import pandas as pd
import numpy as np
from scipy.optimize import minimize
import os
import cma
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
folder = "./country_parameters_"+data_type+"/"+gender+"/"

if not os.path.exists(folder):
    os.makedirs(folder)


# Studied years

# year_list = range(1948, 1953)
year_list = range(1948, 2021)
# year_list = [2011, 2021]

# age at which the data used to fit the model stops
age_stop = 65

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
def fit_cma(x, mortality_array, sigma):
    lx = [10000]
    for i in range(1, mortality_array.shape[0]):
        lx.append(lx[i-1]-lx[i-1]*mortality_array[i])
    lx = np.array(lx)
    surviving_fraction = lx/lx[0]
    surviving_fraction = np.square(surviving_fraction)

    # comment this line to use the computed weights
    surviving_fraction = np.ones_like(surviving_fraction)
    
    global f 
    def f(params):
        return errlearning_gompertz_makaham(params, x, np.log10(mortality_array), surviving_fraction)

    es = cma.CMAEvolutionStrategy([4.76614499e-03, 2.29715390e-02, 1.87466206e-02, 5.58711082e-03, 3.63829644e+01, 1.38484469e-01, 4.34710796e-02, 0.1], sigma)
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
    # Loop over all the data files
    for data_file in data_files:

        # Get country name from file name
        country = os.path.basename(data_file).split("_population")[0]
        print(country)

        # Create the csv file where the parameters for this country will be saved and write the header
        country_csv = folder+country+"_parameters.csv"
        file = open(country_csv, 'a+', newline='')
        file.close()
        header = ['year', 'a', 'b', 'c', 'Lmax', 'k_learning', 'n', 'Gmax', 'growth_speed']
        append_list_as_row(country_csv, header)

        # Loop over the studied years
        for j in range(len(year_list)):
            year = year_list[j]
            print(year)

            # Initialization of the new row of the csv file
            row = [year]

            data_path = data_file

            # read text file into pandas DataFrame
            base_dataframe = pd.read_csv(data_path, sep=" ")

            # Only keep the Year, Age and desired gender columns
            base_dataframe = base_dataframe[['Year', 'Age', gender]]
            # Keep the row corresponding to the current year
            base_dataframe = base_dataframe[base_dataframe['Year'] == year]

            # Check if there is data for this year, pass if there is none
            if base_dataframe.empty:
                pass
            else:
                # Remove the row about mortality for age 110 or more
                base_dataframe = base_dataframe[~(base_dataframe['Age'] == "110+")]
                # Convert columns from string to int
                base_dataframe = base_dataframe.astype({"Year" : int, "Age" : int})

                 # Remove the rows for ages superior to the age where fitting stops
                df1 = base_dataframe[base_dataframe['Age'] <= age_stop]
                # Convert dataframe to float
                df1 = df1.astype(float)

                beg = 0 #starting x 
                fin = age_stop #final x
                steps = age_stop+1
                
                x = np.linspace(beg, fin, steps)

                # Remove the year and gender column
                dataframe = df1[df1['Year'] == year]
                dataframe = dataframe[gender]
                # Convert to numpy array
                mortality_array = dataframe.to_numpy().squeeze()

                # Fit the model
                es = fit_cma(x, mortality_array, 0.005)
                # Retrieve and save the fitted parameters
                params = es.result.xfavorite
                for idx in range(params.shape[0]):
                    row.append(params[idx])

                append_list_as_row(country_csv, row)