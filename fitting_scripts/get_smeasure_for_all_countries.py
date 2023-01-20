import pandas as pd
import numpy as np
import os
import csv

# Male, Female or Total, studied gender of the mortality data
gender = "Total"

# cohort or population, type of study
data_type = "cohort"

# Change to the folder when the mortality data is saved
data_folder = "./mortality_data/"

# Change to the folder where the parameters are saved
params_folder = "./parameters_"+data_type+"/"+gender+"/"

# Change to the folder where you want the csv containing the smeasures to be saved
output_folder = "./country_parameters_"+data_type+"/"
s_measure_csv = output_folder+"smeasure_summary.csv"

# List all the files in this folder
data_files = [os.path.join(data_folder, x) for x in os.listdir(data_folder)]
country_list = []
for data_file in data_files:
    # Get country name from file name
    country = os.path.basename(data_file).split("_population")[0]
    country_list.append(country)

# studied years
year_list = range(1948, 1953)
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

# Utility function, appends a row at the bottom of a csv file
def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

if __name__ == "__main__":
    # Create file and write header
    file = open(s_measure_csv, 'a+', newline='')
    file.close()

    header = ['year'] + country_list
    append_list_as_row(s_measure_csv, header)

    # Loop over all the studied years
    for year in year_list:
        print(year)
        row = [year]
        # Loop over all the studied countries
        for country in country_list:
            params_path = params_folder+country+"_parameters.csv"
            params_dataframe = pd.read_csv(params_path, sep=",")

            # make sure this path leads to where your data is stored
            data_path = data_folder+data_type+"/"+country+"_"+data_type+"_mortality.txt"
            
            # Load the mortality data
            mortality_data_dataframe = pd.read_csv(data_path, sep=" ")
            # Only keep the columns that are needed
            mortality_data_dataframe = mortality_data_dataframe[['Year', 'Age', gender]]

            # Only keep rows for the current studied year
            base_dataframe = mortality_data_dataframe[mortality_data_dataframe['Year'] == year]
            # Remove the row about mortality for age 110 or more
            base_dataframe = base_dataframe[~(base_dataframe['Age'] == "110+")]
            # Convert columns from string to int
            base_dataframe = base_dataframe.astype({"Year" : int, "Age" : int})

            # Load the parameters for the studied year
            params_df = params_dataframe[params_dataframe['year'] == year]
            # Remove the year column
            params_df = params_df.loc[:, params_df.columns != 'year']
            # Convert pandas dataframe to numpy array
            params = params_df.to_numpy().squeeze()
            # Retrive individual parameters from array
            a, b, c, Lmax, k_learning, n, Gmax, growth_speed = params

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

            # Get the predicted mortality curve from the parameters
            pred = log_mortality_gompertz_makaham(x, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7])
            # Compute smeasure between thise curve and the mortality data
            measure = s_measure(pred, np.log10(mortality_array))
            # Append to the row
            row.append(measure)
        # Save row in csv
        append_list_as_row(s_measure_csv, row)