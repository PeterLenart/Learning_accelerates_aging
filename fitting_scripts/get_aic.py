import pandas as pd
import numpy as np
import os
import csv

# Male, Female or Total, studied gender of the mortality data
gender = "Total"

# cohort or population, type of study
data_type = "population"

# folders where the parameters for the GLA and Gompertz-Makeham model are saved
gla_params_folder = "./country_parameters_"+data_type+"/"
gm_params_folder = "./country_parameters_"+data_type+"/"

# Change to the folder when the mortality data is saved
data_folder = "./mortality_data/"

# List all the files in this folder
data_files = [os.path.join(data_folder, x) for x in os.listdir(data_folder)]
country_list = []
for data_file in data_files:
    # Get country name from file name
    country = os.path.basename(data_file).split("_population")[0]
    country_list.append(country)
                        
output_folder = "./country_parameters_"+data_type+"/"

# studied year
# year_list = range(1948, 2016)
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

def rss(params, x, y):
    from numpy.linalg import norm
    a, b, c, Lmax, k_learning, n, Gmax, growth_speed = params

    rss = norm((log_mortality_gompertz_makaham(x, a, b, c, Lmax, k_learning, n, Gmax, growth_speed) - y))
    return rss

def rss_gm(params, x, y):
    from numpy.linalg import norm
    a, b, c = params

    rss = norm((np.log10(aging_gompertz_makaham(x, a, b, c)) - y))
    return rss

def log_likelihood(x, params, rss):
    n = x.shape[0]
    k = params.shape[0]

    ll = -(n * 1/2) * (1 + np.log(2 * np.pi)) - (n / 2) * np.log((rss**2)/ n)

    return ll

def aic(params, x, y):
    residuals = rss(params, x, y)
    ll = log_likelihood(x, params, residuals)
    k = params.shape[0]

    AIC = (-2 * ll) + (2 * k)

    return AIC

def aic_gm(params, x, y):
    rss = rss_gm(params, x, y)
    ll = log_likelihood(x, params, rss)
    k = params.shape[0]

    AIC = (-2 * ll) + (2 * k)

    return AIC

def aicc(params, x, y):
    residuals = rss(params, x, y)
    ll = log_likelihood(x, params, residuals)
    k = params.shape[0]
    n = x.shape[0]

    AICc = (-2 * ll) + (2 * k) + (2*k*(k+1))/(n - k - 1)

    return AICc

def aicc_gm(params, x, y):
    rss = rss_gm(params, x, y)
    ll = log_likelihood(x, params, rss)
    k = params.shape[0]
    n = x.shape[0]

    AICc = (-2 * ll) + (2 * k) + (2*k*(k+1))/(n - k - 1)

    return AICc

# function used for computing the s measure (average distance between the data and the fitted model in %)
def s_measure(pred, data):
    distance = np.absolute(pred - data)
    distance_percentage = (distance/np.absolute(data))*100
    return np.mean(distance_percentage)
    # return np.average(distance_percentage, weights=weights)

def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)
        
if __name__ == "__main__":
    aic_csv = output_folder+"smeasure_summary.csv"
    file = open(aic_csv, 'a+', newline='')
    file.close()

    header = ['year'] + country_list
    append_list_as_row(aic_csv, header)

    aic_gm_csv = output_folder+"aic_summary_gm.csv"

    file = open(aic_gm_csv, 'a+', newline='')
    file.close()

    header = ['year'] + country_list
    append_list_as_row(aic_gm_csv, header)
        

    for year in year_list:
        print(year)
        row = [year]
        row_gm = [year]
        for country in country_list:
            params_path = gla_params_folder+country+"_parameters.csv"
            params_dataframe = pd.read_csv(params_path, sep=",")

            # make sure this is where your parameters are saved
            params_path_gm = gm_params_folder+country+"_parameters_gm.csv"
            params_dataframe_gm = pd.read_csv(params_path_gm, sep=",")

            data_path = data_folder+data_type+"/"+country+"_"+data_type+"_mortality.txt"
            # read text file into pandas DataFrame
            mortality_data_dataframe = pd.read_csv(data_path, sep=" ")
            mortality_data_dataframe = mortality_data_dataframe[['Year', 'Age', gender]]


            base_dataframe = mortality_data_dataframe[mortality_data_dataframe['Year'] == year]
            base_dataframe = base_dataframe[~(base_dataframe['Age'] == "110+")]
            base_dataframe = base_dataframe.astype({"Year" : int, "Age" : int})

            params_df = params_dataframe[params_dataframe['year'] == year]
            params_df = params_df.loc[:, params_df.columns != 'year']
            params = params_df.to_numpy().squeeze()

            params_df_gm = params_dataframe_gm[params_dataframe_gm['year'] == year]
            params_df_gm = params_df_gm.loc[:, params_df_gm.columns != 'year']
            params_gm = params_df_gm.to_numpy().squeeze()

            df1 = base_dataframe[base_dataframe['Age'] <= age_stop]
            df1 = df1.astype(float)

            beg = 0 #starting x 
            fin = age_stop #final x
            steps = age_stop+1
            x = np.linspace(beg, fin, steps)

            dataframe = df1[df1['Year'] == year]
            dataframe = dataframe[gender]
            mortality_array = dataframe.to_numpy().squeeze()

            row.append(aic(params, x, np.log10(mortality_array)))
            row_gm.append(aic_gm(params_gm, x, np.log10(mortality_array)))
        append_list_as_row(aic_csv, row)
        append_list_as_row(aic_gm_csv, row_gm)