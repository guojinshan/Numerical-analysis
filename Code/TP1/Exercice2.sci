function [x1,x2]=SolRacEquation(p,a,c)
    a = a;
    c = c;
    b = 2*p
    if b>0 then
        b = -b;
    end
    theta = b^2-4*a*c;
    if  theta>= 0  then
        x1 = (-b + sqrt(theta))/2*a;  
        x2 = (-b - sqrt(theta))/2*a;  
    else
        error("No racine réelle");
    end
endfunction


Pollution_Jagtvej.csv
Pol_JA = pd.read_csv("Pollution_Jagtvej.csv",sep=',',names= ['DateTime','NO2(myg/m3)', 'NOx(myg/m3)'], header=0, index_col='DateTime', parse_dates=['DateTime'])

# Inspecting data
print_data_info(Pol_JA )

# Complete the missing value
Pol_JA  = fill_missing_plot(Pol_JA , method='ffill', show=False)

# Visualize data
data_plot(Pol_JA , type_name='Pollution')

# Display the statistical characteristics
Pol_JA_resample = statistic_by_column_plot(Pol_JA,type_name='Pollution',frequency='H', list_indicator=['mean', 'min', 'max', 'std'])

# Compare statistical characteristics of each attribute 
Pol_JA_compare =  statistic_plot_comparasion(Pol_JA, type_name='Pollution', frequency='H', list_indicator=['mean', 'min', 'max', 'std'])

================================

Pollution_Ørsted Institutet.csv
Pol_ØI = pd.read_csv("Pollution_Ørsted Institutet.csv",sep=',',names= ['DateTime','NO2(myg/m3)', 'NOx(myg/m3)','O3(myg/m3)', 'CO(mg/m3)'], header=0, parse_dates=['DateTime'], index_col='DateTime')

# Inspecting data
print_data_info(Pol_ØI)

# Complete the missing value
Pol_ØI = fill_missing_plot(Pol_ØI , method='ffill', show=False)

# Visualize data
data_plot(Pol_ØI, type_name='Pollution')

# Display the statistical characteristics
Pol_ØI_resample = statistic_by_column_plot(Pol_ØI,type_name='Pollution',frequency='H', list_indicator=['mean', 'min', 'max', 'std'])

# Compare statistical characteristics of each attribute 
Pol_ØI_compare =  statistic_plot_comparasion(Pol_ØI, type_name='Pollution', frequency='H', list_indicator=['mean', 'min', 'max', 'std'])