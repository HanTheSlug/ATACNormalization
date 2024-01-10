import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
import warnings
import time
warnings.filterwarnings("ignore")

def removeOutliers(raw_tn5_count_matrix, pc_matrix):
    """
    Removes top and bottom 10% of outliers from raw TN5 Count DataFrame.

    Parameters:
    - raw_tn5_count_matrix (pd.DataFrame): A DataFrame containing raw TN5 count data.
    - pc_matrix (pd.DataFrame): Another DataFrame to merge with after removing outliers.

    Returns:
    - outlier_removed_raw_tn5_count_matrix (pd.DataFrame): A DataFrame with outliers removed.
    """
    
    # Drop "header" row if present
    if "header" in raw_tn5_count_matrix.index:
        raw_tn5_count_matrix = raw_tn5_count_matrix.drop("header")
    # Initialize a DataFrame to store results
    outlier_removed_raw_tn5_count_matrix = pd.DataFrame(index=raw_tn5_count_matrix.index)

    # Iterate over each column (peak) and remove the top and bottom 10% of samples
    for peak in raw_tn5_count_matrix.columns:
        values = raw_tn5_count_matrix[peak]
        lower_bound = values.quantile(0.10)
        upper_bound = values.quantile(0.90)
        filtered_values = values[(values >= lower_bound) & (values <= upper_bound)]

        # Insert the filtered values into the new DataFrame
        outlier_removed_raw_tn5_count_matrix[peak] = filtered_values

    # Merge with pc_matrix using an inner join
    outlier_removed_raw_tn5_count_matrix = outlier_removed_raw_tn5_count_matrix.join(pc_matrix, how='inner')

    return outlier_removed_raw_tn5_count_matrix

def optimizeModelParameters(outlier_removed_raw_tn5_count_matrix, pc_dataframe, alpha_values, l1_ratio_values):
    """
    Perform hyperparameter tuning to determine what the best alpha and l1 ratios are for the ElasticNet model.

    Parameters:
    - outlier_removed_raw_tn5_count_matrix (pd.DataFrame): DataFrame containing TN5 count data with outliers removed.
    - pc_dataframe (pd.DataFrame): DataFrame containing PC matrix data.
    - alpha_values (list): List of alpha hyperparameter values to try.
    - l1_ratio_values (list): List of l1_ratio hyperparameter values to try.

    Returns:
    - best_avg_mse (float): Best average MSE obtained during hyperparameter tuning.
    - best_alpha (float): Best alpha hyperparameter.
    - best_l1_ratio (float): Best l1_ratio hyperparameter.
    """
    # Split dataset into training and testing datasets. 80/20 split
    X_train, X_test, y_train, y_test = train_test_split(
        outlier_removed_raw_tn5_count_matrix[pc_dataframe.columns],
        outlier_removed_raw_tn5_count_matrix,
        test_size=0.2,
        random_state=42
    )

    # Initialize best metrics
    best_avg_mse = float('inf')
    best_alpha = None
    best_l1_ratio = None

    for alpha in alpha_values:
        for l1_ratio in l1_ratio_values:
            # Train a separate ElasticNet model for each peak
            # Model is trained on peaks because batch effects may influence each peak differently
            models = {}
            for peak in outlier_removed_raw_tn5_count_matrix.columns:

                # Identify non-NaN indices for the current peak
                non_nan_indices = y_train[peak].dropna().index

                # Check if there are non-NaN samples available
                if len(non_nan_indices) > 0:
                    # Keep only non-NaN entries
                    X_train_subset = X_train.loc[non_nan_indices]
                    y_train_subset = y_train[peak].dropna()

                    # Fit the ElasticNet model
                    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
                    model.fit(X_train_subset, y_train_subset)
                    models[peak] = model
                else:
                    # No valid samples for this peak
                    models[peak] = None
                    print(f"No valid samples for peak {peak}.")

            # Predict and evaluate MSE for each sample
            sample_mse = []
            for index, row in X_test.iterrows():

                # Observed values
                actual_values = y_test.loc[index].dropna()
                predicted_values = pd.Series()
                for peak in outlier_removed_raw_tn5_count_matrix.columns:
                    # Check if a model exists for the current peak and if that peak exists in the observed values
                    if models[peak] is not None and peak in actual_values:
                        prediction = models[peak].predict(row.values.reshape(1, -1))[0]
                        predicted_values[peak] = prediction
                    
                    # Otherwise set that peak prediction to NaN. This is due to no samples existing for that peak. 
                    else:
                        predicted_values[peak] = np.nan
                # Calculate MSE only for the peaks where valid predictions were made
                valid_predictions = predicted_values.dropna()
                valid_actual_values = actual_values[valid_predictions.index]
                if len(valid_predictions) > 0:  # Ensure there are valid predictions
                    mse = mean_squared_error(valid_actual_values, valid_predictions)
                    sample_mse.append((index, mse))
                else:
                    sample_mse.append((index, None))  # No predictions were made 

            # Filter out None values from sample_mse
            filtered_sample_mse = [mse for mse in sample_mse if mse is not None]

            # Ensure all elements in the list are numeric or np.nan
            numeric_sample_mse = []
            for mse in filtered_sample_mse:
                try:
                    numeric_sample_mse.append(float(mse[1]))
                except (ValueError, TypeError):
                    numeric_sample_mse.append(np.nan)  # Convert non-numeric values to np.nan

            # Calculate the average MSE, ignoring np.nan values
            avg_mse = np.nanmean(numeric_sample_mse) if numeric_sample_mse else float('inf')

            # Update best hyperparameters
            if avg_mse < best_avg_mse:
                best_avg_mse = avg_mse
                best_alpha = alpha
                best_l1_ratio = l1_ratio

    return best_avg_mse, best_alpha, best_l1_ratio

def normalizeTN5Counts(raw_tn5_count_matrix, pc_matrix, best_alpha, best_l1_ratio, num_pc):
    """
    Retrain ElasticNet models with the best hyperparameters and adjust TN5 count matrix entries.

    Parameters:
    - raw_tn5_count_matrix (pd.DataFrame): DataFrame containing raw TN5 count data.
    - pc_matrix (pd.DataFrame): DataFrame containing principal component matrix data.
    - best_alpha_values (list): Best alpha value.
    - best_l1_ratio_values (list): Best l1_ratio value.
    - num_pc (int): Number of PCs associated with batch effects subtracted from elasticnet model.

    Returns:
    - adjusted_tn5_counts (pd.DataFrame): DataFrame containing normalized TN5 count data
    """
    # Retrain models with the best hyperparameters
    best_models = {}
    for peak in raw_tn5_count_matrix.columns:
        non_nan_indices = raw_tn5_count_matrix[peak].dropna().index
        if len(non_nan_indices) > 0:
            X_subset = pc_matrix.loc[non_nan_indices, pc_matrix.columns[:num_pc]]
            y_subset = raw_tn5_count_matrix.loc[non_nan_indices, peak]
            model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1_ratio)
            model.fit(X_subset, y_subset)
            best_models[peak] = model
        else:
            best_models[peak] = None

    # Adjust TN5 count matrix entries
    adjusted_tn5_counts = raw_tn5_count_matrix.copy()
    for peak, model in best_models.items():
        if model is not None:
            # Extract coefficients for the first five principal components
            coeffs = model.coef_[:num_pc]

            # Multiply and subtract from the TN5 count for each sample
            for index in adjusted_tn5_counts.index:
                pc_values = pc_matrix.loc[index, pc_matrix.columns[:num_pc]]
                adjustment = np.dot(coeffs, pc_values)
                # Ensure the peak exists in the TN5 count matrix
                if peak in adjusted_tn5_counts.columns:
                    adjusted_tn5_counts.at[index, peak] -= adjustment

    # Iterate over the DataFrame to find NaN and Inf/-Inf values
    for col in adjusted_tn5_counts.columns:
        for index, value in adjusted_tn5_counts[col].items():
            if pd.isna(value):
                print(f"NaN found at Row: {index}, Column: {col}")
            elif np.isinf(value):
                print(f"Infinity found at Row: {index}, Column: {col}")

    return adjusted_tn5_counts


def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process two CSV files.")

    # Add arguments for the file paths
    parser.add_argument("--tn5_counts", type=str, required=True,
                        help="Path to a CSV containing raw tn5 counts. Format should be peaks as columns and samples as rows, with each [row,column] being a tn5 value")
    parser.add_argument("--pc", type=str, required=True,
                        help="Path to a CSV containing the PCs of all samples. Format should be PCs as columns and samples as rows")
    parser.add_argument("--out", type=str, required=True,
                        help="Path to outfile containing the normalized tn5 counts")
    parser.add_argument("--num_pc", type=int, required=True,
                        help="Number of PCs associated with batch effects subtracted from elasticnet model")

    # Parse command line arguments
    args = parser.parse_args()

    # Load raw TN5 count dataframe
    try:
        raw_tn5_count_matrix = pd.read_csv(args.tn5_counts, delimiter=',', index_col=0)
    except FileNotFoundError:
        print(f"Error: The file '{args.tn5_counts}' was not found or not in CSV format.")

    # Load the PC Dataframe
    try:
        pc_matrix = pd.read_csv(args.pc, delimiter=',', index_col=0)
    except FileNotFoundError:
        print(f"Error: The file '{args.pc}' was not found or not in CSV format.")

    start_time = time.time()

    # Call the removeOutliers function to preprocess your data
    outlier_removed_raw_tn5_count_matrix = removeOutliers(raw_tn5_count_matrix, pc_matrix)

    removeOutliers_end_time = time.time()
    elapsed_time = removeOutliers_end_time - start_time
    print(f"The removeOutliers function took {elapsed_time} seconds to run.")
    
    # Define the alpha and l1 ratios you want to test
    alpha_values = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 10]
    l1_ratio_values = [0.05, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0]
    
    # Call the optimizeModelParameters function to find the best hyperparameters
    best_avg_mse, best_alpha, best_l1_ratio = optimizeModelParameters(outlier_removed_raw_tn5_count_matrix, pc_matrix, alpha_values, l1_ratio_values)
    optimizeModelParameters_end_time = time.time()
    elapsed_time = optimizeModelParameters_end_time - removeOutliers_end_time
    print(f"The optimizeModelParameters function took {elapsed_time} seconds to run.")

    # Call the normalizeTN5Counts function to retrain models and adjust counts
    adjusted_tn5_counts = normalizeTN5Counts(raw_tn5_count_matrix, pc_matrix, best_alpha, best_l1_ratio,args.num_pc)
    normalizeTN5Counts_end_time = time.time()
    elapsed_time = normalizeTN5Counts_end_time - optimizeModelParameters_end_time
    print(f"The normalizeTN5Counts function took {elapsed_time} seconds to run.")

    # Save the adjusted TN5 count matrix to the CSV file
    adjusted_tn5_counts.to_csv(args.out, index=True, header=True)

    endtime = time.time()
    total_time = endtime - start_time
    print(f"The total program took {total_time} seconds to run.")

if __name__ == "__main__":
    main()
